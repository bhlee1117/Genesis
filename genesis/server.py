"""Internal web server for Genesis.

Built on the standard library only, so the tool runs on a lab machine with a
plain Python install and no package management. Binds to 127.0.0.1 by
default; the inventory contains unpublished lab constructs, so it should sit
behind the lab network or a VPN rather than on a public host.

API
---
GET  /api/auth-status                   whether a login is needed\nGET  /api/me                            the signed-in user\nPOST /api/login     {username,password} open a session\nPOST /api/logout                        end it\nPOST /api/password  {current,new}       change your own password\nGET  /api/users                         list accounts (admin)\nPOST /api/users     {action,...}        manage accounts (admin)\nGET  /api/stats                         inventory counts
GET  /api/owners                        people with plasmids
GET  /api/plasmids?q=&owner=&limit=     search the inventory
GET  /api/plasmid/<id>                  sequence + annotations
GET  /api/parts?type=&q=&limit=&all=    browse parts (no sequence)
GET  /api/part/<id>                     one part, with its sequence
GET  /api/part/<id>/sources             source plasmids for a part
GET  /api/enzymes                       enzyme catalogue
GET  /api/vocab                         controlled feature vocabulary
GET  /api/config                        effective configuration
POST /api/design        {spec}          run Stages 4-10, return a plan
POST /api/save          {spec, plan}    version a design
GET  /api/designs                       stored designs
POST /api/outcome       {...}           record an experimental result
POST /api/export/<kind> {plan, route_id}  protocol | pdf | primers | genbank
POST /api/import        {limit}         rescan the share
"""

from __future__ import annotations

import json
import sqlite3
import threading
import traceback
import ssl
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import auth
from .config import Config
from .construct import ConstructSpec
from .db import connect
from .design import design, save_design
from .enzymes import CATALOG
from .export import (primer_table_csv, protocol_markdown, protocol_pdf,
                     to_genbank)
from .importer import import_tree

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

_local = threading.local()


class Handler(BaseHTTPRequestHandler):
    cfg: Config = Config()
    db_path: str = ""

    server_version = "Genesis"

    # -- plumbing --------------------------------------------------------
    def log_message(self, fmt, *args):  # quieter console
        if "/api/" in str(args):
            super().log_message(fmt, *args)

    COOKIE = "genesis_session"

    #: Endpoints reachable without a session. Everything else needs one.
    PUBLIC_PATHS = {"/login", "/login.html", "/api/login", "/api/auth-status",
                    "/style.css", "/login.js"}

    #: Minimum role for the endpoints that change something.
    WRITE_ENDPOINTS = {
        "/api/save": "member", "/api/outcome": "member",
        "/api/import": "admin", "/api/users": "admin",
    }

    def cookie_token(self) -> str:
        raw = self.headers.get("Cookie")
        if not raw:
            return ""
        try:
            jar = SimpleCookie(raw)
        except Exception:
            return ""
        morsel = jar.get(self.COOKIE)
        return morsel.value if morsel else ""

    def current_user(self):
        if getattr(self, "_user_cached", False):
            return self._user
        self._user_cached = True
        self._user = auth.user_for_token(self.conn, self.cookie_token())
        return self._user

    def set_session_cookie(self, token: str, max_age: int) -> list[tuple]:
        flags = [f"{self.COOKIE}={token}", "Path=/", "HttpOnly",
                 "SameSite=Strict", f"Max-Age={max_age}"]
        if getattr(self.server, "is_tls", False):
            flags.append("Secure")
        return [("Set-Cookie", "; ".join(flags))]

    def clear_session_cookie(self) -> list[tuple]:
        return [("Set-Cookie",
                 f"{self.COOKIE}=; Path=/; HttpOnly; SameSite=Strict; Max-Age=0")]

    @property
    def conn(self) -> sqlite3.Connection:
        if not hasattr(_local, "conn"):
            _local.conn = connect(self.db_path)
        return _local.conn

    def _send(self, code: int, body: bytes, ctype: str, extra: dict | None = None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def json(self, data, code: int = 200):
        self._send(code, json.dumps(data, default=str).encode(), "application/json")

    def text(self, data: str, ctype="text/plain; charset=utf-8",
             filename: str | None = None):
        extra = {}
        if filename:
            extra["Content-Disposition"] = f'attachment; filename="{filename}"'
        self._send(200, data.encode(), ctype, extra)

    def body_json(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        if not length:
            return {}
        return json.loads(self.rfile.read(length) or b"{}")

    # -- routing ---------------------------------------------------------
    #: Set by serve(); False when no account exists yet, since a login with
    #: no accounts would simply lock everyone out.
    auth_active = False

    def guard(self, path: str) -> bool:
        """Return True if the request may proceed."""
        if not self.auth_active:
            return True
        if path in self.PUBLIC_PATHS:
            return True
        if self.current_user() is not None:
            return True
        if path.startswith("/api/"):
            self.json({"error": "not signed in", "login_required": True}, 401)
        else:
            self.send_response(302)
            self.send_header("Location", "/login")
            self.send_header("Content-Length", "0")
            self.end_headers()
        return False

    def do_GET(self):
        url = urlparse(self.path)
        path, query = url.path, parse_qs(url.query)
        try:
            if not self.guard(path):
                return
            if path.startswith("/api/"):
                return self.api_get(path, query)
            return self.static(path)
        except Exception as exc:  # pragma: no cover - surfaced to the client
            traceback.print_exc()
            return self.json({"error": str(exc)}, 500)

    def do_POST(self):
        url = urlparse(self.path)
        try:
            if not self.guard(url.path):
                return
            required = self.WRITE_ENDPOINTS.get(url.path)
            if required and self.auth_active:
                user = self.current_user()
                if user is None or not user.can(required):
                    return self.json(
                        {"error": f"this needs the {required} role; you are "
                                  f"signed in as "
                                  f"{user.role if user else 'nobody'}"}, 403)
            return self.api_post(url.path, self.body_json())
        except Exception as exc:
            traceback.print_exc()
            return self.json({"error": str(exc),
                              "trace": traceback.format_exc()[-1500:]}, 500)

    # -- static files ----------------------------------------------------
    def static(self, path: str):
        if path == "/login":
            path = "/login.html"
        rel = "index.html" if path in ("/", "") else path.lstrip("/")
        target = (WEB_DIR / rel).resolve()
        if not str(target).startswith(str(WEB_DIR.resolve())) or not target.is_file():
            return self._send(404, b"not found", "text/plain")
        ctype = {".html": "text/html; charset=utf-8",
                 ".js": "application/javascript; charset=utf-8",
                 ".css": "text/css; charset=utf-8",
                 ".svg": "image/svg+xml"}.get(target.suffix, "application/octet-stream")
        self._send(200, target.read_bytes(), ctype)

    # -- GET endpoints ---------------------------------------------------
    def api_get(self, path: str, q: dict):
        c = self.conn
        one = lambda s, *a: c.execute(s, a).fetchone()
        rows = lambda s, *a: [dict(r) for r in c.execute(s, a).fetchall()]

        if path == "/api/auth-status":
            return self.json({
                "require_auth": self.auth_active,
                "any_users": auth.count_users(c) > 0,
                "tls": bool(getattr(self.server, "is_tls", False)),
            })

        if path == "/api/me":
            user = self.current_user()
            return self.json({"user": user.to_dict() if user else None})

        if path == "/api/users":
            user = self.current_user()
            if self.auth_active and (user is None or not user.can("admin")):
                return self.json({"error": "admin role required"}, 403)
            return self.json(auth.list_users(c))

        if path == "/api/stats":
            return self.json({
                "plasmids": one("SELECT COUNT(*) n FROM plasmid")["n"],
                "parts": one("SELECT COUNT(*) n FROM part")["n"],
                "annotations": one("SELECT COUNT(*) n FROM annotation")["n"],
                "files": one("SELECT COUNT(*) n FROM source_file")["n"],
                "reads": one("SELECT COUNT(*) n FROM source_file WHERE role='read'")["n"],
                "owners": one("SELECT COUNT(DISTINCT owner) n FROM plasmid")["n"],
                "designs": one("SELECT COUNT(*) n FROM design")["n"],
                "by_type": rows("""SELECT feature_type, COUNT(*) n FROM part
                                   GROUP BY 1 ORDER BY n DESC"""),
            })

        if path == "/api/owners":
            return self.json(rows("""SELECT owner, COUNT(*) n FROM plasmid
                                     GROUP BY 1 ORDER BY n DESC"""))

        if path == "/api/plasmids":
            term = (q.get("q") or [""])[0].strip()
            owner = (q.get("owner") or [""])[0].strip()
            limit = min(int((q.get("limit") or ["100"])[0]), 500)
            sql = """SELECT p.id, p.name, p.owner, p.seq_length, p.circular,
                            p.verified, p.annotation_status, p.description,
                            (SELECT COUNT(*) FROM annotation a
                             WHERE a.plasmid_id = p.id) AS n_features
                     FROM plasmid p WHERE 1=1"""
            args: list = []
            if term:
                sql += """ AND (p.name LIKE ? COLLATE NOCASE
                             OR p.description LIKE ? COLLATE NOCASE
                             OR EXISTS (SELECT 1 FROM annotation a
                                        WHERE a.plasmid_id = p.id
                                          AND a.name LIKE ? COLLATE NOCASE))"""
                args += [f"%{term}%"] * 3
            if owner:
                sql += " AND p.owner = ?"
                args.append(owner)
            sql += " ORDER BY p.name LIMIT ?"
            args.append(limit)
            return self.json(rows(sql, *args))

        if path.startswith("/api/plasmid/"):
            pid = int(path.rsplit("/", 1)[1])
            p = one("SELECT * FROM plasmid WHERE id = ?", pid)
            if not p:
                return self.json({"error": "not found"}, 404)
            data = dict(p)
            data["annotations"] = rows(
                """SELECT * FROM annotation WHERE plasmid_id = ?
                   ORDER BY start""", pid)
            src = one("SELECT rel_path, format, role FROM source_file WHERE id = ?",
                      p["source_file_id"]) if p["source_file_id"] else None
            data["source"] = dict(src) if src else None
            data["alt_sources"] = rows(
                """SELECT s.rel_path FROM plasmid_alt_source a
                   JOIN source_file s ON s.id = a.source_file_id
                   WHERE a.plasmid_id = ?""", pid)
            return self.json(data)

        if path == "/api/parts":
            ftype = (q.get("type") or [""])[0]
            term = (q.get("q") or [""])[0].strip()
            limit = min(int((q.get("limit") or ["120"])[0]), 300)
            # Hide annotations too short to be a believable instance of their
            # type, so the picker cannot offer something the planner would
            # never choose itself (a 36 bp "promoter" is a primer site).
            # ?all=1 shows everything.
            show_all = (q.get("all") or ["0"])[0] in ("1", "true")
            # The sequence itself is deliberately left out: a browse query
            # returns hundreds of parts and the full DNA makes the response
            # hundreds of kilobytes per keystroke in a typeahead. Fetch one
            # part's sequence from /api/part/<id> when it is actually needed.
            sql = """SELECT p.id, p.feature_type, p.name, p.seq_length,
                            p.translated, p.exact_required, p.confidence,
                            p.contains, COUNT(o.id) AS n_sources
                     FROM part p
                     LEFT JOIN part_occurrence o ON o.part_id = p.id
                     WHERE 1=1"""
            args: list = []
            if ftype:
                sql += " AND p.feature_type = ?"
                args.append(ftype)
                if not show_all:
                    from .vocab import min_part_length
                    sql += " AND p.seq_length >= ?"
                    args.append(min_part_length(ftype))
            if term:
                sql += " AND p.name LIKE ? COLLATE NOCASE"
                args.append(f"%{term}%")
            sql += " GROUP BY p.id ORDER BY n_sources DESC, p.name LIMIT ?"
            args.append(limit)
            return self.json(rows(sql, *args))

        if path.startswith("/api/part/") and path.count("/") == 3:
            pid = int(path.rsplit("/", 1)[1])
            row = one("SELECT * FROM part WHERE id = ?", pid)
            if row is None:
                return self.json({"error": "not found"}, 404)
            return self.json(dict(row))

        if path.startswith("/api/part/") and path.endswith("/sources"):
            pid = int(path.split("/")[3])
            return self.json(rows(
                """SELECT pl.id, pl.name, pl.owner, pl.verified, pl.seq_length,
                          o.start, o.end, o.strand
                   FROM part_occurrence o JOIN plasmid pl ON pl.id = o.plasmid_id
                   WHERE o.part_id = ? ORDER BY pl.name""", pid))

        if path == "/api/enzymes":
            return self.json([{
                "name": e.name, "site": e.site, "end_type": e.end_type,
                "overhang": e.overhang_len, "stocked": e.name in self.cfg.stocked_enzymes,
                "star": e.star_activity, "dam": e.dam_sensitive,
                "type_iis": e.is_type_iis, "notes": e.notes,
            } for e in sorted(CATALOG.values(), key=lambda x: x.name)])

        if path == "/api/vocab":
            from .vocab import FEATURE_TYPES, PROTECTED_TYPES
            from .operations import OPERATIONS
            return self.json({
                "feature_types": FEATURE_TYPES,
                "protected_types": sorted(PROTECTED_TYPES),
                "operations": [{
                    "code": o.code, "meaning": o.meaning, "output": o.output,
                    "mandatory": o.mandatory,
                    "minutes": o.default_minutes,
                    "hands_on": o.hands_on_minutes, "notes": o.notes,
                } for o in OPERATIONS.values()],
            })

        if path == "/api/designs":
            return self.json(rows(
                """SELECT d.id, d.version, d.owner, d.created_at, c.name,
                          c.id AS construct_id
                   FROM design d JOIN construct c ON c.id = d.construct_id
                   ORDER BY d.created_at DESC LIMIT 100"""))

        if path.startswith("/api/design/"):
            did = int(path.rsplit("/", 1)[1])
            d = one("SELECT * FROM design WHERE id = ?", did)
            if not d:
                return self.json({"error": "not found"}, 404)
            out = dict(d)
            out["plan"] = json.loads(d["plan"])
            # Return the spec too, so a loaded design can be edited and saved
            # as a new version rather than being read-only in the client.
            con = one("SELECT name, spec, owner FROM construct WHERE id = ?",
                      d["construct_id"])
            if con:
                out["construct_name"] = con["name"]
                try:
                    out["spec"] = json.loads(con["spec"])
                except (json.JSONDecodeError, TypeError):
                    out["spec"] = None
            out["outcomes"] = rows(
                "SELECT * FROM outcome WHERE design_id = ? ORDER BY recorded_at", did)
            return self.json(out)

        if path == "/api/config":
            from dataclasses import asdict
            return self.json(asdict(self.cfg))

        return self.json({"error": "unknown endpoint"}, 404)

    # -- POST endpoints --------------------------------------------------
    def api_post(self, path: str, body: dict):
        c = self.conn

        if path == "/api/login":
            try:
                token, user = auth.login(
                    c, body.get("username", ""), body.get("password", ""),
                    user_agent=self.headers.get("User-Agent", ""),
                    session_hours=self.cfg.session_hours)
            except auth.AuthError as exc:
                return self.json({"error": str(exc)}, 401)
            payload = json.dumps({"user": user.to_dict()}).encode()
            return self._send(200, payload, "application/json",
                              dict(self.set_session_cookie(
                                  token, int(self.cfg.session_hours * 3600))))

        if path == "/api/logout":
            auth.logout(c, self.cookie_token())
            return self._send(200, b'{"ok":true}', "application/json",
                              dict(self.clear_session_cookie()))

        if path == "/api/password":
            user = self.current_user()
            if user is None:
                return self.json({"error": "not signed in"}, 401)
            try:
                # Re-authenticate before changing anything.
                auth.login(c, user.username, body.get("current", ""),
                           session_hours=0.01)
                auth.set_password(c, user.username, body.get("new", ""))
            except auth.AuthError as exc:
                return self.json({"error": str(exc)}, 400)
            return self._send(200, b'{"ok":true}', "application/json",
                              dict(self.clear_session_cookie()))

        if path == "/api/users":
            user = self.current_user()
            action = body.get("action", "create")
            try:
                if action == "create":
                    created = auth.create_user(
                        c, body["username"], body["password"],
                        role=body.get("role", "member"),
                        display_name=body.get("display_name", ""),
                        must_change=bool(body.get("must_change", True)))
                    return self.json({"user": created.to_dict()})
                if action == "role":
                    auth.set_role(c, body["username"], body["role"])
                elif action == "disable":
                    auth.set_disabled(c, body["username"], True)
                elif action == "enable":
                    auth.set_disabled(c, body["username"], False)
                elif action == "reset":
                    auth.set_password(c, body["username"], body["password"])
                else:
                    return self.json({"error": f"unknown action {action}"}, 400)
            except auth.AuthError as exc:
                return self.json({"error": str(exc)}, 400)
            except KeyError as exc:
                return self.json({"error": f"missing field {exc}"}, 400)
            return self.json({"ok": True})

        if path == "/api/design":
            spec_in = body.get("spec")
            if isinstance(spec_in, str):
                spec = ConstructSpec.from_yaml(spec_in)
            else:
                spec = ConstructSpec.from_dict(spec_in or {})
            plan = design(c, spec, self.cfg)
            return self.json(plan)

        if path == "/api/save":
            spec = ConstructSpec.from_dict(body["spec"]) if isinstance(
                body.get("spec"), dict) else ConstructSpec.from_yaml(body["spec"])
            did = save_design(c, spec, body["plan"], body.get("owner", ""))
            return self.json({"design_id": did})

        if path == "/api/outcome":
            c.execute(
                """INSERT INTO outcome(design_id, route_id, step_index, operation,
                        result, notes, recorded_by)
                   VALUES (?,?,?,?,?,?,?)""",
                (body.get("design_id"), body.get("route_id"), body.get("step_index"),
                 body.get("operation"), body.get("result"), body.get("notes"),
                 body.get("recorded_by", "")))
            c.commit()
            return self.json({"ok": True})

        if path.startswith("/api/export/"):
            kind = path.rsplit("/", 1)[1]
            plan = body["plan"]
            route = body.get("route_id")
            name = (plan.get("construct", {}).get("name") or "construct").replace(" ", "_")
            if kind == "protocol":
                return self.text(protocol_markdown(plan, self.cfg, route),
                                 "text/markdown; charset=utf-8", f"{name}_protocol.md")
            if kind == "primers":
                return self.text(primer_table_csv(plan, route),
                                 "text/csv; charset=utf-8", f"{name}_primers.csv")
            if kind == "pdf":
                data = protocol_pdf(plan, self.cfg, route)
                return self._send(200, data, "application/pdf",
                                  {"Content-Disposition":
                                   f'attachment; filename="{name}_protocol.pdf"'})
            if kind == "genbank":
                return self.text(to_genbank(plan), "chemical/x-genbank",
                                 f"{name}.gb")
            return self.json({"error": f"unknown export {kind}"}, 400)

        if path == "/api/import":
            limit = body.get("limit")
            stats = import_tree(c, self.cfg, limit=limit,
                                motif_scan=body.get("motif_scan", True))
            return self.json(stats.as_dict())

        if path.startswith("/api/plasmid/") and path.endswith("/annotate"):
            if not self.cfg.allow_annotation_edits:
                return self.json({"error": "annotation editing is disabled"}, 403)
            pid = int(path.split("/")[3])
            c.execute(
                """INSERT INTO annotation(plasmid_id, feature_type, name, start,
                        end, strand, evidence) VALUES (?,?,?,?,?,?,'manual')""",
                (pid, body["feature_type"], body.get("name", ""), body["start"],
                 body["end"], body.get("strand", 1)))
            c.commit()
            return self.json({"ok": True})

        return self.json({"error": "unknown endpoint"}, 404)


LOOPBACK = {"127.0.0.1", "::1", "localhost"}


def serve(cfg: Config, db_path: str | None = None) -> None:
    Handler.cfg = cfg
    Handler.db_path = str(db_path or cfg.resolved_db)
    conn = connect(Handler.db_path)
    auth.ensure_schema(conn)
    n_users = auth.count_users(conn)
    auth.purge_expired(conn)
    conn.close()

    reachable = cfg.host not in LOOPBACK

    # Refuse to put an unauthenticated inventory on the network. This is the
    # one configuration that would expose every unpublished construct in the
    # lab to anyone who can route to this machine.
    if reachable and not cfg.require_auth:
        raise SystemExit(
            f"refusing to bind {cfg.host}: require_auth is off, so anyone who "
            f"can reach this machine could read the whole inventory.\n"
            f"Set require_auth to true, or bind 127.0.0.1.")
    if reachable and cfg.require_auth and n_users == 0:
        raise SystemExit(
            f"refusing to bind {cfg.host}: no accounts exist yet, so nobody "
            f"could sign in and the login would be the only thing served.\n"
            f"Create the first account:  python3 -m genesis useradd <name> --admin")

    # Auth cannot work before anyone has an account. On loopback that would
    # lock the owner out of their own machine, so run open and say so; a
    # reachable bind in the same state was already refused above.
    effective_auth = cfg.require_auth and n_users > 0
    Handler.auth_active = effective_auth

    ssl_context = None
    if cfg.tls_cert and cfg.tls_key:
        ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ssl_context.load_cert_chain(cfg.tls_cert, cfg.tls_key)
        ssl_context.minimum_version = ssl.TLSVersion.TLSv1_2

    httpd = ThreadingHTTPServer((cfg.host, cfg.port), Handler)
    httpd.is_tls = ssl_context is not None
    if ssl_context is not None:
        httpd.socket = ssl_context.wrap_socket(httpd.socket, server_side=True)

    scheme = "https" if httpd.is_tls else "http"
    print(f"Genesis listening on {scheme}://{cfg.host}:{cfg.port}")
    print(f"  database : {Handler.db_path}")
    print(f"  DNA root : {cfg.dna_root}")
    if effective_auth:
        print(f"  accounts : {n_users}  (login required)")
    elif n_users == 0:
        print("  accounts : none yet - running without a login, local only.")
        print("             Create one to turn the login on:")
        print("               python3 -m genesis useradd <name> --admin")
    else:
        print(f"  accounts : {n_users}  (LOGIN DISABLED in config)")
    if reachable and not httpd.is_tls:
        print("  WARNING  : plain HTTP on a reachable address - passwords and")
        print("             session cookies cross the network unencrypted.")
        print("             Configure tls_cert / tls_key, or put this behind")
        print("             a reverse proxy that terminates TLS.")
    print("  Ctrl-C to stop.")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
        httpd.shutdown()
