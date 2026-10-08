"""Genesis command line.

    python3 -m genesis init                  write a config file
    python3 -m genesis scan [--limit N]      index the shared DNA folder
    python3 -m genesis scan --force          re-read every file, not just changed ones
    python3 -m genesis reclassify             re-label what is already indexed
    python3 -m genesis infer-orfs             annotate unannotated payload CDSs
    python3 -m genesis stats                 what is in the inventory
    python3 -m genesis find <term>           search plasmids
    python3 -m genesis parts [--type CDS]    browse reusable parts
    python3 -m genesis design spec.yaml      plan a construct
    python3 -m genesis serve                 start the web interface
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import auth
from .config import Config
from .construct import ConstructSpec
from .db import connect
from .design import design as run_design, save_design
from .export import (primer_table_csv, protocol_markdown, protocol_pdf,
                     to_genbank)
from .importer import import_tree, reclassify
from .orfs import infer_all

CONFIG_PATH = Path("~/.genesis/config.json").expanduser()


def _conn(cfg: Config, override: str | None = None):
    return connect(override or cfg.resolved_db)


def cmd_init(args, cfg: Config):
    if args.root:
        cfg.dna_root = args.root
    cfg.save(CONFIG_PATH)
    print(f"wrote {CONFIG_PATH}")
    print(f"  dna_root : {cfg.dna_root}")
    print(f"  database : {cfg.resolved_db}")
    root = Path(cfg.dna_root)
    print(f"  root reachable: {'yes' if root.exists() else 'NO - mount the share first'}")


def cmd_scan(args, cfg: Config):
    conn = _conn(cfg, args.db)
    seen = {"n": 0}

    def progress(st):
        print(f"  {st.scanned:6d} files | {st.plasmids:5d} plasmids | "
              f"{st.reads:5d} reads | {st.duplicates:4d} duplicates",
              file=sys.stderr)

    print(f"scanning {cfg.dna_root}", file=sys.stderr)
    if args.force:
        print("  --force: re-reading every file, not just changed ones",
              file=sys.stderr)
    stats = import_tree(conn, cfg, limit=args.limit, progress=progress,
                        motif_scan=not args.no_motifs, force=args.force)
    d = stats.as_dict()
    print(json.dumps(d, indent=2))
    if d["errors"]:
        print(f"\n{len(d['errors'])} file(s) could not be read; first few:",
              file=sys.stderr)
        for e in d["errors"][:5]:
            print(f"  {e}", file=sys.stderr)


def cmd_reclassify(args, cfg: Config):
    conn = _conn(cfg, args.db)

    def progress(done, total):
        print(f"  {done}/{total} relabelled annotations applied", file=sys.stderr)

    print("re-deriving feature types from stored labels (the share is not read)",
          file=sys.stderr)
    print(json.dumps(reclassify(conn, progress=progress), indent=2))


def cmd_infer_orfs(args, cfg: Config):
    conn = _conn(cfg, args.db)

    def progress(done, total, found, skipped):
        note = f" | {skipped} already done" if skipped else ""
        print(f"  {done}/{total} plasmids | {found} CDS inferred{note}",
              file=sys.stderr)

    print("looking for payload CDSs nobody annotated", file=sys.stderr)
    print(json.dumps(infer_all(conn, limit=args.limit, refresh=args.refresh,
                               progress=progress), indent=2))


def cmd_stats(args, cfg: Config):
    conn = _conn(cfg, args.db)
    q = lambda s: conn.execute(s).fetchall()
    n = lambda s, *a: conn.execute(s, a).fetchone()[0]
    print(f"files indexed   : {n('SELECT COUNT(*) FROM source_file')}")
    print(f"  sequencing reads: "
          f"{n('SELECT COUNT(*) FROM source_file WHERE role = ?', 'read')}")
    print(f"plasmids        : {n('SELECT COUNT(*) FROM plasmid')}")
    print(f"annotations     : {n('SELECT COUNT(*) FROM annotation')}")
    print(f"reusable parts  : {n('SELECT COUNT(*) FROM part')}")
    print(f"people          : {n('SELECT COUNT(DISTINCT owner) FROM plasmid')}")
    print(f"designs saved   : {n('SELECT COUNT(*) FROM design')}")
    print("\nparts by type:")
    for r in q("""SELECT feature_type, COUNT(*) n FROM part
                  GROUP BY 1 ORDER BY n DESC"""):
        print(f"  {r['feature_type']:20s} {r['n']}")


def cmd_find(args, cfg: Config):
    conn = _conn(cfg, args.db)
    term = f"%{args.term}%"
    rows = conn.execute(
        """SELECT p.id, p.name, p.owner, p.seq_length, p.circular, p.verified,
                  (SELECT COUNT(*) FROM annotation a WHERE a.plasmid_id = p.id) nf
           FROM plasmid p
           WHERE p.name LIKE ? COLLATE NOCASE
              OR EXISTS (SELECT 1 FROM annotation a WHERE a.plasmid_id = p.id
                         AND a.name LIKE ? COLLATE NOCASE)
           ORDER BY p.name LIMIT ?""", (term, term, args.limit)).fetchall()
    if not rows:
        print("no matches")
        return
    print(f"{'id':>5}  {'name':40s} {'owner':20s} {'bp':>7} {'feat':>5}  topology")
    for r in rows:
        print(f"{r['id']:>5}  {r['name'][:40]:40s} {(r['owner'] or '')[:20]:20s} "
              f"{r['seq_length']:>7} {r['nf']:>5}  "
              f"{'circular' if r['circular'] else 'linear'}")


def cmd_parts(args, cfg: Config):
    conn = _conn(cfg, args.db)
    sql = """SELECT p.id, p.name, p.feature_type, p.seq_length,
                    COUNT(o.id) n FROM part p
             LEFT JOIN part_occurrence o ON o.part_id = p.id WHERE 1=1"""
    a: list = []
    if args.type:
        sql += " AND p.feature_type = ?"
        a.append(args.type)
    if args.term:
        sql += " AND p.name LIKE ? COLLATE NOCASE"
        a.append(f"%{args.term}%")
    sql += " GROUP BY p.id ORDER BY n DESC, p.name LIMIT ?"
    a.append(args.limit)
    rows = conn.execute(sql, a).fetchall()
    print(f"{'id':>5}  {'type':18s} {'name':34s} {'bp':>6} {'sources':>8}")
    for r in rows:
        print(f"{r['id']:>5}  {r['feature_type']:18s} {r['name'][:34]:34s} "
              f"{r['seq_length']:>6} {r['n']:>8}")


def cmd_design(args, cfg: Config):
    conn = _conn(cfg, args.db)
    text = Path(args.spec).read_text()
    spec = (ConstructSpec.from_dict(json.loads(text))
            if args.spec.endswith(".json") else ConstructSpec.from_yaml(text))

    plan = run_design(conn, spec, cfg)
    v = plan["validation"]

    if args.json:
        print(json.dumps(plan, indent=2, default=str))
        return

    print(protocol_markdown(plan, cfg, args.route))

    if args.primers:
        Path(args.primers).write_text(primer_table_csv(plan, args.route))
        print(f"\n<!-- primer table written to {args.primers} -->")
    if args.genbank:
        Path(args.genbank).write_text(to_genbank(plan))
        print(f"<!-- GenBank written to {args.genbank} -->")
    if args.pdf:
        Path(args.pdf).write_bytes(protocol_pdf(plan, cfg, args.route))
        print(f"<!-- PDF protocol written to {args.pdf} -->")
    if args.save:
        did = save_design(conn, spec, plan, cfg.design_visibility)
        print(f"<!-- saved as design {did} -->")

    if v["errors"]:
        sys.exit(2)


def _prompt_password(confirm: bool = True) -> str:
    import getpass
    while True:
        pw = getpass.getpass("password: ")
        problem = auth.password_problem(pw)
        if problem:
            print(f"  {problem}", file=sys.stderr)
            continue
        if confirm and getpass.getpass("confirm : ") != pw:
            print("  passwords did not match", file=sys.stderr)
            continue
        return pw


def cmd_setup(args, cfg: Config):
    """Get a new machine running in one command.

    Copying a prepared database is much faster than re-scanning the share,
    and it means everyone sees the same inventory.
    """
    import shutil
    target = cfg.resolved_db
    target.parent.mkdir(parents=True, exist_ok=True)

    if args.source:
        source = Path(args.source).expanduser()
        if not source.exists():
            sys.exit(f"error: no database at {source}")
        if target.exists() and not args.force:
            sys.exit(f"error: {target} already exists; pass --force to replace it")
        size_mb = source.stat().st_size / 1024 / 1024
        print(f"copying {size_mb:.0f} MB inventory from {source}")
        shutil.copy2(source, target)
        # A copied database may carry stale write-ahead logs.
        for suffix in ("-wal", "-shm"):
            stray = Path(str(target) + suffix)
            if stray.exists():
                stray.unlink()
    elif not target.exists():
        print("no database given and none exists yet; scanning the share "
              "(this takes a few minutes)", file=sys.stderr)
        conn = connect(target)
        stats = import_tree(conn, cfg)
        print(f"  indexed {stats.plasmids} plasmids")

    conn = connect(target)
    n = conn.execute("SELECT COUNT(*) FROM plasmid").fetchone()[0]
    parts = conn.execute("SELECT COUNT(*) FROM part").fetchone()[0]
    cfg.save(CONFIG_PATH)

    print(f"\nready: {n:,} plasmids, {parts:,} parts")
    print(f"  database : {target}")
    print("\nStart it with:")
    print("  python3 -m genesis serve")
    print(f"then open  http://127.0.0.1:{cfg.port}  in a browser on this "
          f"computer.")


def cmd_useradd(args, cfg: Config):
    conn = _conn(cfg, args.db)
    auth.ensure_schema(conn)
    password = _prompt_password()
    try:
        user = auth.create_user(conn, args.username, password,
                                role="admin" if args.admin else args.role,
                                display_name=args.name or "",
                                must_change=False)
    except auth.AuthError as exc:
        sys.exit(f"error: {exc}")
    print(f"created {user.username} ({user.role})")


def cmd_users(args, cfg: Config):
    conn = _conn(cfg, args.db)
    auth.ensure_schema(conn)
    rows = auth.list_users(conn)
    if not rows:
        print("no accounts yet; create one with: "
              "python3 -m genesis useradd <name> --admin")
        return
    print(f"{'username':20s} {'role':8s} {'last login':20s} status")
    for r in rows:
        print(f"{r['username'][:20]:20s} {r['role']:8s} "
              f"{(r['last_login'] or 'never')[:19]:20s} "
              f"{'disabled' if r['disabled'] else 'active'}")


def cmd_passwd(args, cfg: Config):
    conn = _conn(cfg, args.db)
    auth.ensure_schema(conn)
    password = _prompt_password()
    try:
        auth.set_password(conn, args.username, password)
    except auth.AuthError as exc:
        sys.exit(f"error: {exc}")
    print(f"password updated for {args.username}; their sessions were ended")


def cmd_role(args, cfg: Config):
    conn = _conn(cfg, args.db)
    auth.ensure_schema(conn)
    try:
        if args.disable:
            auth.set_disabled(conn, args.username, True)
            print(f"{args.username} disabled")
        elif args.enable:
            auth.set_disabled(conn, args.username, False)
            print(f"{args.username} enabled")
        else:
            auth.set_role(conn, args.username, args.role)
            print(f"{args.username} is now {args.role}")
    except auth.AuthError as exc:
        sys.exit(f"error: {exc}")


def cmd_serve(args, cfg: Config):
    from .server import serve
    if args.port:
        cfg.port = args.port
    if args.host:
        cfg.host = args.host
    serve(cfg, args.db)


def main(argv=None):
    cfg = Config.load(CONFIG_PATH)
    p = argparse.ArgumentParser(
        prog="genesis",
        description="Genesis - shared cloning designer for the lab's plasmid inventory")
    p.add_argument("--db", help="override the database path")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="write a configuration file")
    s.add_argument("--root", help="path to the shared DNA folder")
    s.set_defaults(fn=cmd_init)

    s = sub.add_parser("scan", help="index the shared DNA folder")
    s.add_argument("--limit", type=int, help="stop after N files")
    s.add_argument("--no-motifs", action="store_true",
                   help="skip de-novo motif annotation (faster)")
    s.add_argument("--force", action="store_true",
                   help="re-read every file, including ones whose timestamp "
                        "has not changed (use after a classifier change)")
    s.set_defaults(fn=cmd_scan)

    s = sub.add_parser(
        "reclassify",
        help="re-apply the labelling rules to the inventory already indexed")
    s.set_defaults(fn=cmd_reclassify)

    s = sub.add_parser(
        "infer-orfs",
        help="annotate payload CDSs that the source files never marked")
    s.add_argument("--limit", type=int, help="stop after N plasmids")
    s.add_argument("--refresh", action="store_true",
                   help="discard previously inferred CDSs and redo them")
    s.set_defaults(fn=cmd_infer_orfs)

    s = sub.add_parser("stats", help="summarise the inventory")
    s.set_defaults(fn=cmd_stats)

    s = sub.add_parser("find", help="search plasmids by name or feature")
    s.add_argument("term")
    s.add_argument("--limit", type=int, default=40)
    s.set_defaults(fn=cmd_find)

    s = sub.add_parser("parts", help="browse reusable parts")
    s.add_argument("term", nargs="?", default="")
    s.add_argument("--type", help="filter by feature type, e.g. CDS")
    s.add_argument("--limit", type=int, default=40)
    s.set_defaults(fn=cmd_parts)

    s = sub.add_parser("design", help="plan a construct from a spec file")
    s.add_argument("spec", help="YAML or JSON construct specification")
    s.add_argument("--route", help="route id to detail (default: best ranked)")
    s.add_argument("--json", action="store_true", help="emit the raw plan")
    s.add_argument("--primers", help="also write a primer CSV here")
    s.add_argument("--genbank", help="also write the construct GenBank here")
    s.add_argument("--pdf", help="also write a printable PDF protocol here")
    s.add_argument("--save", action="store_true", help="store this design")
    s.set_defaults(fn=cmd_design)

    s = sub.add_parser(
        "setup", help="prepare this machine (copy an inventory, or scan)")
    s.add_argument("--from", dest="source",
                   help="path to a prepared genesis.db to copy")
    s.add_argument("--force", action="store_true",
                   help="replace an existing database")
    s.set_defaults(fn=cmd_setup)

    s = sub.add_parser("useradd", help="create an account")
    s.add_argument("username")
    s.add_argument("--admin", action="store_true", help="make them an admin")
    s.add_argument("--role", default="member", choices=auth.ROLES)
    s.add_argument("--name", help="display name")
    s.set_defaults(fn=cmd_useradd)

    s = sub.add_parser("users", help="list accounts")
    s.set_defaults(fn=cmd_users)

    s = sub.add_parser("passwd", help="set someone's password")
    s.add_argument("username")
    s.set_defaults(fn=cmd_passwd)

    s = sub.add_parser("role", help="change a role, or disable an account")
    s.add_argument("username")
    s.add_argument("--role", default="member", choices=auth.ROLES)
    s.add_argument("--disable", action="store_true")
    s.add_argument("--enable", action="store_true")
    s.set_defaults(fn=cmd_role)

    s = sub.add_parser("serve", help="start the internal web interface")
    s.add_argument("--port", type=int)
    s.add_argument("--host")
    s.set_defaults(fn=cmd_serve)

    args = p.parse_args(argv)
    return args.fn(args, cfg)


if __name__ == "__main__":
    main()
