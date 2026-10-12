"""Searching Addgene for a part the lab has not got.

This module only ever *searches*. It does not order, it does not add to a
cart, and it does not download sequences in bulk. When a design needs
something no plasmid in the share carries, the useful answer is "here are
the catalogue numbers you could order", and the ordering is a decision for
a person with a budget code.

Nothing found here is ever offered as a donor or a PCR template. A plasmid
on Addgene is weeks and a purchase order away from being in a freezer, and
letting it compete on cost against something already in the lab would make
every route estimate a lie. It is reported beside the plan, not inside it.

The API is Addgene's official read-only one, documented at
https://docs.developers.addgene.org. It needs a free token, which goes in
the GENESIS_ADDGENE_TOKEN environment variable or in the configuration
file. Without one this module says so and does nothing else: the website
is explicitly not to be scraped in its place.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

#: Addgene's documented API root.
BASE = "https://api.developers.addgene.org"

#: Where a search for plasmids goes. `q` is free text and uses the same
#: relevance ranking as the website's own search box.
SEARCH_PATH = "/catalog/plasmid/"

#: Results are cached on disk, because the same few parts get looked up
#: over and over while somebody tries different architectures, and there
#: is no reason to ask Addgene again for an answer that changes monthly.
CACHE_TTL_SECONDS = 7 * 24 * 60 * 60

TIMEOUT_SECONDS = 20


class AddgeneUnavailable(RuntimeError):
    """Raised when Addgene cannot be asked, with a reason worth printing."""


def token(cfg=None) -> str:
    """The API token, from the environment or the configuration file."""
    from_env = os.environ.get("GENESIS_ADDGENE_TOKEN", "").strip()
    if from_env:
        return from_env
    return (getattr(cfg, "addgene_token", "") or "").strip()


def _cache_path(cfg=None) -> str:
    base = os.path.expanduser(
        os.path.dirname(getattr(cfg, "db_path", "~/.genesis/genesis.db"))
        or "~/.genesis")
    return os.path.join(base, "addgene_cache.json")


def _load_cache(cfg=None) -> dict:
    try:
        with open(_cache_path(cfg), "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _save_cache(cache: dict, cfg=None) -> None:
    path = _cache_path(cfg)
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh)
    except OSError:
        pass          # a cache that cannot be written is not an error


def _request(path: str, params: dict, key: str) -> dict:
    url = BASE + path + "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(url, headers={
        "Authorization": f"Token {key}",
        "Accept": "application/json",
        "User-Agent": "Genesis/1.0 (lab cloning planner)",
    })
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 403):
            raise AddgeneUnavailable(
                "Addgene refused the token; check GENESIS_ADDGENE_TOKEN and "
                "that it carries the catalog scope") from exc
        if exc.code == 429:
            raise AddgeneUnavailable(
                "Addgene is rate-limiting this token; try again later") from exc
        raise AddgeneUnavailable(
            f"Addgene returned HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise AddgeneUnavailable(
            f"could not reach Addgene ({exc.reason})") from exc
    except (ValueError, TimeoutError) as exc:
        raise AddgeneUnavailable(f"Addgene sent no usable answer ({exc})") from exc


def _tidy(row: dict) -> dict:
    """The few fields somebody deciding what to order actually reads."""
    inserts = row.get("inserts") or []
    genes = [i.get("name") for i in inserts if i.get("name")]
    promoters = [i.get("promoter") for i in inserts if i.get("promoter")]
    article = row.get("article") or {}
    ident = row.get("id")
    return {
        "id": ident,
        "name": row.get("name") or "",
        "url": row.get("url") or (f"https://www.addgene.org/{ident}/"
                                  if ident else ""),
        "depositor": ", ".join(row.get("depositor") or []),
        "genes": genes,
        "promoters": [p for p in promoters if p],
        "backbone": (row.get("cloning") or {}).get("backbone") or "",
        "resistance": row.get("bacterial_resistance") or "",
        "type": row.get("plasmid_type") or "",
        "article": article.get("title") or "",
        "pmid": article.get("pubmed_id") or article.get("pmid") or "",
    }


def search(query: str, limit: int = 8, cfg=None,
           use_cache: bool = True, **filters) -> list[dict]:
    """Plasmids on Addgene matching `query`, best match first.

    `filters` passes through the catalogue's own fields - `genes`,
    `promoters`, `expression`, `plasmid_type` - for a narrower question
    than free text can ask.
    """
    query = (query or "").strip()
    if not query and not filters:
        return []
    key = token(cfg)
    if not key:
        raise AddgeneUnavailable(
            "no Addgene API token; register free at "
            "https://developers.addgene.org/ and put the token in "
            "GENESIS_ADDGENE_TOKEN")

    params = {"q": query, "page_size": max(1, min(limit, 50))}
    params.update({k: v for k, v in filters.items() if v})

    slug = json.dumps([query, sorted(params.items())], sort_keys=True)
    cache = _load_cache(cfg) if use_cache else {}
    hit = cache.get(slug)
    if hit and time.time() - hit.get("at", 0) < CACHE_TTL_SECONDS:
        return hit["results"][:limit]

    payload = _request(SEARCH_PATH, params, key)
    rows = payload.get("results")
    if rows is None:
        rows = payload.get("plasmids") or []
    results = [_tidy(r) for r in rows][:limit]

    cache[slug] = {"at": time.time(), "results": results}
    _save_cache(cache, cfg)
    return results


def search_url(query: str) -> str:
    """The website search a person can open when there is no token."""
    return ("https://www.addgene.org/search/catalog/plasmids/?q="
            + urllib.parse.quote(query))


def suggest_for_missing(missing: list[str], limit: int = 5,
                        cfg=None) -> dict:
    """What could be ordered for each part the lab has not got.

    Takes the names a scan reported as absent and asks Addgene about each.
    Returns one entry per name, carrying either the hits or the reason
    there are none, so a caller can always say something rather than
    quietly showing an empty list.
    """
    out: dict = {}
    for name in dict.fromkeys(n for n in missing if n):
        entry = {"query": name, "browse": search_url(name)}
        try:
            entry["results"] = search(name, limit=limit, cfg=cfg)
        except AddgeneUnavailable as exc:
            entry["results"] = []
            entry["unavailable"] = str(exc)
        out[name] = entry
    return out
