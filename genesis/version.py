"""When the code that is running was last changed.

The lab runs Genesis from a checkout on one Mac, and everyone else
reaches it over the network. So "which version am I looking at" is not
answerable from the browser: there is no release, no tag, and nothing on
the page that moves when the server is restarted on newer code. Somebody
reading a protocol has no way to tell whether it came from the build
that had the Tm bug in it.

Git knows, so ask git. Where there is no checkout - a copy taken to
another machine, say - fall back to the newest source file on disk,
which answers the same question less precisely and never lies about it.
The result says which of the two it is.
"""

from __future__ import annotations

import subprocess
from datetime import datetime, timezone
from pathlib import Path

#: Files whose change would mean a different Genesis. The database moves
#: constantly and says nothing about the code.
_SOURCE_GLOBS = ("genesis/*.py", "web/*.js", "web/*.css", "web/*.html")

_CACHE: dict = {}


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _from_git(root: Path) -> tuple[str, str] | None:
    """(ISO date, short sha) of the commit this checkout is on."""
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "log", "-1",
             "--format=%cs %h", "--date=short"],
            capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    parts = out.stdout.strip().split()
    if len(parts) != 2:
        return None
    return parts[0], parts[1]


def _dirty(root: Path) -> bool:
    """Whether the checkout has edits the commit date does not cover."""
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return False
    return out.returncode == 0 and bool(out.stdout.strip())


def _from_files(root: Path) -> str | None:
    newest = 0.0
    for pattern in _SOURCE_GLOBS:
        for path in root.glob(pattern):
            try:
                newest = max(newest, path.stat().st_mtime)
            except OSError:
                continue
    if not newest:
        return None
    return datetime.fromtimestamp(newest, timezone.utc).strftime("%Y-%m-%d")


def last_updated(refresh: bool = False) -> dict:
    """When this running copy of Genesis was last changed.

    Cached, because it shells out and the answer cannot change while the
    process runs - except for the working tree being edited underneath a
    dev server, which is why `refresh` exists.
    """
    if _CACHE and not refresh:
        return _CACHE

    root = _repo_root()
    git = _from_git(root)
    if git:
        date, sha = git
        dirty = _dirty(root)
        info = {"date": date, "commit": sha, "source": "git",
                "modified": dirty}
        if dirty:
            # Uncommitted edits are a different Genesis from the commit,
            # and saying only the commit date would overstate the case.
            newer = _from_files(root)
            if newer and newer > date:
                info["date"] = newer
    else:
        date = _from_files(root)
        info = {"date": date, "commit": "", "source": "files",
                "modified": False}

    _CACHE.clear()
    _CACHE.update(info)
    return info
