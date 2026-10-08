"""Accounts, passwords and sessions.

Design notes, since these choices are the whole point of the module:

* Passwords are never stored. Each one is put through PBKDF2-HMAC-SHA256
  with a fresh 16-byte salt and a high iteration count, and only the salt,
  the iteration count and the derived key are kept. The iteration count is
  stored per user, so it can be raised later without invalidating anyone.
* Verification uses ``hmac.compare_digest`` so a wrong password cannot be
  narrowed down by timing the response.
* Session tokens are 256 bits from ``secrets.token_urlsafe``. The database
  holds only their SHA-256, so a stolen database copy does not hand over
  live sessions.
* Failed logins are throttled per account with an increasing delay, which
  makes online guessing impractical without locking anyone out permanently.

What this module cannot do is protect a password in transit. Over plain
HTTP on a lab network, a password and a session cookie are both readable by
anything that can see the traffic. That needs TLS, not a better hash.
"""

from __future__ import annotations

import hashlib
import hmac
import re
import secrets
import sqlite3
import time
from dataclasses import dataclass

#: OWASP's 2023 floor for PBKDF2-HMAC-SHA256. Stored per user so it can rise.
PBKDF2_ITERATIONS = 600_000
SALT_BYTES = 16
TOKEN_BYTES = 32
DEFAULT_SESSION_HOURS = 12

#: Roles, least to most capable.
ROLES = ("viewer", "member", "admin")
ROLE_RANK = {role: i for i, role in enumerate(ROLES)}

SCHEMA = """
CREATE TABLE IF NOT EXISTS user (
    id            INTEGER PRIMARY KEY,
    username      TEXT UNIQUE NOT NULL COLLATE NOCASE,
    display_name  TEXT,
    role          TEXT NOT NULL DEFAULT 'member',
    salt          BLOB NOT NULL,
    password_hash BLOB NOT NULL,
    iterations    INTEGER NOT NULL,
    created_at    TEXT DEFAULT CURRENT_TIMESTAMP,
    last_login    TEXT,
    disabled      INTEGER NOT NULL DEFAULT 0,
    must_change    INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS session (
    token_hash  TEXT PRIMARY KEY,
    user_id     INTEGER NOT NULL REFERENCES user(id) ON DELETE CASCADE,
    created_at  REAL NOT NULL,
    expires_at  REAL NOT NULL,
    user_agent  TEXT
);
CREATE INDEX IF NOT EXISTS idx_session_user ON session(user_id);

CREATE TABLE IF NOT EXISTS login_attempt (
    username     TEXT PRIMARY KEY COLLATE NOCASE,
    failures     INTEGER NOT NULL DEFAULT 0,
    last_failure REAL
);
"""


class AuthError(Exception):
    """Raised for any failed authentication. The message is safe to show."""


@dataclass
class User:
    id: int
    username: str
    display_name: str
    role: str
    must_change: bool = False

    def can(self, required: str) -> bool:
        return ROLE_RANK.get(self.role, -1) >= ROLE_RANK.get(required, 99)

    def to_dict(self) -> dict:
        return {"id": self.id, "username": self.username,
                "display_name": self.display_name, "role": self.role,
                "must_change": self.must_change}


def ensure_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


# --------------------------------------------------------------------------
# Passwords
# --------------------------------------------------------------------------

def hash_password(password: str, salt: bytes | None = None,
                  iterations: int = PBKDF2_ITERATIONS) -> tuple[bytes, bytes, int]:
    salt = salt or secrets.token_bytes(SALT_BYTES)
    derived = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt,
                                  iterations)
    return salt, derived, iterations


def check_password(password: str, salt: bytes, expected: bytes,
                   iterations: int) -> bool:
    derived = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt,
                                  iterations)
    return hmac.compare_digest(derived, expected)


#: Deliberately mild. A long passphrase beats a short one full of symbols,
#: and rules that force symbols mostly produce "Password1!".
MIN_PASSWORD_LENGTH = 12


def password_problem(password: str, username: str = "") -> str | None:
    if len(password) < MIN_PASSWORD_LENGTH:
        return (f"use at least {MIN_PASSWORD_LENGTH} characters; a short "
                f"phrase you can remember is fine")
    if username and password.lower() == username.lower():
        return "the password cannot be the username"
    # Strip trailing digits and punctuation before comparing: "password1234"
    # and "Password!" are the same guess as "password".
    stem = re.sub(r"[^a-z]", "", password.lower())
    if password.lower() in _COMMON or stem in _COMMON_STEMS:
        return "that password is too close to a commonly used one"
    return None


#: Base words that remain weak however they are decorated.
_COMMON_STEMS = {
    "password", "passwd", "qwerty", "qwertyuiop", "letmein", "welcome",
    "changeme", "iloveyou", "admin", "administrator", "genesis", "cohenlab",
    "cohen", "plasmid", "cloning", "harvard", "secret", "monkey", "dragon",
    "sunshine", "princess", "football", "baseball", "abcdefgh", "asdfghjkl",
}

_COMMON = {
    "password", "password1", "password123", "123456789012", "qwertyuiop",
    "letmein12345", "administrator", "iloveyou1234", "welcome12345",
    "changeme1234", "genesis12345", "cohenlab1234",
}


# --------------------------------------------------------------------------
# Users
# --------------------------------------------------------------------------

def create_user(conn: sqlite3.Connection, username: str, password: str,
                role: str = "member", display_name: str = "",
                must_change: bool = False) -> User:
    username = username.strip()
    if not username:
        raise AuthError("a username is required")
    if role not in ROLES:
        raise AuthError(f"role must be one of {', '.join(ROLES)}")
    problem = password_problem(password, username)
    if problem:
        raise AuthError(problem)

    salt, derived, iterations = hash_password(password)
    try:
        cur = conn.execute(
            """INSERT INTO user(username, display_name, role, salt,
                    password_hash, iterations, must_change)
               VALUES (?,?,?,?,?,?,?)""",
            (username, display_name or username, role, salt, derived,
             iterations, int(must_change)))
    except sqlite3.IntegrityError:
        raise AuthError(f"user {username!r} already exists") from None
    conn.commit()
    return User(cur.lastrowid, username, display_name or username, role,
                must_change)


def set_password(conn: sqlite3.Connection, username: str, password: str) -> None:
    problem = password_problem(password, username)
    if problem:
        raise AuthError(problem)
    salt, derived, iterations = hash_password(password)
    cur = conn.execute(
        """UPDATE user SET salt=?, password_hash=?, iterations=?, must_change=0
           WHERE username=? COLLATE NOCASE""",
        (salt, derived, iterations, username))
    if cur.rowcount == 0:
        raise AuthError(f"no such user: {username}")
    # Changing a password ends every existing session for that account.
    conn.execute("""DELETE FROM session WHERE user_id IN
                    (SELECT id FROM user WHERE username=? COLLATE NOCASE)""",
                 (username,))
    conn.commit()


def set_role(conn: sqlite3.Connection, username: str, role: str) -> None:
    if role not in ROLES:
        raise AuthError(f"role must be one of {', '.join(ROLES)}")
    if role != "admin" and _is_last_admin(conn, username):
        raise AuthError("this is the only admin; promote someone else first")
    cur = conn.execute("UPDATE user SET role=? WHERE username=? COLLATE NOCASE",
                       (role, username))
    if cur.rowcount == 0:
        raise AuthError(f"no such user: {username}")
    conn.commit()


def set_disabled(conn: sqlite3.Connection, username: str, disabled: bool) -> None:
    if disabled and _is_last_admin(conn, username):
        raise AuthError("this is the only admin; promote someone else first")
    conn.execute("UPDATE user SET disabled=? WHERE username=? COLLATE NOCASE",
                 (int(disabled), username))
    if disabled:
        conn.execute("""DELETE FROM session WHERE user_id IN
                        (SELECT id FROM user WHERE username=? COLLATE NOCASE)""",
                     (username,))
    conn.commit()


def _is_last_admin(conn: sqlite3.Connection, username: str) -> bool:
    row = conn.execute("SELECT role FROM user WHERE username=? COLLATE NOCASE",
                       (username,)).fetchone()
    if row is None or row["role"] != "admin":
        return False
    n = conn.execute(
        "SELECT COUNT(*) c FROM user WHERE role='admin' AND disabled=0"
    ).fetchone()["c"]
    return n <= 1


def list_users(conn: sqlite3.Connection) -> list[dict]:
    return [dict(r) for r in conn.execute(
        """SELECT id, username, display_name, role, created_at, last_login,
                  disabled FROM user ORDER BY username""")]


def count_users(conn: sqlite3.Connection) -> int:
    try:
        return conn.execute("SELECT COUNT(*) c FROM user").fetchone()["c"]
    except sqlite3.OperationalError:
        return 0


# --------------------------------------------------------------------------
# Login throttling
# --------------------------------------------------------------------------

#: Delay after N consecutive failures, in seconds. Caps so an account cannot
#: be locked out indefinitely by someone else guessing at it.
def _lockout_seconds(failures: int) -> float:
    if failures < 3:
        return 0.0
    return min(2.0 ** (failures - 2), 300.0)


def _throttle_remaining(conn: sqlite3.Connection, username: str) -> float:
    row = conn.execute(
        "SELECT failures, last_failure FROM login_attempt WHERE username=?",
        (username,)).fetchone()
    if not row or not row["last_failure"]:
        return 0.0
    wait = _lockout_seconds(row["failures"])
    remaining = (row["last_failure"] + wait) - time.time()
    return max(remaining, 0.0)


def _record_failure(conn: sqlite3.Connection, username: str) -> None:
    conn.execute(
        """INSERT INTO login_attempt(username, failures, last_failure)
           VALUES (?, 1, ?)
           ON CONFLICT(username) DO UPDATE SET
             failures = failures + 1, last_failure = excluded.last_failure""",
        (username, time.time()))
    conn.commit()


def _clear_failures(conn: sqlite3.Connection, username: str) -> None:
    conn.execute("DELETE FROM login_attempt WHERE username=?", (username,))
    conn.commit()


# --------------------------------------------------------------------------
# Sessions
# --------------------------------------------------------------------------

def _token_hash(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


def login(conn: sqlite3.Connection, username: str, password: str,
          user_agent: str = "",
          session_hours: float = DEFAULT_SESSION_HOURS) -> tuple[str, User]:
    """Verify a password and open a session. Returns ``(token, user)``."""
    username = (username or "").strip()
    wait = _throttle_remaining(conn, username)
    if wait > 0:
        raise AuthError(f"too many failed attempts; try again in "
                        f"{int(wait) + 1}s")

    row = conn.execute(
        """SELECT id, username, display_name, role, salt, password_hash,
                  iterations, disabled, must_change
           FROM user WHERE username=? COLLATE NOCASE""", (username,)).fetchone()

    if row is None:
        # Spend comparable time on an unknown user so the response does not
        # reveal which usernames exist.
        hash_password(password, b"\x00" * SALT_BYTES, PBKDF2_ITERATIONS)
        _record_failure(conn, username)
        raise AuthError("incorrect username or password")

    if not check_password(password, row["salt"], row["password_hash"],
                          row["iterations"]):
        _record_failure(conn, username)
        raise AuthError("incorrect username or password")

    if row["disabled"]:
        raise AuthError("this account is disabled")

    _clear_failures(conn, username)
    token = secrets.token_urlsafe(TOKEN_BYTES)
    now = time.time()
    conn.execute(
        """INSERT INTO session(token_hash, user_id, created_at, expires_at,
                               user_agent)
           VALUES (?,?,?,?,?)""",
        (_token_hash(token), row["id"], now, now + session_hours * 3600,
         (user_agent or "")[:200]))
    conn.execute("UPDATE user SET last_login=CURRENT_TIMESTAMP WHERE id=?",
                 (row["id"],))
    conn.commit()

    return token, User(row["id"], row["username"], row["display_name"],
                       row["role"], bool(row["must_change"]))


def user_for_token(conn: sqlite3.Connection, token: str) -> User | None:
    if not token:
        return None
    row = conn.execute(
        """SELECT u.id, u.username, u.display_name, u.role, u.disabled,
                  u.must_change, s.expires_at
           FROM session s JOIN user u ON u.id = s.user_id
           WHERE s.token_hash = ?""", (_token_hash(token),)).fetchone()
    if row is None:
        return None
    if row["expires_at"] < time.time() or row["disabled"]:
        conn.execute("DELETE FROM session WHERE token_hash=?",
                     (_token_hash(token),))
        conn.commit()
        return None
    return User(row["id"], row["username"], row["display_name"], row["role"],
                bool(row["must_change"]))


def logout(conn: sqlite3.Connection, token: str) -> None:
    if token:
        conn.execute("DELETE FROM session WHERE token_hash=?",
                     (_token_hash(token),))
        conn.commit()


def purge_expired(conn: sqlite3.Connection) -> int:
    cur = conn.execute("DELETE FROM session WHERE expires_at < ?",
                       (time.time(),))
    conn.commit()
    return cur.rowcount
