"""SQLite index of called genomes, keyed by content hash. One row per genome ever called."""
import json
import sqlite3
import threading
import time

SCHEMA = """
CREATE TABLE IF NOT EXISTS genomes (
    hash TEXT PRIMARY KEY,
    job_id TEXT NOT NULL,
    name TEXT NOT NULL,
    contigs INTEGER NOT NULL,
    bases INTEGER NOT NULL,
    created_at REAL NOT NULL,
    provenance TEXT NOT NULL,
    visible INTEGER NOT NULL DEFAULT 0,
    visible_at REAL,
    visible_by TEXT
)
"""


def _row(cursor_row):
    if cursor_row is None:
        return None
    row = dict(cursor_row)
    row["provenance"] = json.loads(row["provenance"])
    row["visible"] = bool(row["visible"])
    return row


class GenomeIndex:
    """One sqlite connection, one lock. Reads are guarded too: sqlite3 connections
    aren't safe for concurrent use even with check_same_thread=False."""

    def __init__(self, db_path, storage):
        self.db_path = db_path
        self.storage = storage
        self.lock = threading.Lock()
        self.storage.download_db(self.db_path)
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute(SCHEMA)
        self.conn.commit()

    def lookup(self, genome_hash):
        with self.lock:
            cur = self.conn.execute("SELECT * FROM genomes WHERE hash = ?", (genome_hash,))
            return _row(cur.fetchone())

    def try_insert(self, genome_hash, *, job_id, name, contigs, bases, provenance):
        """Claim the canonical slot for this hash. True if this call won the race."""
        with self.lock:
            cur = self.conn.execute(
                "INSERT OR IGNORE INTO genomes (hash, job_id, name, contigs, bases, created_at, provenance) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (genome_hash, job_id, name, contigs, bases, time.time(), json.dumps(provenance)),
            )
            self.conn.commit()
            won = cur.rowcount == 1
            if won:
                self.storage.upload_db(self.db_path)
            return won

    def set_visible(self, genome_hash, visible, by=None):
        with self.lock:
            cur = self.conn.execute(
                "UPDATE genomes SET visible = ?, visible_at = ?, visible_by = ? WHERE hash = ?",
                (int(visible), time.time() if visible else None, by if visible else None, genome_hash),
            )
            self.conn.commit()
            changed = cur.rowcount == 1
            if changed:
                self.storage.upload_db(self.db_path)
            return changed

    def list_visible(self):
        with self.lock:
            cur = self.conn.execute("SELECT * FROM genomes WHERE visible = 1 ORDER BY created_at DESC")
            return [_row(row) for row in cur.fetchall()]
