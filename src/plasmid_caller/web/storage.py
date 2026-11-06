"""Durable storage for the genome index and called genomes, backed by a GCS bucket.

With no bucket configured, NullStorage keeps everything local-only -- today's
behaviour, and what local dev/tests use. google-cloud-storage is imported lazily
so the app runs without that dependency installed as long as no bucket is set.
"""
import logging

log = logging.getLogger(__name__)


class NullStorage:
    def download_db(self, local_path):
        pass

    def upload_db(self, local_path):
        pass

    def download_genome(self, genome_hash, local_dir):
        return local_dir.is_dir()

    def upload_genome(self, genome_hash, local_dir):
        pass


class GCSStorage:
    DB_BLOB = "index.sqlite3"

    def __init__(self, bucket_name, client=None):
        from google.cloud import storage as gcs

        self.client = client or gcs.Client()
        self.bucket = self.client.bucket(bucket_name)

    def download_db(self, local_path):
        blob = self.bucket.blob(self.DB_BLOB)
        if blob.exists():
            local_path.parent.mkdir(parents=True, exist_ok=True)
            blob.download_to_filename(str(local_path))

    def upload_db(self, local_path):
        self.bucket.blob(self.DB_BLOB).upload_from_filename(str(local_path))

    def _prefix(self, genome_hash):
        return f"genomes/{genome_hash}/output/"

    def download_genome(self, genome_hash, local_dir):
        prefix = self._prefix(genome_hash)
        found = False
        for blob in self.bucket.list_blobs(prefix=prefix):
            relative = blob.name[len(prefix):]
            if not relative:
                continue
            target = local_dir / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            blob.download_to_filename(str(target))
            found = True
        return found

    def upload_genome(self, genome_hash, local_dir):
        prefix = self._prefix(genome_hash)
        for path in local_dir.rglob("*"):
            if path.is_file():
                blob_name = prefix + path.relative_to(local_dir).as_posix()
                self.bucket.blob(blob_name).upload_from_filename(str(path))


def build_storage(bucket_name):
    if not bucket_name:
        return NullStorage()
    try:
        return GCSStorage(bucket_name)
    except Exception:
        log.exception("Could not initialise GCS storage for bucket %s; falling back to local-only", bucket_name)
        return NullStorage()
