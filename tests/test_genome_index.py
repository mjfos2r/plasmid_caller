from plasmid_caller.web.genome_index import GenomeIndex
from plasmid_caller.web.storage import NullStorage


def new_index(tmp_path):
    return GenomeIndex(tmp_path / "index.sqlite3", NullStorage())


class FakeStorage:
    def __init__(self):
        self.downloads = 0
        self.uploads = 0

    def download_db(self, local_path):
        self.downloads += 1

    def upload_db(self, local_path):
        self.uploads += 1

    def download_genome(self, genome_hash, local_dir):
        return False

    def upload_genome(self, genome_hash, local_dir):
        pass


def test_try_insert_wins_once(tmp_path):
    index = new_index(tmp_path)
    assert index.try_insert("h1", job_id="j1", name="a.fna", contigs=1, bases=10, provenance={"v": 1})
    assert not index.try_insert("h1", job_id="j2", name="b.fna", contigs=2, bases=20, provenance={"v": 2})
    row = index.lookup("h1")
    assert row["job_id"] == "j1" and row["provenance"] == {"v": 1} and row["visible"] is False


def test_lookup_missing_is_none(tmp_path):
    assert new_index(tmp_path).lookup("nope") is None


def test_set_visible_toggles_and_is_reflected_in_list(tmp_path):
    index = new_index(tmp_path)
    index.try_insert("h1", job_id="j1", name="a.fna", contigs=1, bases=10, provenance={})
    assert index.list_visible() == []
    assert index.set_visible("h1", True, by="pi@example.org")
    [row] = index.list_visible()
    assert row["hash"] == "h1" and row["visible_by"] == "pi@example.org" and row["visible_at"]
    assert index.set_visible("h1", False)
    assert index.list_visible() == []


def test_set_visible_unknown_hash_is_noop(tmp_path):
    assert not new_index(tmp_path).set_visible("nope", True)


def test_storage_is_consulted_on_load_and_on_writes_only(tmp_path):
    storage = FakeStorage()
    index = GenomeIndex(tmp_path / "index.sqlite3", storage)
    assert storage.downloads == 1 and storage.uploads == 0
    index.try_insert("h1", job_id="j1", name="a", contigs=1, bases=1, provenance={})
    assert storage.uploads == 1
    index.try_insert("h1", job_id="j2", name="b", contigs=1, bases=1, provenance={})  # ignored duplicate
    assert storage.uploads == 1
    index.set_visible("h1", True)
    assert storage.uploads == 2
