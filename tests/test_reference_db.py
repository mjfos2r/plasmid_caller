from plasmid_caller.web import reference_db

MANIFEST_HEADER = "species\tstrain\tassembly_accession\tsource_bioproject\tnotes\n"


def write_manifest(tmp_path, rows):
    path = tmp_path / "assemblies.tsv"
    path.write_text(MANIFEST_HEADER + "".join(rows))
    return path


def test_strain_accessions_skips_rows_without_one(tmp_path):
    manifest = write_manifest(tmp_path, [
        "Borrelia burgdorferi\tB31\tGCA_000008685.2\tPRJ\t\n",
        "Borrelia garinii\tNoAcc\t\tPRJ\t\n",
    ])
    assert reference_db.strain_accessions(manifest) == {"B31": "GCA_000008685.2"}


def test_strain_accessions_missing_manifest_is_empty(tmp_path):
    assert reference_db.strain_accessions(tmp_path / "nope.tsv") == {}


def test_strain_job_links_joins_on_accession_prefix(tmp_path):
    manifest = write_manifest(tmp_path, [
        "Borrelia burgdorferi\tB31\tGCA_000008685.2\tPRJ\t\n",
        "Borrelia afzelii\tPKo\tGCA_000222115.1\tPRJ\t\n",
    ])
    jobs = [
        {"id": "job1", "name": "GCA_000008685.2_ASM868v2_genomic.fna"},
        {"id": "job2", "name": "unrelated.fna"},
    ]
    assert reference_db.strain_job_links(jobs, manifest) == {"B31": "job1"}


def test_strain_job_links_empty_manifest_links_nothing(tmp_path):
    jobs = [{"id": "job1", "name": "GCA_000008685.2_ASM868v2_genomic.fna"}]
    assert reference_db.strain_job_links(jobs, tmp_path / "nope.tsv") == {}
