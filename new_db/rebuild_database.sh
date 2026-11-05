#!/usr/bin/env bash
# Run from the plasmid_caller repository root. No conda required.
set -Eeuo pipefail

usage() {
  cat <<'HELP'
Usage: bash rebuild_database.sh --seed-db PATH [options]

Required source files in --seed-db:
  ref_seqs/wp_v2.fna, ref_seqs/pf32_v4.faa, blast_parsing_dict.pkl

Options:
  --seed-db PATH        Curated seed database, with matching FASTAs/dictionary
  --seed-wp PATH        WP FASTA relative to seed-db (default: ref_seqs/wp_v2.fna)
  --seed-pf32 PATH      PF32 FASTA relative to seed-db (default: ref_seqs/pf32_v4.faa)
  --assemblies PATH     Default: new_db/assemblies.tsv
  --pfam32 PATH         Default: new_db/pfam32_references.tsv
  --renames PATH        Default: new_db/plasmid_name_corrections.tsv
  --no-renames          Do not apply a plasmid-name corrections table
  --tools PATH          Directory containing the four Python helpers (default: tools)
  --blast-bin PATH      BLAST binary directory; otherwise use bundled BLAST or PATH
  --python PATH         Python executable; default: uv run python
  --build PATH          New output directory; must not already exist
  --localized PATH      Copy existing localized inputs instead of downloading:
                        genomes/ncbi_dataset/data/, pfam32_proteins.faa,
                        pfam32_proteins.gp, optional pfam32_cds.fna
  --overrides PATH      Taxon override TSV: database, reference_id, genospecies
                        Only explicit reviewed overrides, not a whole old taxa TSV
  --min-wp-bp N         Drop WP references shorter than this (default: 1000)
  --exclude PATH        WP references to remove, with reasons
                        (default: new_db/excluded_references.tsv)
  --pf32-corrections PATH
                        Seed PF32 records to rename, with reasons
                        (default: new_db/pf32_name_corrections.tsv)
  --threads N           Self-test threads (default: 15)
  --skip-test           Build/check indexes without running the caller self-test
  -h, --help            Show this help

Required helpers:
  format_plasmid_headers.py, build_reference_taxa.py, sync_reference_taxa.py,
  filter_short_references.py

Conflicts are quarantined. Unresolved labels are reported. All output goes to
--build. The installed database is NOT replaced. Existing manual taxonomy edits
are backed up when available, but must be reapplied as reviewed overrides.
HELP
}
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
seed_db=; seed_wp=ref_seqs/wp_v2.fna; seed_pf32=ref_seqs/pf32_v4.faa
assemblies=new_db/assemblies.tsv; pfam32=new_db/pfam32_references.tsv
renames=new_db/plasmid_name_corrections.tsv; tools_dir=tools
blast_bin=; python_bin=; localized=; overrides=; excluded=new_db/excluded_references.tsv
pf32_corrections=new_db/pf32_name_corrections.tsv
build="new_db/rebuild_$(date +%Y%m%d_%H%M%S)"; threads=15; run_test=1; min_wp_bp=1000
while (($#)); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --skip-test) run_test=0; shift; continue ;;
    --no-renames) renames=; shift; continue ;;
  esac
  (($# >= 2)) || die "Missing value for $1"
  case "$1" in
    --seed-db) seed_db=$2 ;; --seed-wp) seed_wp=$2 ;;
    --seed-pf32) seed_pf32=$2 ;; --assemblies) assemblies=$2 ;;
    --pfam32) pfam32=$2 ;; --renames) renames=$2 ;;
    --tools) tools_dir=$2 ;; --blast-bin) blast_bin=$2 ;;
    --python) python_bin=$2 ;; --build) build=$2 ;;
    --localized) localized=$2 ;; --overrides) overrides=$2 ;;
    --threads) threads=$2 ;; --min-wp-bp) min_wp_bp=$2 ;; --exclude) excluded=$2 ;;
    --pf32-corrections) pf32_corrections=$2 ;;
    *) die "Unknown option: $1" ;;
  esac
  shift 2
done
[[ -n "$seed_db" ]] || die '--seed-db is required'
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die '--threads must be a positive integer'
[[ "$min_wp_bp" =~ ^[0-9]+$ ]] || die '--min-wp-bp must be a non-negative integer'
[[ ! -e "$build" ]] || die "Build directory already exists: $build"
if [[ -n "$python_bin" ]]; then PY=("$python_bin"); else
  command -v uv >/dev/null || die 'uv not found; use --python /path/to/python'
  PY=(uv run python)
fi
"${PY[@]}" -c 'import sys; assert sys.version_info >= (3, 9)'
[[ "$seed_wp" = /* ]] || seed_wp="$seed_db/$seed_wp"
[[ "$seed_pf32" = /* ]] || seed_pf32="$seed_db/$seed_pf32"
for file in "$assemblies" "$pfam32" "$seed_wp" "$seed_pf32" "$excluded" "$pf32_corrections" \
  "$seed_db/blast_parsing_dict.pkl" \
  "$tools_dir/format_plasmid_headers.py" "$tools_dir/build_reference_taxa.py" \
  "$tools_dir/sync_reference_taxa.py" "$tools_dir/filter_short_references.py"; do
  [[ -s "$file" ]] || die "Missing or empty file: $file"
done
[[ -z "$renames" || -s "$renames" ]] || die "Missing corrections table: $renames (or use --no-renames)"
[[ -z "$overrides" || -s "$overrides" ]] || die "Missing overrides: $overrides"
if [[ -z "$blast_bin" && -x src/plasmid_caller/vendor/blast/bin/makeblastdb ]]; then
  blast_bin="$PWD/src/plasmid_caller/vendor/blast/bin"
fi
if [[ -n "$blast_bin" ]]; then export PATH="$blast_bin:$PATH"; fi
command -v makeblastdb >/dev/null || die 'makeblastdb not found'
command -v blastdbcmd >/dev/null || die 'blastdbcmd not found'
if [[ -z "$localized" ]]; then
  command -v datasets >/dev/null || die 'NCBI datasets CLI not found'
  command -v curl >/dev/null || die 'curl not found'
else
  [[ -d "$localized/genomes/ncbi_dataset/data" ]] || die 'Invalid --localized genome directory'
fi
mkdir -p "$build"/{manifests,localized,seed/ref_seqs,helpers,previous_metadata}
build=$(cd "$build" && pwd)
export REBUILD_ROOT="$build"
exec > >(tee -a "$build/rebuild.log") 2>&1
trap 'rc=$?; printf "FAILED at line %s (exit %s). See %s/rebuild.log\n" "$LINENO" "$rc" "$build" >&2; exit "$rc"' ERR
printf 'Build directory: %s\n' "$build"
makeblastdb -version
cp "$assemblies" "$build/manifests/assemblies.tsv"
cp "$pfam32" "$build/manifests/pfam32_references.tsv"
cp "$excluded" "$build/manifests/excluded_references.tsv"
cp "$pf32_corrections" "$build/manifests/pf32_name_corrections.tsv"
cp "$seed_wp" "$build/seed/ref_seqs/wp_v2.fna"
cp "$seed_pf32" "$build/seed/ref_seqs/pf32_v4.faa"
cp "$seed_db/blast_parsing_dict.pkl" "$build/seed/"
for helper in format_plasmid_headers build_reference_taxa sync_reference_taxa filter_short_references; do
  cp "$tools_dir/$helper.py" "$build/helpers/"
done
rename_args=(); override_args=()
if [[ -n "$renames" ]]; then
  cp "$renames" "$build/manifests/plasmid_name_corrections.tsv"
  rename_args=(--renames "$build/manifests/plasmid_name_corrections.tsv")
fi
if [[ -n "$overrides" ]]; then
  cp "$overrides" "$build/manifests/taxon_overrides.tsv"
  override_args=(--overrides "$build/manifests/taxon_overrides.tsv")
fi
if [[ -f src/plasmid_caller/db/reference_taxa.tsv ]]; then
  cp src/plasmid_caller/db/reference_taxa.tsv "$build/previous_metadata/"
fi
printf '\n1. Validate source tables and write accession lists\n'
"${PY[@]}" - <<'PY'
import csv, os
from pathlib import Path
root = Path(os.environ['REBUILD_ROOT'])
def read(name, required):
    with (root/'manifests'/name).open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{name}: missing columns {required - set(reader.fieldnames or [])}')
        rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f'{name}: inconsistent number of TSV columns')
    return rows
assemblies = read('assemblies.tsv', {'species','strain','assembly_accession'})
proteins = read('pfam32_references.tsv', {'pfam32_type','source','protein_accession'})
seen = {}
for row in assemblies:
    acc = row['assembly_accession'].strip()
    if not acc:
        print('UNRESOLVED assembly:', row['species'], row['strain']); continue
    if acc in seen:
        raise ValueError(f'Duplicate assembly accession in manifest: {acc}')
    seen[acc] = row
protein_ids = sorted({r['protein_accession'].strip() for r in proteins
                      if r['source']=='GenBank' and r['protein_accession'].strip()})
if not seen or not protein_ids: raise ValueError('Empty assembly/PF32 accession list')
(root/'manifests/assembly_accessions.txt').write_text('\n'.join(sorted(seen))+'\n')
(root/'manifests/pf32_accessions.csv').write_text(','.join(protein_ids))
for name in ['wp_v2.fna','pf32_v4.faa']:
    with (root/'seed/ref_seqs'/name).open() as f:
        first = next((line.strip() for line in f if line.strip()), '')
    if not first.startswith('>'):
        raise ValueError(f'{name}: not a FASTA; check for a Git LFS pointer')
print(f'{len(seen)} assemblies; {len(protein_ids)} manifest-selected PF32 proteins')
PY
printf '\n2. Localize genomes and PF32 metadata\n'
if [[ -n "$localized" ]]; then
  cp -a "$localized/genomes" "$build/localized/"
  cp "$localized/pfam32_proteins.faa" "$localized/pfam32_proteins.gp" "$build/localized/"
  if [[ -f "$localized/pfam32_cds.fna" ]]; then cp "$localized/pfam32_cds.fna" "$build/localized/"; fi
else
  datasets download genome accession \
    --inputfile "$build/manifests/assembly_accessions.txt" \
    --include genome,gbff,protein,seq-report \
    --dehydrated \
    --no-progressbar \
    --filename "$build/genomes.dehydrated.zip"

  "${PY[@]}" -m zipfile -t "$build/genomes.dehydrated.zip"

  "${PY[@]}" -m zipfile -e \
    "$build/genomes.dehydrated.zip" \
    "$build/localized/genomes"

  datasets rehydrate \
    --directory "$build/localized/genomes" \
    --max-workers 3 \
    --no-progressbar
    
  for kind in fasta gp; do
    if [[ "$kind" = fasta ]]; then output=pfam32_proteins.faa; else output=pfam32_proteins.gp; fi
    curl --fail --show-error --silent --retry 3 --connect-timeout 30 --max-time 300 \
      --data-urlencode 'db=protein' \
      --data-urlencode "id@$build/manifests/pf32_accessions.csv" \
      --data-urlencode "rettype=$kind" --data-urlencode 'retmode=text' \
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' \
      > "$build/localized/$output"
    sleep 0.4
  done
fi
# The formatter otherwise reports missing assemblies and continues; stop here instead.
"${PY[@]}" - <<'PY'
import os, re
from pathlib import Path
r = Path(os.environ['REBUILD_ROOT'])
for acc in (r/'manifests/assembly_accessions.txt').read_text().splitlines():
    folder = r/'localized/genomes/ncbi_dataset/data'/acc
    fasta = [p for p in folder.glob('*_genomic.fna') if p.name != 'cds_from_genomic.fna']
    if len(fasta)!=1 or not fasta[0].stat().st_size:
        raise ValueError(f'{acc}: expected one nonempty genomic FASTA')
    if not (folder/'genomic.gbff').is_file() or not (folder/'sequence_report.jsonl').is_file():
        raise ValueError(f'{acc}: missing genomic.gbff or sequence_report.jsonl')
ids = set((r/'manifests/pf32_accessions.csv').read_text().split(','))
gp = (r/'localized/pfam32_proteins.gp').read_text()
missing = ids - set(re.findall(r'^VERSION\s+(\S+)', gp, re.M))
if missing: raise ValueError(f'Missing PF32 GenPept records: {sorted(missing)}')
if not (r/'localized/pfam32_proteins.faa').read_text().lstrip().startswith('>'):
    raise ValueError('PF32 response is not FASTA')
PY
printf '\n3. Format, deduplicate, and quarantine conflicts\n'
"${PY[@]}" "$build/helpers/format_plasmid_headers.py" \
  --existing-db "$build/seed" --localized "$build/localized" \
  --assemblies "$build/manifests/assemblies.tsv" \
  --pfam32 "$build/manifests/pfam32_references.tsv" "${rename_args[@]}" \
  --pf32-corrections "$build/manifests/pf32_name_corrections.tsv" \
  --label-conflicts separate --out "$build/deduplicated_refs"
printf '\n4. Rebuild taxonomy and synchronize the dictionary\n'
"${PY[@]}" "$build/helpers/build_reference_taxa.py" \
  --wp "$build/deduplicated_refs/wp.fna" --pf32 "$build/deduplicated_refs/pf32.faa" \
  --lookup "$build/deduplicated_refs/blast_parsing_dict.pkl" \
  --assemblies "$build/manifests/assemblies.tsv" \
  --audit "$build/deduplicated_refs/headers.tsv" \
  --genpept "$build/localized/pfam32_proteins.gp" "${override_args[@]}" \
  --out "$build/taxon_metadata"
"${PY[@]}" "$build/helpers/sync_reference_taxa.py" \
  --edited "$build/taxon_metadata/reference_taxa.tsv" \
  --db "$build/taxon_metadata" --out "$build/taxon_metadata_final"
printf '\n5. Stage final database and build indexes\n'
mkdir -p "$build/db"/{ref_seqs,wp,pf32}
cp "$build/deduplicated_refs/pf32.faa" "$build/db/ref_seqs/"
"${PY[@]}" "$build/helpers/filter_short_references.py" \
  --wp "$build/deduplicated_refs/wp.fna" \
  --taxa "$build/taxon_metadata_final/reference_taxa.tsv" \
  --lookup "$build/taxon_metadata_final/blast_parsing_dict.pkl" \
  --min-bp "$min_wp_bp" --exclude "$build/manifests/excluded_references.tsv" --out "$build/db"
if [[ -f src/plasmid_caller/db/__init__.py ]]; then cp src/plasmid_caller/db/__init__.py "$build/db/"; fi
makeblastdb -in "$build/db/ref_seqs/wp.fna" -dbtype nucl -parse_seqids -title wp -out "$build/db/wp/wp"
makeblastdb -in "$build/db/ref_seqs/pf32.faa" -dbtype prot -parse_seqids -title pf32 -out "$build/db/pf32/pf32"
blastdbcmd -db "$build/db/wp/wp" -info
blastdbcmd -db "$build/db/pf32/pf32" -info
if ((run_test)); then
  printf '\n6. Run caller self-test with fresh BLAST alignments\n'
  "${PY[@]}" -c 'import sys; from plasmid_caller.main import main; sys.exit(main())' \
    -i "$build/db/ref_seqs/wp.fna" -o "$build/self_test" \
    --database "$build/db" --threads "$threads"
fi
printf '\nBUILD COMPLETE: %s/db\n' "$build"
printf 'Review before installing:\n  %s/deduplicated_refs/conflicts.tsv\n  %s/deduplicated_refs/unresolved.tsv\n  %s/taxon_metadata_final/taxa_unresolved.tsv\n' "$build" "$build" "$build"
printf 'Original database unchanged. Manual taxonomy backup: %s/previous_metadata/\n' "$build"
