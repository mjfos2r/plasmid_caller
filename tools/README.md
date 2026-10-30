# Genospecies annotation for plasmid_caller 6.2.4

Adds `genospecies` and `genospecies_candidates` beside strain metadata in hit tables,
`genospecies_wp` and `genospecies_pf32` in merged summaries, a final `genospecies`
column, and `genospecies_composition.tsv`. The patch preserves the existing plasmid
scoring rules. It also preserves the versioned-filename fix, or applies it if absent.

## Install the code

Save the supplied files together in a folder named `plasmid_caller_taxa`, then run from your repository root, in the caller's Python environment:

```bash
python3 /path/to/plasmid_caller_taxa/apply_taxonomy.py --repo . --check
python3 /path/to/plasmid_caller_taxa/apply_taxonomy.py --repo .
```

The first command prints the proposed changes. The second validates every edit anchor,
backs up `main.py` as `main.py.before-taxonomy`, and adds `taxonomy.py`.
If your current code no longer matches the relevant 6.2.4 anchors, it stops without
changing files. It preserves other local edits, including your earlier filename fix.

If you run a non-editable installed copy of the package, reinstall it after applying
the source changes. The source paths in your logs indicate that your current caller
runs from `src/plasmid_caller`.

## Build reference metadata locally

Use the **exact final FASTAs used to build your BLAST databases**, including any
CD-HIT filtering. Use the matching whole-genome lookup dictionary, not the obsolete
one from the original database. Replace the illustrative paths below with your paths.

```bash
python3 /path/to/plasmid_caller_taxa/build_reference_taxa.py \
  --wp src/plasmid_caller/db/ref_seqs/wp_v3.fna \
  --pf32 /path/to/final_pf32.faa \
  --lookup src/plasmid_caller/db/blast_parsing_dict.pkl \
  --assemblies new_db/assemblies.tsv \
  --audit new_db/deduplicated_refs/headers.tsv \
  --genpept new_db/localized/pfam32_proteins.gp \
  --out taxon_metadata
```

This only reads local files. It writes:

- `blast_parsing_dict.pkl`: exactly the retained WP IDs, with corrected sequence lengths
  and added `genospecies` / `genospecies_candidates` fields; existing plasmid and strain
  labels are preserved. Missing WP lookup entries cause a clear error, not a guessed label.
- `reference_taxa.tsv`: both WP and PF32 identifiers and their taxonomy.
- `taxa_unresolved.tsv`: references with unknown or ambiguous taxonomy.

Species epithets (e.g. `afzelii`, `burgdorferi`) are used as genospecies labels.
The builder reads organisms from WP descriptions and GenPept records and matches
strains against the supplied assembly manifest. Older PF32 IDs are matched to complete
strain tokens, including IDs that start with locus tags. It does not assume that
unresolved strains belong to burgdorferi.

Supply the formatter's `headers.tsv` with `--audit` to preserve taxon associations from
removed duplicate aliases. A reference shared across genospecies becomes `ambiguous`.
Without this audit, the builder can only describe the retained representatives.
The builder never restores entries separated into your formatter's `review/` directory.

To resolve metadata manually, use `--overrides corrections.tsv`, whose columns are:

```text
database	reference_id	genospecies
wp	YOUR_ACCESSION_VERSION	afzelii
pf32	YOUR_EXACT_PROTEIN_FASTA_ID	garinii
```

Semicolons can preserve multiple possible genospecies, e.g. `afzelii;garinii`. Rerun
into a new output directory after making corrections.

Install the two generated metadata files in the DB directory your caller uses:

```bash
cp -n src/plasmid_caller/db/blast_parsing_dict.pkl \
      src/plasmid_caller/db/blast_parsing_dict.before-taxonomy.pkl
cp taxon_metadata/blast_parsing_dict.pkl taxon_metadata/reference_taxa.tsv \
   src/plasmid_caller/db/
```

No BLAST database rebuild is necessary for this metadata-only change. Rerun your usual
caller command. To reparse existing XML, use the same output directory and `--skip_blast`.

## Meaning of the composition output

```text
assembly_id  genospecies   contigs  assigned_bp  total_query_bp  percent_composition
sample       afzelii      12       900000       1000000         90.000000
sample       ambiguous    1         50000       1000000          5.000000
sample       unclassified 3         50000       1000000          5.000000
```

`percent_composition = 100 * sum(full contig lengths assigned to this label) / sum(all query contig lengths)`.

Each contig counts once. The percentage is **assembly-base composition under the
reference assignments**, not read abundance, cell abundance, or a count of BLAST hits.
It assigns a whole contig to a bin; it does not partition mosaic contigs into taxon-specific segments.

Assignments use WP hits passing the caller's existing `WP_MIN_COV_PCT` or
`WP_MIN_COV_BP` gates, ranked by coverage times identity. Equally scoring top hits
(with numerical tolerance 1e-9) with different genospecies are `ambiguous`.
Missing taxonomy is not guessed. A contig below the existing `MIN_CALL_BP` length
threshold becomes `unclassified`. No new species-specific identity cutoff is imposed.
Only hits returned by BLAST can be considered; the existing `max_target_seqs` limit
still applies.

PF32 reference genospecies is reported as hit metadata. A PF32-only match does not
assign the full contig to a genospecies in the composition report. The final taxon
assignment uses WP evidence independently of the plasmid-type decision.

Contigs without qualifying WP hits remain in the denominator. Multiple assemblies
are grouped separately even when they reuse contig IDs. In XML-only mode, the denominator
includes query records present in the XML files; it cannot recover queries absent
from all supplied XML files. Batch reparsing additionally writes
`ASSEMBLY_genospecies_composition.tsv` for each assembly.

## Checks

No extra runtime dependencies beyond the caller's existing pandas dependency.
The metadata builder and patcher use only the Python standard library.

```bash
python3 -m unittest discover -s /path/to/plasmid_caller_taxa -p 'test_taxonomy.py' -v
```

Checked with local fixtures: reference ID formats, cross-genospecies duplicate aliases,
co-best hit disagreements, no-hit queries, per-assembly denominators, and final FASTA /
lookup alignment. Both patched summary paths were also exercised with deterministic
alignment fixtures. BLAST itself was not run for these tests.
