# PF32 additions: *B. afzelii* strains BO23, K78, Tom3107

**Date:** 2026-10-30  
**Author:** mjfos2r  
**Context:** Pre-rebuild database inventory. Adding PF32 partition protein sequences for the three *B. afzelii* strains present in `wp.fna` but absent from `pf32.faa`. PKo already had 6 entries in `pf32.faa` and was left unchanged.

---

## What was done

Fetched GenBank records for every plasmid accession belonging to BO23, K78, and Tom3107 (nucleotide accessions from `wp.fna` headers) via NCBI Entrez. Extracted all CDS features whose `product`, `gene`, or `note` fields matched the pattern `partition|ParA|plasmid.partition`.

**Filtered to 230–290 aa** to isolate true PF32 partition proteins.

### Why the length filter?

Every Borrelia plasmid carries a Type Ib partition module: a full-length **ParA ATPase** (~245–271 aa, the PF32-family protein actually used by the classifier) plus a shorter **BpaB/ParB-family** centromere-binding partner (~170–210 aa). Both annotations say "partition protein", so keyword matching alone returns 2 hits per plasmid. Length clearly separates them:

| Class | Range observed | Action |
|---|---|---|
| PF32 / ParA ATPase | 245–271 aa | **kept** |
| BpaB / ParB partner | 172–209 aa | dropped |
| Unrelated large ATPase (lp28-7) | 885 aa | dropped |
| K78/lp28-2 fragments (pseudogenes?) | 52–179 aa | dropped |

Lower bound 230 aa provides headroom below the smallest observed PF32 (245 aa for BO23/K78 lp38/lp17). Upper bound 290 aa provides headroom above the largest (271 aa for lp28-8). When multiple sequences passed (very rare, only K78/lp28-2), the one nearest to 253 aa (median PKo PF32 length) was kept.

---

## Results

**26 sequences added** to `src/plasmid_caller/db/ref_seqs/pf32.faa` (was 125, now 151).

| Strain | Replicons added |
|---|---|
| BO23 | cp26, cp32-3, cp32-4, cp32-5, cp9, lp17, lp28-3, lp28-4, lp28-7, lp28-8, lp38, lp54 |
| K78  | cp26, cp32-3, cp32-4, cp32-5, cp32-9, lp17, lp28-1, lp28-2, lp28-3, lp28-4, lp28-8, lp38, lp54 |
| Tom3107 | lp54 |

---

## Known gap: Tom3107 cp26 partition protein

Tom3107 cp26 (`CP009213.1`, 26,561 bp) has a broken annotation: the partition protein locus is split into two fragments of 172 aa and 189 aa — both below the filter threshold and neither resembling a full PF32 protein. This is an NCBI annotation artifact; the genomic sequence is full-length but the gene call is fragmented.

**Resolution options (not yet applied):**
1. `tblastn` a known cp26 ParA (e.g., B31 AAC66318.1, 253 aa) against CP009213.1 to recover the correct ORF coordinates, then extract and translate manually.
2. Omit Tom3107 from the PF32 database entirely (Tom3107 is the sparsest afzelii assembly — only chromosome + cp26 + lp54).

Tom3107/lp54 was successfully added (250 aa, `AIK19203.1`).

---

## What this does NOT affect

- `wp.fna` — unchanged; all four afzelii strains (BO23, K78, PKo, Tom3107) were already present.
- `blast_parsing_dict.pkl` — unchanged; this maps accession IDs to strain/plasmid name metadata and will need updating when the database is formally rebuilt via `rebuild_database.sh`.
- `reference_taxa.tsv` — currently absent (deleted from tracked db); these strains will need entries mapping to `genospecies = afzelii` when the taxonomy table is regenerated.
- `pfam32_references.tsv` — this manifest drives the formal PF32 rebuild. The protein accessions fetched here (`XDA*`, `AJY*`, `AIK*`) should be added to that TSV before the next full rebuild so the additions are reproducible.

---

## Reproducibility

Source script: `fetch_afzelii_pf32.py` was run from the session scratchpad. The protein accessions and nucleotide accession–protein accession pairs are recorded in the sequence headers:

```
>BO23_cp26_ParA_XDA25918.1 CP161025.1-XDA25918.1
>K78_lp54_ParA_AJY72836.1 CP009059.1-AJY72836.1
...
```

To reproduce: re-run `fetch_afzelii_pf32.py` with the same accession list and apply the 230–290 aa length filter.
