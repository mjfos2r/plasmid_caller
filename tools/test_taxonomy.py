import argparse
import pickle
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

import build_reference_taxa as builder
import taxonomy as tx


class TaxonomyTests(unittest.TestCase):
    def test_hit_ids_and_shared_reference(self):
        lookup = {'CP123.1': {'name': 'cp26', 'strain': 'B31', 'genospecies': 'burgdorferi'}}
        for ident in ['CP123.1', 'gb|CP123.1|', 'lcl|CP123.1']:
            self.assertEqual(tx.lookup_reference(ident, lookup)['strain'], 'B31')
        self.assertEqual(tx.lookup_reference('gnl|BL_ORD_ID|0', lookup, 'CP123.1 title')['name'], 'cp26')
        taxa = {('wp', 'CP123.1'): {'genospecies': 'ambiguous', 'genospecies_candidates': 'afzelii;burgdorferi'}}
        label, candidates = tx.hit_taxonomy('wp', SimpleNamespace(hit_id='gb|CP123.1|'), taxa, lookup)
        self.assertEqual(label, 'ambiguous')
        self.assertEqual(candidates, 'afzelii;burgdorferi')
        with self.assertRaises(KeyError):
            tx.lookup_reference('CP123.2', lookup)

    def test_top_score_ties_and_no_hit_queries(self):
        def hit(contig, species, length, pid, coverage, covered):
            return dict(assembly_id='a', contig_id=contig, contig_len=length,
                        query_length=length, genospecies=species,
                        genospecies_candidates=species or '', overall_percent_identity=pid,
                        query_coverage_percent=coverage, query_covered_length=covered)
        full = pd.DataFrame([
            hit('one', 'afzelii', 6000, 99, 100, 6000),
            hit('one', 'burgdorferi', 6000, 98, 100, 6000),
            hit('two', 'afzelii', 2000, 99, 100, 2000),
            hit('two', 'burgdorferi', 2000, 99, 100, 2000),
            hit('three', None, 1000, None, None, None),
            hit('four', 'afzelii', 1000, 99, 10, 100),
        ])
        best = full.iloc[[0, 2]].copy()
        result = tx.annotate_best_hits(full, best, 'wp').set_index('contig_id')
        self.assertEqual(result.loc['one', 'genospecies'], 'afzelii')
        self.assertEqual(result.loc['two', 'genospecies'], 'ambiguous')
        self.assertEqual(result.loc['three', 'genospecies'], 'unclassified')
        self.assertEqual(result.loc['four', 'genospecies'], 'unclassified')
        joined = result.reset_index().rename(columns={'genospecies': 'genospecies_wp',
                                                     'genospecies_candidates': 'genospecies_candidates_wp'})
        final = tx.add_genospecies_calls(joined)
        with tempfile.TemporaryDirectory() as tmp:
            composition = tx.write_genospecies_composition(final, Path(tmp) / 'comp.tsv').set_index('genospecies')
            self.assertEqual(composition.loc['afzelii', 'percent_composition'], 75)
            self.assertEqual(composition.loc['ambiguous', 'percent_composition'], 25)
            self.assertEqual(composition.loc['unclassified', 'percent_composition'], 0)
            self.assertEqual(composition['percent_composition'].sum(), 100)

    def test_batch_denominator_and_pf32_only(self):
        frame = pd.DataFrame([
            dict(assembly_id='A', contig_id='same', contig_len=1000, genospecies='afzelii'),
            dict(assembly_id='B', contig_id='same', contig_len=2000, genospecies='burgdorferi'),
        ])
        with tempfile.TemporaryDirectory() as tmp:
            result = tx.write_genospecies_composition(frame, Path(tmp) / 'comp.tsv')
            self.assertEqual(list(result['percent_composition']), [100, 100])
            with self.assertRaises(ValueError):
                tx.write_genospecies_composition(pd.concat([frame, frame]), Path(tmp) / 'bad.tsv')
        pf_only = tx.add_genospecies_calls(frame.drop(columns='genospecies'))
        self.assertEqual(set(pf_only['genospecies']), {'unclassified'})

    def test_local_metadata_final_dictionary_and_duplicate_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'wp.fna').write_text('>CP1.1 Borreliella afzelii strain A plasmid cp26\nACGT\n')
            (root / 'pf.faa').write_text('>locus_A_cp26_ParA_X\nMABCDE\n')
            old = {'CP1.1': dict(ID='CP1.1', name='cp26', strain='A', length=99),
                   'stale': dict(ID='stale', name='cp26', strain='C', length=2)}
            with open(root / 'old.pkl', 'wb') as handle:
                pickle.dump(old, handle)
            (root / 'assemblies.tsv').write_text('species\tstrain\nBorrelia afzelii\tA\nBorrelia garinii\tB\n')
            (root / 'headers.tsv').write_text('file\tstatus\tformatted_id\tretained_id\tstrain\n'
                'wp.fna\tretained\tCP1.1\tCP1.1\tA\n'
                'wp.fna\tduplicate\tCP2.1\tCP1.1\tB\n')
            args = argparse.Namespace(wp=root/'wp.fna', pf32=root/'pf.faa', lookup=root/'old.pkl',
                assemblies=root/'assemblies.tsv', audit=root/'headers.tsv', genpept=None, overrides=None, out=root/'out')
            builder.build(args)
            with open(root/'out/blast_parsing_dict.pkl', 'rb') as handle:
                result = pickle.load(handle)
            self.assertEqual(set(result), {'CP1.1'})
            self.assertEqual(result['CP1.1']['length'], 4)
            self.assertEqual(result['CP1.1']['genospecies'], 'ambiguous')
            self.assertEqual(result['CP1.1']['genospecies_candidates'], ['afzelii', 'garinii'])
            taxa = tx.load_reference_taxa(root/'out')
            self.assertEqual(taxa['pf32', 'locus_A_cp26_ParA_X']['genospecies'], 'afzelii')

    def test_species_normalization(self):
        self.assertEqual(builder.genospecies('Borreliella burgdorferi B31'), 'burgdorferi')
        self.assertEqual(builder.genospecies('Borrelia burgdorferi'), 'burgdorferi')
        self.assertIsNone(builder.genospecies('Borrelia burgdorferi sensu lato'))
        self.assertIsNone(builder.genospecies('Borrelia sp. isolate'))


if __name__ == '__main__':
    unittest.main()
