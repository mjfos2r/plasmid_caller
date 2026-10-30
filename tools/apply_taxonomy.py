#!/usr/bin/env python3
"""Apply the genospecies extension to the 6.2.4 caller, preserving other edits.
Run from the repository root: python3 /path/to/apply_taxonomy.py --repo .
Validates all edit anchors first; backs up main.py before writing. No downloads.
"""
import argparse
import ast
import difflib
import re
from pathlib import Path

MARKER = '# Genospecies extension v1'


def replace_once(text, old, new):
    count = text.count(old)
    if count != 1:
        raise ValueError(f'Expected exactly one edit anchor, found {count}: {old[:100]!r}. '
                         'Nothing changed; adapt this patch to the current source.')
    return text.replace(old, new, 1)


def edit_function(source, name, callback):
    tree = ast.parse(source)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    lines = source.splitlines(keepends=True)
    old = ''.join(lines[node.lineno - 1:node.end_lineno])
    return ''.join(lines[:node.lineno - 1]) + callback(old) + ''.join(lines[node.end_lineno:])


def patched_source(source):
    if MARKER in source:
        raise ValueError('The genospecies extension is already installed')
    source = replace_once(source, 'from .blast_manager import BlastManager', '''# Genospecies extension v1
from .taxonomy import (
    load_reference_taxa, hit_taxonomy, lookup_reference, annotate_best_hits,
    add_genospecies_calls, write_genospecies_composition,
)
from .scoring import WP_MIN_COV_PCT, WP_MIN_COV_BP
from .blast_manager import BlastManager''')

    def edit_parser(block):
        block = replace_once(block, '    have_fasta = args.input is not None',
                             '    reference_taxa = load_reference_taxa(dbs_dir)\n\n    have_fasta = args.input is not None')
        block = replace_once(block, '        "strain",',
                             '        "strain",\n        "genospecies",\n        "genospecies_candidates",')
        block = replace_once(block, '                ref_length = aln.length', '''                ref_length = aln.length
                genospecies, genospecies_candidates = hit_taxonomy(
                    parsing_type, aln, reference_taxa, parsing_dict
                )''')
        block = replace_once(block, '                    strain, plasmid_name = get_name_from_acc(plasmid_id, parsing_dict)', '''                    reference_info = lookup_reference(
                        plasmid_id, parsing_dict, getattr(aln, "hit_def", ""),
                        getattr(aln, "accession", "")
                    )
                    strain, plasmid_name = reference_info["strain"], reference_info["name"]''')
        return replace_once(block, '                    "strain": strain,', '''                    "strain": strain,
                    "genospecies": genospecies,
                    "genospecies_candidates": genospecies_candidates,''')

    source = edit_function(source, 'parse_blast_xml', edit_parser)
    source = edit_function(source, 'parse_to_tsv', lambda block: replace_once(
        block, '    best_hits_df.to_csv(best_table_path, sep="\\t", index=False)',
        '''    best_hits_df = annotate_best_hits(
        full_hits_df, best_hits_df, parsing_type, WP_MIN_COV_PCT, WP_MIN_COV_BP
    )
    best_hits_df.to_csv(best_table_path, sep="\\t", index=False)'''))

    pattern = r'(?m)^( +)summary_df\["final_call"\] = summary_df\.apply\(choose_final_call, axis=1\)$'
    matches = list(re.finditer(pattern, source))
    if len(matches) != 2:
        raise ValueError('Expected both single-file and batch summary anchors; nothing changed')

    def add_summary(match):
        indent = match[1]
        target = ('args.output / f"{file_id}_genospecies_composition.tsv"'
                  if len(indent) == 16 else 'args.output / "genospecies_composition.tsv"')
        return (match[0] + '\n' + indent + 'summary_df = add_genospecies_calls(summary_df, MIN_CALL_BP)\n' +
                indent + 'write_genospecies_composition(summary_df, ' + target + ')')

    source = re.sub(pattern, add_summary, source)
    source = replace_once(source,
        '            combined_df = pandas.concat(all_summary_frames, ignore_index=True)',
        '''            combined_df = pandas.concat(all_summary_frames, ignore_index=True)
            write_genospecies_composition(combined_df, args.output / "genospecies_composition.tsv")''')
    # Preserve a previously applied filename fix; repair it if still on the old source.
    def fix_filename(block):
        old = '    file_id = Path(input_file).stem'
        if old in block:
            return replace_once(block, old,
                "    file_id = re.sub(r'\\.(fa|fna|fasta)(\\.gz)?$', '', Path(input_file).name, flags=re.IGNORECASE)")
        return block
    source = edit_function(source, 'get_blast_command', fix_filename)
    ast.parse(source)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', default='.')
    parser.add_argument('--check', action='store_true', help='Print the proposed diff without changing files')
    args = parser.parse_args()
    try:
        main_path = Path(args.repo) / 'src/plasmid_caller/main.py'
        original = main_path.read_text()
        updated = patched_source(original)
        module_path = Path(__file__).with_name('taxonomy.py')
        module = module_path.read_text()
        ast.parse(module)
        backup = main_path.with_name('main.py.before-taxonomy')
        target = main_path.with_name('taxonomy.py')
        if backup.exists() or target.exists():
            raise FileExistsError('Backup or taxonomy.py already exists; no files changed')
        if args.check:
            print(''.join(difflib.unified_diff(original.splitlines(True), updated.splitlines(True),
                                             fromfile=str(main_path), tofile=str(main_path))))
            return
        backup.write_text(original)
        target.write_text(module)
        main_path.write_text(updated)
        print(f'Updated {main_path}; original saved as {backup}')
        print(f'Added {target}. Install reference_taxa.tsv and the updated parsing dictionary in your DB directory.')
    except (OSError, ValueError, StopIteration, SyntaxError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
