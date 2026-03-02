#!/usr/bin/env python3
"""
Iterative NGS-PrimerPlex runner with priority-based escalation.

Runs NGS_primerplex.py multiple times with escalating --primers-number1 values,
filtering successful primers into a draft file between iterations.
When priority support is enabled (column 8 in TSV), runs iterations per priority level.

Usage:
    python3 iterative_primerplex.py \
        --regions-file input.tsv \
        --reference-genome ref.fna \
        --primer-nums 10,30,100,300 \
        --run-name my_panel \
        --do-blast --threads 4 \
        [all other NGS_primerplex.py args passed through]

The script will:
  1. For each priority level (1, then 2, then 3):
     a. For each --primers-number1 value in the escalating list:
        - Run NGS_primerplex.py
        - Filter output to keep only multiplexed primers (Designed_Multiplex != '')
        - Use filtered output as draft for next iteration
        - If all regions at current priority are covered, move to next priority
  2. Report progress at each step
"""

import argparse
import subprocess
import sys
import os
import zipfile
import xml.etree.ElementTree as ET

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def filter_multiplexed_primers(info_xls, output_xls):
    """
    Filter info XLS to keep only rows where Designed_Multiplex (column 18, 0-indexed 17) is not empty.
    The filtered file can be used directly as draft input (readDraftPrimers reads info format).
    Returns (kept_count, total_count).
    """
    import shutil
    from lxml import etree

    # Read shared strings and identify rows to keep
    rows_to_keep = set()
    total_data_rows = 0
    with zipfile.ZipFile(info_xls) as z:
        with z.open('xl/sharedStrings.xml') as f:
            tree = ET.parse(f)
            ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
            shared_strings = [t.text if t.text else '' for t in tree.findall(f'.//{ns}t')]

        with z.open('xl/worksheets/sheet1.xml') as f:
            tree = ET.parse(f)
            ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
            rows = tree.findall(f'.//{ns}row')

            for row in rows:
                row_num = row.get('r')
                if row_num == '1':  # Header
                    rows_to_keep.add('1')
                    continue

                total_data_rows += 1
                cells = row.findall(f'.//{ns}c')
                # Build cell values indexed by column letter
                cell_values = {}
                for cell in cells:
                    ref = cell.get('r', '')
                    col_letter = ''.join(c for c in ref if c.isalpha())
                    v = cell.find(f'{ns}v')
                    if v is not None and v.text:
                        if cell.get('t') == 's':
                            idx = int(v.text)
                            cell_values[col_letter] = shared_strings[idx] if idx < len(shared_strings) else v.text
                        else:
                            cell_values[col_letter] = v.text
                    else:
                        cell_values[col_letter] = ''

                # Column R = Designed_Multiplex (18th column = 'R')
                designed_mult = cell_values.get('R', '')
                if designed_mult != '' and designed_mult != ' ':
                    rows_to_keep.add(row_num)

    if len(rows_to_keep) <= 1:  # Only header
        print(f'  WARNING: No multiplexed primers found in {info_xls}')
        return 0, total_data_rows

    # Copy and filter
    shutil.copy(info_xls, output_xls)

    with zipfile.ZipFile(output_xls, 'r') as zin:
        files_to_copy = {name: zin.read(name) for name in zin.namelist()}

    ns_map = {'main': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    sheet_xml = files_to_copy['xl/worksheets/sheet1.xml']
    root = etree.fromstring(sheet_xml)
    sheet_data = root.find('.//main:sheetData', ns_map)

    for row in sheet_data.findall('main:row', ns_map):
        row_num = row.get('r')
        if row_num not in rows_to_keep:
            sheet_data.remove(row)

    # Renumber rows
    new_row_num = 1
    for row in sheet_data.findall('main:row', ns_map):
        row.set('r', str(new_row_num))
        for cell in row.findall('main:c', ns_map):
            old_ref = cell.get('r')
            col_letter = ''.join(c for c in old_ref if c.isalpha())
            cell.set('r', f'{col_letter}{new_row_num}')
        new_row_num += 1

    files_to_copy['xl/worksheets/sheet1.xml'] = etree.tostring(root, xml_declaration=True, encoding='UTF-8')

    with zipfile.ZipFile(output_xls, 'w', zipfile.ZIP_DEFLATED) as zout:
        for name, content in files_to_copy.items():
            zout.writestr(name, content)

    kept = len(rows_to_keep) - 1  # Exclude header
    return kept, total_data_rows


def count_regions_by_priority(regions_file):
    """Count regions per priority level in the TSV file."""
    counts = {1: set(), 2: set(), 3: set()}
    with open(regions_file) as f:
        for line in f:
            line = line.strip()
            if not line or 'Chrom' in line or 'Start\tEnd' in line:
                continue
            cols = line.split('\t')
            name = cols[3]
            priority = 1
            if len(cols) > 7 and cols[7].strip():
                try:
                    priority = int(cols[7].strip())
                    if priority not in [1, 2, 3]:
                        priority = 1
                except ValueError:
                    priority = 1
            counts[priority].add(name)
    return {p: len(names) for p, names in counts.items() if names}


def create_priority_subset(regions_file, output_file, max_priority):
    """Create a subset TSV containing only regions with priority <= max_priority."""
    with open(regions_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            stripped = line.strip()
            if not stripped or 'Chrom' in stripped or 'Start\tEnd' in stripped:
                fout.write(line)
                continue
            cols = stripped.split('\t')
            priority = 1
            if len(cols) > 7 and cols[7].strip():
                try:
                    priority = int(cols[7].strip())
                except ValueError:
                    priority = 1
            if priority <= max_priority:
                fout.write(line)


def run_primerplex(regions_file, ref_genome, primer_num, run_name,
                   draft_file=None, extra_args=None):
    """Run NGS_primerplex.py with given parameters. Returns exit code."""
    cmd = [
        sys.executable,
        os.path.join(SCRIPT_DIR, 'NGS_primerplex.py'),
        '--regions-file', regions_file,
        '--reference-genome', ref_genome,
        '--primers-number1', str(primer_num),
        '--run-name', run_name,
    ]
    if draft_file:
        cmd.extend(['--draft-primers', draft_file])
    if extra_args:
        cmd.extend(extra_args)

    print(f'\n{"=" * 70}')
    print(f'Running NGS-PrimerPlex: --primers-number1 {primer_num}, --run-name {run_name}')
    if draft_file:
        print(f'  Draft file: {draft_file}')
    print(f'  Command: {" ".join(cmd)}')
    print(f'{"=" * 70}\n')

    result = subprocess.run(cmd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description='Iterative NGS-PrimerPlex with priority-based escalation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--regions-file', '-regions', required=True,
                        help='Input TSV with regions (with optional priority column 8)')
    parser.add_argument('--reference-genome', '-ref', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--primer-nums', '-nums', default='10,30,100,300,600',
                        help='Comma-separated escalating primer numbers (default: 10,30,100,300,600)')
    parser.add_argument('--run-name', '-run', default='iterative',
                        help='Base run name for output files')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (default: same as regions file)')
    parser.add_argument('--by-priority', action='store_true', default=True,
                        help='Run iterations separately per priority level (default: True)')
    parser.add_argument('--no-by-priority', action='store_false', dest='by_priority',
                        help='Run all regions together in each iteration')

    args, extra_args = parser.parse_known_args()

    primer_nums = [int(x) for x in args.primer_nums.split(',')]
    regions_file = os.path.abspath(args.regions_file)
    ref_genome = os.path.abspath(args.reference_genome)
    base_dir = args.output_dir or os.path.dirname(regions_file)

    # Count regions by priority
    priority_counts = count_regions_by_priority(regions_file)
    max_priority = max(priority_counts.keys())
    print(f'\n{"=" * 70}')
    print(f'Iterative NGS-PrimerPlex')
    print(f'{"=" * 70}')
    print(f'Regions file: {regions_file}')
    print(f'Reference genome: {ref_genome}')
    print(f'Primer number escalation: {primer_nums}')
    for p in sorted(priority_counts):
        print(f'  Priority {p}: {priority_counts[p]} regions')
    print(f'Mode: {"by-priority" if args.by_priority and max_priority > 1 else "all-together"}')
    print(f'Extra args: {" ".join(extra_args)}')
    print(f'{"=" * 70}\n')

    draft_file = None
    final_info_file = None

    if args.by_priority and max_priority > 1:
        # Run iterations per priority level
        for pri in sorted(priority_counts.keys()):
            print(f'\n{"#" * 70}')
            print(f'# PRIORITY {pri}: {priority_counts[pri]} regions')
            print(f'{"#" * 70}')

            # Create subset TSV for current priority level and below
            subset_file = os.path.join(base_dir,
                                       f'{os.path.splitext(os.path.basename(regions_file))[0]}_pri{pri}.tsv')
            create_priority_subset(regions_file, subset_file, pri)
            input_base = subset_file[:-4]

            for iter_idx, pnum in enumerate(primer_nums):
                run_name = f'{args.run_name}_pri{pri}_iter{iter_idx + 1}_pn{pnum}'

                rc = run_primerplex(
                    subset_file, ref_genome, pnum, run_name,
                    draft_file=draft_file,
                    extra_args=extra_args
                )

                if rc != 0:
                    print(f'WARNING: Iteration exited with code {rc}')

                # Find output info file
                info_file = f'{input_base}_{run_name}_primers_combination_1_info.xls'
                if not os.path.exists(info_file):
                    print(f'WARNING: Expected output not found: {info_file}')
                    print(f'  Trying next primer number...')
                    continue

                final_info_file = info_file

                # Filter to multiplexed primers and use as draft for next iteration
                new_draft = f'{input_base}_{run_name}_draft.xls'
                kept, total = filter_multiplexed_primers(info_file, new_draft)

                print(f'\n  Iteration result: {kept}/{total} amplicons placed in multiplex')

                if kept > 0:
                    draft_file = new_draft
                    print(f'  Draft for next iteration: {draft_file}')
                else:
                    print(f'  No primers placed - keeping previous draft')

                # Check if all regions covered (kept == total means all placed)
                if kept == total and total > 0:
                    print(f'  All amplicons placed! Moving to next priority level.')
                    break

    else:
        # Run all regions together
        input_base = regions_file[:-4]
        for iter_idx, pnum in enumerate(primer_nums):
            run_name = f'{args.run_name}_iter{iter_idx + 1}_pn{pnum}'

            rc = run_primerplex(
                regions_file, ref_genome, pnum, run_name,
                draft_file=draft_file,
                extra_args=extra_args
            )

            if rc != 0:
                print(f'WARNING: Iteration exited with code {rc}')

            info_file = f'{input_base}_{run_name}_primers_combination_1_info.xls'
            if not os.path.exists(info_file):
                print(f'WARNING: Expected output not found: {info_file}')
                continue

            final_info_file = info_file

            new_draft = f'{input_base}_{run_name}_draft.xls'
            kept, total = filter_multiplexed_primers(info_file, new_draft)

            print(f'\n  Iteration result: {kept}/{total} amplicons placed in multiplex')

            if kept > 0:
                draft_file = new_draft
            if kept == total and total > 0:
                print(f'  All amplicons placed!')
                break

    print(f'\n{"=" * 70}')
    print(f'Iterative NGS-PrimerPlex finished!')
    if final_info_file:
        print(f'Final output: {final_info_file}')
    if draft_file:
        print(f'Final draft: {draft_file}')
    print(f'{"=" * 70}')


if __name__ == '__main__':
    main()
