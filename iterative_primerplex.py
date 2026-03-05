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

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def filter_multiplexed_primers(info_tsv, output_tsv):
    """
    Filter info TSV to keep only rows where Designed_Multiplex (column 18, 0-indexed 17) is not empty.
    The filtered file can be used directly as draft input (readDraftPrimers reads info format).
    Returns (kept_count, total_count).
    """
    kept_rows = []
    header = None
    total_data_rows = 0

    # Find the actual TSV file - could be single-sheet or multi-sheet naming
    actual_file = info_tsv
    if not os.path.exists(actual_file):
        # Try multi-sheet naming: base_ngs_primerplex_internal_primers.tsv
        base = info_tsv
        for ext in ('.tsv', '.xls'):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        alt = base + '_ngs_primerplex_internal_primers.tsv'
        if os.path.exists(alt):
            actual_file = alt
        else:
            print(f'  WARNING: File not found: {info_tsv}')
            return 0, 0

    with open(actual_file) as f:
        for i, line in enumerate(f):
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if i == 0:
                header = line
                continue
            total_data_rows += 1
            # Column 17 (0-indexed) = Designed_Multiplex
            if len(cols) > 17 and cols[17].strip() not in ('', ' '):
                kept_rows.append(line)

    if not kept_rows:
        print(f'  WARNING: No multiplexed primers found in {actual_file}')
        return 0, total_data_rows

    with open(output_tsv, 'w') as f:
        f.write(header + '\n')
        for row in kept_rows:
            f.write(row + '\n')

    return len(kept_rows), total_data_rows


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


def _find_intermediate_draft(input_base):
    """Find the best intermediate draft file produced by NGS_primerplex.py.
    These files are named based on the input regions file (not --run-name):
      *_all_draft_primers_after_specificity_draft_internal.tsv  (best - passed specificity)
      *_all_draft_primers_after_SNPs_draft_internal.tsv         (passed SNP check)
      *_all_draft_primers_after_joinment_draft_internal.tsv     (after block joining)
      *_all_draft_primers_draft_internal.tsv                    (all candidates from primer3)
    Also checks legacy .xls names for backwards compatibility.
    Returns the path to the best available file, or None.
    """
    candidates = [
        f'{input_base}_all_draft_primers_after_specificity_draft_internal.tsv',
        f'{input_base}_all_draft_primers_after_SNPs_draft_internal.tsv',
        f'{input_base}_all_draft_primers_after_joinment_draft_internal.tsv',
        f'{input_base}_all_draft_primers_draft_internal.tsv',
        # Legacy .xls names
        f'{input_base}_all_draft_primers_after_specificity.xls',
        f'{input_base}_all_draft_primers_after_SNPs.xls',
        f'{input_base}_all_draft_primers_after_joinment.xls',
        f'{input_base}_all_draft_primers.xls',
    ]
    for path in candidates:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return path
    return None


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

                # Find output info file (try TSV multi-sheet naming, then single, then legacy .xls)
                info_base = f'{input_base}_{run_name}_primers_combination_1_info'
                info_file = None
                for candidate in [
                    f'{info_base}_ngs_primerplex_internal_primers.tsv',
                    f'{info_base}.tsv',
                    f'{info_base}.xls',
                ]:
                    if os.path.exists(candidate):
                        info_file = candidate
                        break

                if info_file is None:
                    print(f'WARNING: Expected output not found: {info_base}.*')
                    # Use intermediate draft files as fallback (best available)
                    intermediate_draft = _find_intermediate_draft(input_base)
                    if intermediate_draft and intermediate_draft != draft_file:
                        draft_file = intermediate_draft
                        print(f'  Using intermediate draft: {draft_file}')
                    print(f'  Trying next primer number...')
                    continue

                final_info_file = info_file

                # Filter to multiplexed primers and use as draft for next iteration
                new_draft = f'{input_base}_{run_name}_draft.tsv'
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

            info_base = f'{input_base}_{run_name}_primers_combination_1_info'
            info_file = None
            for candidate in [
                f'{info_base}_ngs_primerplex_internal_primers.tsv',
                f'{info_base}.tsv',
                f'{info_base}.xls',
            ]:
                if os.path.exists(candidate):
                    info_file = candidate
                    break

            if info_file is None:
                print(f'WARNING: Expected output not found: {info_base}.*')
                intermediate_draft = _find_intermediate_draft(input_base)
                if intermediate_draft and intermediate_draft != draft_file:
                    draft_file = intermediate_draft
                    print(f'  Using intermediate draft: {draft_file}')
                continue

            final_info_file = info_file

            new_draft = f'{input_base}_{run_name}_draft.tsv'
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
