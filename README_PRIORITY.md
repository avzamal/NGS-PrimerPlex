# NGS-PrimerPlex: Priority-Based Multiplex Primer Design

Extended fork of [NGS-PrimerPlex](https://github.com/aakechin/NGS-PrimerPlex) with priority-based region handling, iterative selection, smart primer re-selection, and force-include support.

## Installation

### Prerequisites

```bash
# BWA (for --do-blast specificity checks)
sudo apt install bwa

# Python dependencies
pip3 install biopython primer3-py pysam networkx numpy xlrd==1.2.0 xlsxwriter lxml openpyxl
```

### Clone

```bash
git clone https://github.com/avzamal/NGS-PrimerPlex.git
```

### BWA Index (if not already done)

```bash
bwa index reference_genome.fna
```

---

## Input File Format

Tab-separated file with 7 or 8 columns:

```
Chrom	Start	End	Name	MinAmpLen	MaxAmpLen	EndShift	Priority
chr1	1000	1200	marker_A	120	200	10	1
chr2	5000	5150	marker_B	120	200	10	2
chr3	8000	8100	marker_C	120	200	10	2
```

| Column | Description |
|--------|-------------|
| Chrom | Chromosome name (must match reference genome) |
| Start | Region start (0-based) |
| End | Region end |
| Name | Unique region name |
| MinAmpLen | Minimum amplicon length |
| MaxAmpLen | Maximum amplicon length |
| EndShift | How far primer can extend beyond region |
| Priority | **Optional.** 1 (highest), 2, or 3 (lowest). Default: 1 |

**Backwards compatible:** Existing 7-column files work unchanged (all regions default to priority 1).

---

## Priority System

When regions have mixed priorities, the multiplex assignment works in phases:

1. **Phase 1:** Place priority-1 amplicons using the full clique-based algorithm (guaranteed best placement)
2. **Phase 2:** Insert priority-2 amplicons into existing multiplex pools (only if compatible with all pool members)
3. **Phase 3:** Insert priority-3 amplicons the same way

This ensures high-priority markers are never displaced by lower-priority ones.

---

## Usage

### Basic Run (Single Iteration)

```bash
python3 NGS_primerplex.py \
    --regions-file regions.tsv \
    --reference-genome genome.fna \
    --primers-number1 30 \
    --run-name my_panel \
    --do-blast \
    --threads 8
```

### With Smart Re-selection

Identifies primers that are incompatible with many others and tries alternative primer pairs:

```bash
python3 NGS_primerplex.py \
    --regions-file regions.tsv \
    --reference-genome genome.fna \
    --primers-number1 30 \
    --run-name my_panel \
    --do-blast \
    --threads 8 \
    --smart-reselection
```

### With Force-Include

Force specific regions into the multiplex. The tool reports blockers and suggests removals if placement fails:

```bash
# Create a force-include file (one region name per line)
cat > force_include.txt << 'EOF'
critical_marker_1
critical_marker_2
EOF

python3 NGS_primerplex.py \
    --regions-file regions.tsv \
    --reference-genome genome.fna \
    --primers-number1 30 \
    --run-name my_panel \
    --do-blast \
    --threads 8 \
    --force-include force_include.txt \
    --force-include-max-primers 1000
```

Force-include regions are internally set to priority 0 (highest). If they cannot be placed, the output shows:
- Which existing amplicons are blocking them
- Which blockers are lowest-priority (candidates for removal)
- Whether manual review is needed (all blockers are high-priority)

---

## Iterative Selection

The `iterative_primerplex.py` script automates running NGS-PrimerPlex multiple times with escalating `--primers-number1` values. Between iterations, it filters successfully multiplexed primers and uses them as draft input for the next run.

### Basic Iterative Run

```bash
python3 iterative_primerplex.py \
    --regions-file regions_with_priority.tsv \
    --reference-genome genome.fna \
    --primer-nums 10,30,100,300,600 \
    --run-name panel_v1 \
    --do-blast \
    --threads 8
```

### Iterative Run with Priority Levels

When regions have mixed priorities, `--by-priority` (default) runs iterations per priority level:

```bash
python3 iterative_primerplex.py \
    --regions-file regions_with_priority.tsv \
    --reference-genome genome.fna \
    --primer-nums 10,30,100,300 \
    --run-name panel_v1 \
    --by-priority \
    --smart-reselection \
    --do-blast \
    --threads 8
```

This will:
1. Run iterations for priority-1 regions only (escalating primer numbers)
2. Once all priority-1 regions are placed, add priority-2 regions
3. Continue escalating until all regions are placed or primer numbers exhausted

### Without Priority Separation

```bash
python3 iterative_primerplex.py \
    --regions-file regions.tsv \
    --reference-genome genome.fna \
    --primer-nums 10,30,100,300 \
    --run-name panel_v1 \
    --no-by-priority \
    --do-blast \
    --threads 8
```

### Iterative Script Arguments

| Argument | Description |
|----------|-------------|
| `--regions-file` | Input TSV with regions |
| `--reference-genome` | Reference genome FASTA |
| `--primer-nums` | Comma-separated escalating values (default: `10,30,100,300,600`) |
| `--run-name` | Base name for output files |
| `--output-dir` | Output directory (default: same as input) |
| `--by-priority` | Run iterations per priority level (default: True) |
| `--no-by-priority` | Run all regions together |

All other arguments are passed through to `NGS_primerplex.py` (e.g., `--do-blast`, `--threads`, `--smart-reselection`, `--force-include`).

---

## New Command-Line Arguments

| Argument | Description |
|----------|-------------|
| `--smart-reselection` | Enable smart primer re-selection (replace problematic primers with better alternatives) |
| `--force-include FILE` | File with region names that must be in the multiplex |
| `--force-include-max-primers N` | Max primers to try for force-include regions (default: 1000) |

---

## Output

Standard NGS-PrimerPlex output files are produced:
- `*_primers_combination_N.xls` — Primer sequences
- `*_primers_combination_N_info.xls` — Detailed info including multiplex assignment
- `*_primers_combination_N.fa` — FASTA of primer sequences
- `*_primers_combination_N_internal_amplicons.fa` — Internal amplicon sequences

The `Designed_Multiplex` column in the info file indicates which multiplex pool each primer pair was assigned to. Empty means the primer was not placed in any pool.

When using priority-based assignment, the console output shows placement statistics per priority level.

---

## Example Workflow

```bash
# 1. Prepare input with priorities
#    42 critical markers (priority 1) + 159 additional markers (priority 2)

# 2. First iterative run with priority separation
python3 iterative_primerplex.py \
    --regions-file combined_markers.tsv \
    --reference-genome genome.fna \
    --primer-nums 10,30,100,300 \
    --run-name panel_v1 \
    --by-priority \
    --smart-reselection \
    --do-blast \
    --threads 8

# 3. Review output — check how many markers were placed per priority level

# 4. If a new critical marker needs to be added later:
echo "new_critical_marker" > force_include.txt
python3 NGS_primerplex.py \
    --regions-file combined_markers.tsv \
    --reference-genome genome.fna \
    --primers-number1 100 \
    --run-name panel_v2 \
    --draft-primers panel_v1_draft.xls \
    --force-include force_include.txt \
    --smart-reselection \
    --do-blast \
    --threads 8
```

---

## Changes from Original NGS-PrimerPlex

1. **networkx 2.x compatibility** — Fixed deprecated `attr_dict`, `nodes()` NodeView, `neighbors()` iterator, frozen subgraph views
2. **xlrd + Python 3.12 compatibility** — Fixed `Element_has_iter` issue with defusedxml
3. **Priority column (column 8)** — Optional region priority (1-3) in input TSV
4. **Priority-aware multiplex assignment** — Higher-priority regions placed first
5. **Smart primer re-selection** (`--smart-reselection`) — Replaces low-compatibility primers
6. **Force-include** (`--force-include`) — Guarantees specific regions in multiplex or reports blockers
7. **Iterative runner** (`iterative_primerplex.py`) — Automated escalating primer number runs
