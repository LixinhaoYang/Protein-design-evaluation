# Protein Design Triage

A reusable, staged triage workflow for protein designs. It is organized to remove implausible candidates before the expensive AlphaFold 3 (AF3) step.

This repository intentionally contains **no** real sequences, structures, inputs, job logs, model files, result tables, or test cases. Keep your own data outside the repository or in an ignored directory such as `data/`, `inputs/`, or `results/`.

## Workflow order

**Candidate sequences and generated candidate structures → Stage 1: sequence screening → Stage 2: fast physical screening → PASS candidates only → AF3 prediction or refolding → Stage 3: fold and motif recovery**

AF3 is deliberately **after** the fast physical screen. Do not submit a candidate to AF3 until it passes the Stage 1 + Stage 2 pre-AF3 triage table.

| Stage | Directory | Purpose |
| --- | --- | --- |
| 1 | `stage1_sequence_screening/` | Sequence validity, composition, low-complexity, hydropathy, charge, and cysteine-pattern checks. |
| 2 | `stage2_fast_physics/` | Restrained relaxation, clash checks, packing and surface proxies, and optional PyRosetta metrics for candidate structures. |
| 3 | `stage3_fold_recovery/` | AF3/refolded-structure confidence, global recovery, and motif recovery metrics. |
| Aggregate | `triage_table.py` | Joins sequential stage tables by `sequence_id`. |

The reported signals support ranking and quality control. They do not, on their own, establish or rule out function, binding, catalysis, or specificity.

## Environments and dependencies

Python 3.10 or later is required. The core environment runs Stage 1, Stage 3, `triage_table.py`, and the geometry-only Stage 2 helpers:

```bash
conda env create -f environment.yml
conda activate protein-design-evaluation
```

Alternatively, create a virtual environment and install the same core dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

`numpy` is required for Stage 3 structural alignment. `biopython` provides the full Stage 1 pI and charge calculations; without it, Stage 1 falls back to simpler charge metrics.

The full Stage 2 path requires **PyRosetta** for coordinate-restrained relaxation and official Rosetta metrics. PyRosetta is intentionally not listed in `requirements.txt`: install it separately in an environment that is licensed and configured for your site, then run the PyRosetta commands in that environment. For example:

```bash
conda run -n pyrosetta python stage2_fast_physics/restrained_relax.py --help
```

AF3 execution is also external to this repository. It requires your site's AF3 software or container, model weights, databases, GPU allocation, and scheduler configuration. This repository only consumes the resulting mmCIF and confidence JSON files in Stage 3.

## Usage

Every path below is a placeholder. Replace it with a path to your own local data. Write outputs to `results/`, which is ignored by Git.

### Stage 1: sequence screening

```bash
python stage1_sequence_screening/run.py /path/to/sequences.fasta \
  --csv results/stage1_sequence.csv

python stage1_sequence_screening/triage_features.py /path/to/sequences.fasta \
  --csv results/stage1_triage.csv
```

### Stage 2: fast physical screening

Stage 2 operates on the candidate PDB structures produced by the design method, not on AF3 structures. With PyRosetta available, relax first and score the relaxed structures. Replace `pyrosetta` below with the name of your site-managed PyRosetta environment:

```bash
conda run -n pyrosetta python stage2_fast_physics/restrained_relax.py /path/to/candidate_models/*.pdb \
  --outdir results/relaxed_structures --csv results/stage2_relax.csv

conda run -n pyrosetta python stage2_fast_physics/rosetta_official_metrics.py results/relaxed_structures/*.pdb \
  --csv results/stage2_rosetta.csv

python stage2_fast_physics/triage_features.py \
  --relax-csv results/stage2_relax.csv \
  --rosetta-csv results/stage2_rosetta.csv \
  --csv results/stage2_triage.csv
```

The geometry-only scripts, such as `clash_score.py`, remain useful when PyRosetta is unavailable. The full PyRosetta path above produces the compact Stage 2 table used by `triage_table.py`.

### Build the pre-AF3 triage table

```bash
python triage_table.py \
  --sequence-csv results/stage1_triage.csv \
  --physics-csv results/stage2_triage.csv \
  --csv results/pre_af3_triage.csv
```

Submit only rows with `triage_status=PASS` to AF3. Your AF3 input-preparation and scheduler setup remain environment-specific and are intentionally not included in this repository.

### Stage 3: fold and motif recovery after AF3

After AF3 finishes, create a case CSV with `sequence_id`, `design_structure`, `predicted_structure`, `summary_confidences`, and `motif_mapping`. Paths may be relative to the case CSV location.

```bash
python stage3_fold_recovery/structure_recovery.py \
  --case-csv /path/to/af3_cases.csv \
  --csv results/stage3_metrics.csv

python stage3_fold_recovery/triage_features.py \
  --metrics-csv results/stage3_metrics.csv \
  --csv results/stage3_triage.csv
```

`predicted_structure` currently uses an AlphaFold 3-style mmCIF parser, and `motif_mapping` requires a `candidate_position_1based` column. Stage 3 emits reference features and does not impose calibrated hard cutoffs by default.

### Build the final triage table

```bash
python triage_table.py \
  --sequence-csv results/stage1_triage.csv \
  --physics-csv results/stage2_triage.csv \
  --fold-csv results/stage3_triage.csv \
  --csv results/final_triage.csv
```

## Pre-publication check

`.gitignore` covers common real-input, structure, log, and output formats. Stage only the source code and repository metadata when committing:

```bash
git add README.md requirements.txt .gitignore triage_table.py \
  stage1_sequence_screening stage2_fast_physics stage3_fold_recovery
git status
```

Confirm that no data file is staged before connecting a remote repository or pushing. Add the license of your choice before making the repository public.
