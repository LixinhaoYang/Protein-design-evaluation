# Stage 4: MD preparation and run setup — experimental

> **Works on my machine; not guaranteed to work on other machines or clusters.**

The original workflow was used in the author's environment. This version extracts its reusable preparation and execution logic, removes site-specific settings, and fixes several portability issues. It has synthetic regression tests, but has not been validated end to end on real systems across VMD/NAMD installations or clusters.

This stage prepares selected designs for optional MD simulations after fold/motif recovery. It provides system building, NAMD configuration, and Slurm script generation. Trajectory analysis and integration of MD metrics into the triage table are not implemented yet. A prepared system or a completed simulation is not assigned a stability `PASS`.

## Requirements and scope

- Python 3.10+; the Stage 4 commands use the standard library.
- VMD with PSFGEN, Solvate, and Autoionize plugins for system preparation. Supply its executable with `--vmd` or `VMD_BIN`.
- Your own compatible CHARMM topology and parameter files. Files are not bundled, and parameter loading order is explicit.
- A single-node multicore NAMD executable on the simulation host. The Slurm launcher uses one process with `+p` CPU threads; MPI/multi-node launchers are outside this template's scope.
- Slurm only if using the provided job-script generator.

The preparer accepts one selected chain of standard amino acids in **PDB** format. Convert AF3 mmCIF files to PDB first and verify chain IDs and residue numbering. Select whether the starting structure is from the design method or AF3 using `--source design|af3`; it is never inferred from the filename.

Other chains and water are excluded. Selected-chain ligands, modified residues, alternate residue identities, and multiple models are rejected. Disulfide patching, glycans, membranes, cofactors, and task-specific protonation are not supported. Inputs containing `SSBOND` records are rejected; cysteines without those records are treated as unpatched cysteines. Supply only systems for which that chemistry is appropriate.

The selected chain is rewritten as chain A with consecutive residue numbers and built as PSF segment `PRO`. `residue_map.csv` records the original chain, residue number, and insertion code for future motif mapping. Multiple fragments separated by `TER` in the selected chain are rejected. PSFGEN uses the original workflow's NTER/CTER patches and maps HIS to HSE; no pH-based protonation assignment is performed. Input hydrogens are removed before PSFGEN rebuilds them.

## 1. Prepare a system

Run from the repository root. All paths below are placeholders for your own inputs and installed force field:

```bash
python stage4_md_stability/prepare.py inputs/af3_pdbs/candidate_001.pdb \
  --source af3 --sequence-id candidate_001 --chain A \
  --topology /path/to/top_all36_prot.rtf \
  --vmd /path/to/vmd \
  --outdir results/md_batch
```

This runs VMD preparation, not MD or job submission. The output directory must be new or empty; a failed preparation leaves its logs for inspection and must be retried in a fresh directory. It never reuses existing structures as evidence of a successful build.

Defaults are 25 Å water padding on every bounding-box face and 0.15 mol/L salt via Autoionize with potassium/chloride ions. Change these with `--padding` and `--salt-molar`. HMR is opt-in: add `--hmr` to repartition solute hydrogen masses to 3.024 Da, leaving water masses unchanged. The Tcl code verifies each hydrogen has one heavy-atom partner and checks total mass conservation.

For batches, use multiple PDB paths with `--source`, or a case CSV:

```csv
sequence_id,input_structure,structure_source,chain_id
candidate_001,af3_pdbs/candidate_001.pdb,af3,A
candidate_002,models/candidate_002.pdb,design,A
```

```bash
python stage4_md_stability/prepare.py \
  --case-csv inputs/md_cases.csv \
  --topology /path/to/top_all36_prot.rtf \
  --vmd /path/to/vmd --outdir results/md_batch
```

Paths in the case CSV are relative to that CSV. Candidate IDs must be unique, start with a letter or digit, and contain at most 80 letters, digits, dots, underscores, or hyphens. Use the same IDs as Stages 1–3. `chain_id` is optional and defaults to `--chain A`; use `_` for a blank PDB chain.

```text
results/md_batch/
  systems.json
  candidate_001/
    case.json                 # source, box, HMR, topology provenance
    protein.pdb               # selected, renumbered heavy-atom protein
    residue_map.csv
    system.psf
    system.pdb                # solvated/ionized system
    conf/
    out/
    logs/prepare.log
```

Intermediate preparation files remain inside each case directory. Only successful builds receive `case.json` with `state=prepared`; `systems.json` is written after the full batch succeeds. These metadata describe preparation, not a scientific stability verdict.

## 2. Generate NAMD configurations

Transfer the prepared batch to the simulation host if needed, then generate configurations there. Supply your actual compatible parameter files, in the required loading order:

```bash
python stage4_md_stability/configure.py \
  --manifest results/md_batch/systems.json \
  --segment equilibration --run-steps 250000 --minimize-steps 2000 \
  --parameter /path/to/protein.prm \
  --parameter /path/to/water_ions.prm \
  --parameter /path/to/nbfix.prm
```

This writes `candidate_001/conf/equilibration.conf` for each system, without starting NAMD. Use repeated `--sequence-id` options to select a subset. Parameter paths must exist on the NAMD host; configuration generation does not inspect their contents. Generated configurations check their presence at run time.

The example step counts illustrate configuration syntax; they are **not** a validated equilibration schedule. Choose minimization, equilibration, restraints, replicates, and production duration for your system. The supplied NPT template has no positional restraints or heating schedule; naming a segment `equilibration` does not establish equilibration.

| Setting | Default |
| --- | --- |
| Timestep | 2 fs without HMR; 4 fs with HMR |
| Minimization | 200 steps for a fresh start; 0 for a restart |
| Run length | Required `--run-steps`, measured as additional steps |
| Temperature / pressure | 300 K / 1.01325 bar |
| DCD / restart / cell interval | 25,000 steps: 50 ps at 2 fs, 100 ps at 4 fs |
| Energy / pressure reporting | Every 5,000 steps |

`--timestep-fs`, `--temperature`, `--pressure-bar`, `--output-every`, and `--energy-every` override these values. Timesteps above 2 fs require recorded HMR preparation and cannot exceed 4 fs. This check is not evidence that 4 fs is suitable for every system. Step counts and output intervals must be multiples of the template's 10-step cycle.

Each segment writes a distinct `out/<segment>` prefix. Existing configurations and output files are not overwritten. NAMD configuration paths are quoted for Tcl, including spaces and special characters; the generated configuration derives its case directory at run time.

### Continue from a restart

After a previous run has produced a consistent `equilibration.restart.coor`, `.vel`, and `.xsc` set:

```bash
python stage4_md_stability/configure.py \
  --manifest results/md_batch/systems.json \
  --segment production --restart-from equilibration --run-steps 5000000 \
  --parameter /path/to/protein.prm \
  --parameter /path/to/water_ions.prm \
  --parameter /path/to/nbfix.prm
```

Restart mode loads coordinates, velocities, and periodic cell state, and resumes timestep numbering from the XSC file. It does not minimize or reinitialize velocities. `--run-steps` adds time to the restart step. If a job stopped partway through production, continue under a new segment name, for example `production_2 --restart-from production`.

For one selected system, `--restart-prefix /path/to/existing.restart` can reference an external restart set instead. Binary restart/topology compatibility and scientific equilibration remain the user's responsibility. See the [NAMD restart example](https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node27.html).

## 3. Generate Slurm scripts

```bash
python stage4_md_stability/slurm.py results/md_batch/*/conf/equilibration.conf \
  --namd /path/on/cluster/to/namd3 \
  --time 24:00:00 --cpus 6 --gpus 1 \
  --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

By default this only writes `.slurm` files and prints submission commands. It does not call `sbatch`. Optional flags include `--qos`, `--gpu-type`, and `--dependency JOB_ID[:JOB_ID]` (using `afterok`). There are no embedded accounts, email addresses, executable paths, or GPU models. Regenerate scripts after moving the batch, because scheduler script paths are absolute.

Submit a generated script yourself using its printed `sbatch` command. Alternatively, repeat the generation command with `--submit` to submit explicitly. Identical scripts can be reused; existing scripts with different settings are not overwritten. A scheduler dependency orders jobs; it does not create restart inputs or configure continuation. Completed jobs and existing output prefixes must not be resubmitted as fresh starts.

## Validation and remaining work

```bash
python -m unittest discover -s tests -v
```

Tests use synthetic PDB records, mocked VMD preparation, Tcl helpers, and generated configuration/job scripts. They do not run VMD plugins, NAMD simulations, or real scheduler submissions. Test the complete pipeline with your own installation and inspect preparation and simulation logs before scaling up.

Trajectory imaging/alignment, RMSD/RMSF and motif analysis, convergence/replicate assessment, MD summary CSVs, calibrated selection rules, and a `--md-csv` aggregation interface remain future work. No real project structures, trajectories, restart files, or force-field files are included in this repository.
