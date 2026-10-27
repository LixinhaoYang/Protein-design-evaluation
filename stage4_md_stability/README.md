# Stage 4: MD preparation and run setup — experimental

> **Works on my machine; not guaranteed to work on other machines or clusters.**

The original workflow was used in the author's environment. This version retains its CHARMM36m parameter combination and NAMD 3 CUDA/Slurm defaults, while exposing command-line overrides. It has synthetic regression tests, but has not been validated end to end on real systems.

> **Update (2026-10-27):** early end-to-end testing has now been done, and the pipeline runs under Slurm with NAMD 3.0.2 CUDA. Twenty-one ~67,000-atom designed-protein systems were built with `prepare.py`, configured with `configure.py`, and run as NPT MD on A100 and V100 GPUs: minimization, stable 300 K dynamics, and trajectory/restart output all behave. One fix came out of that testing and is in the current code: the HMR mass-conservation tolerance now scales with system size. A reproducible segfault at the first trajectory frame remains under investigation and is not yet explained. Long production runs and the analysis layer are still unproven, so treat the numbers as untested and read your startup logs.

This stage prepares selected designs for optional MD simulations after fold/motif recovery. It provides system building, NAMD configuration, and Slurm script generation. Trajectory analysis and integration of MD metrics into the triage table are not implemented yet. A prepared system or a completed simulation is not assigned a stability `PASS`.

## Requirements and scope

- Python 3.10+; the Stage 4 commands use the standard library.
- VMD with PSFGEN, Solvate, and Autoionize plugins for system preparation. Supply its executable with `--vmd` or `VMD_BIN`.
- The local default force field is `top_all36_prot.rtf` plus `toppar_water_ions.str` for PSF building, and `toppar_water_ions_modified.prm`, `par_all36m_prot.prm`, and `NBFIX.str` for NAMD. They are read from `/storage/home/hcoda1/9/lyang403/r-jcgumbart-0/ParamFiles`, or from `STAGE4_PARAMDIR` when set.
- The local default NAMD executable is `/storage/project/r-jcgumbart-0/jcgumbart/software/NAMD_3.0.2_Linux-x86_64-multicore-CUDA/namd3`. It is a single-node CUDA build. The Slurm launcher runs one process with `+p` CPU threads; `+setcpuaffinity` and explicit NAMD GPU devices are opt-in (`--cpu-affinity`, `--gpu-devices`), since pinning across a non-contiguous cgroup core set is a common source of trouble. MPI/multi-node launchers are outside this template's scope.
- Slurm only if using the provided job-script generator.

The preparer accepts one selected chain of standard amino acids in **PDB** format. Convert AF3 mmCIF files to PDB first with `cif_to_pdb.py` (§0) and verify chain IDs and residue numbering. Select whether the starting structure is from the design method or AF3 using `--source design|af3`; it is never inferred from the filename.

Other chains and water are excluded. Selected-chain ligands, modified residues, alternate residue identities, and multiple models are rejected. Disulfide patching, glycans, membranes, cofactors, and task-specific protonation are not supported. Inputs containing `SSBOND` records are rejected; cysteines without those records are treated as unpatched cysteines. Supply only systems for which that chemistry is appropriate.

The selected chain is rewritten as chain A with consecutive residue numbers and built as PSF segment `PRO`. `residue_map.csv` records the original chain, residue number, and insertion code for future motif mapping. Multiple fragments separated by `TER` in the selected chain are rejected. PSFGEN uses the original workflow's NTER/CTER patches and maps HIS to HSE; no pH-based protonation assignment is performed. Input hydrogens are removed before PSFGEN rebuilds them.

## 0. Convert AF3 mmCIF to PDB

AF3 writes mmCIF, and the preparer reads PDB. `cif_to_pdb.py` performs that one conversion:

```bash
python stage4_md_stability/cif_to_pdb.py af3_output/*/*_model.cif \
  --chain A --outdir inputs/af3_pdbs
```

Each file is written as `<sequence_id>.pdb`, taking the ID from the mmCIF stem with a trailing `_model` removed; use `--sequence-id` to name a single output. Existing files are never overwritten.

Only the requested `auth_asym_id` is kept. Hydrogens are dropped because PSFGEN rebuilds them, and `OXT` is dropped because the CTER patch builds `OT1`/`OT2` itself. A HETATM record in the selected chain, a second model, an alternate conformer, a non-standard residue, or a missing `_atom_site` column is rejected rather than silently dropped, so an unexpected input fails here instead of inside PSFGEN. This tool converts; it does not select, score, or repair structures.

## 1. Prepare a system

Run from the repository root. All paths below are placeholders for your own inputs and installed force field:

```bash
python stage4_md_stability/prepare.py inputs/af3_pdbs/candidate_001.pdb \
  --source af3 --sequence-id candidate_001 --chain A \
  --outdir results/md_batch
```

This runs VMD preparation, not MD or job submission. The output directory must be new or empty; a failed preparation leaves its logs for inspection and must be retried in a fresh directory. It never reuses existing structures as evidence of a successful build.

Defaults are 25 Å water padding on every bounding-box face and 0.15 mol/L salt via Autoionize with potassium/chloride ions. Change these with `--padding` and `--salt-molar`. HMR is opt-in: add `--hmr` to repartition solute hydrogen masses to 3.024 Da, leaving water masses unchanged. The Tcl code verifies each hydrogen has one heavy-atom partner and checks total mass conservation. VMD stores masses as single-precision floats, so summing a solvated box accumulates rounding error that grows with system size; the conservation tolerance therefore scales as `max(0.0001, 1e-7 * total_mass)` rather than being a fixed absolute value. For a ~400,000 amu box that is ~0.04 amu, still far below a single hydrogen mass.

For batches, use multiple PDB paths with `--source`, or a case CSV:

```csv
sequence_id,input_structure,structure_source,chain_id
candidate_001,af3_pdbs/candidate_001.pdb,af3,A
candidate_002,models/candidate_002.pdb,design,A
```

```bash
python stage4_md_stability/prepare.py \
  --case-csv inputs/md_cases.csv \
  --outdir results/md_batch
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

Generate configurations on the simulation host. By default this uses the local water/ion, CHARMM36m protein, and NBFIX parameter combination in the order used by the original workflow:

```bash
python stage4_md_stability/configure.py \
  --manifest results/md_batch/systems.json \
  --segment equilibration --run-steps 250000
```

This writes `candidate_001/conf/equilibration.conf` for each system, without starting NAMD. Use repeated `--sequence-id` options to select a subset. Pass repeated `--parameter` options only to replace the complete default parameter list; their order is preserved. Generated configurations check every input at run time.

The MD settings in `templates/namd.tcl` reproduce the proven `base.conf` from the previous version of this workflow (`testrun/automaticMD`) line for line, and the configure defaults match its values. `tests/test_md_setup.py` pins them, so a change there should be a deliberate decision. Keep minimization at the 200-step default: every run configured with `--minimize-steps 10000` segfaulted with no error message at its first trajectory write (step 25,000), while 200-step runs on the same systems did not. Why a longer minimization does this is not understood.

`--gpu-resident` enables NAMD 3 GPU-resident integration. It is off by default because the proven configuration does not use it; on a V100 it measured only about 5–10% faster for a ~67,000-atom system. It also silently disables features it cannot support, so read the startup log before relying on it. GPU-resident (CUDA-accelerated) runs with the proven settings are planned for testing next.

The example step counts illustrate configuration syntax; they are **not** a validated equilibration schedule. Choose minimization, equilibration, restraints, replicates, and production duration for your system. The supplied NPT template has no positional restraints or heating schedule; naming a segment `equilibration` does not establish equilibration.

| Setting | Default |
| --- | --- |
| Timestep | 2 fs without HMR; 4 fs with HMR |
| GPU-resident integration | off; `--gpu-resident` enables it |
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
  --segment production --restart-from equilibration --run-steps 5000000
```

Restart mode loads coordinates, velocities, and periodic cell state, and resumes timestep numbering from the XSC file. It does not minimize or reinitialize velocities. `--run-steps` adds time to the restart step. If a job stopped partway through production, continue under a new segment name, for example `production_2 --restart-from production`.

For one selected system, `--restart-prefix /path/to/existing.restart` can reference an external restart set instead. Binary restart/topology compatibility and scientific equilibration remain the user's responsibility. See the [NAMD restart example](https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node27.html).

## 3. Generate Slurm scripts

```bash
python stage4_md_stability/slurm.py results/md_batch/*/conf/equilibration.conf \
  --time 24:00:00
```

By default this writes scripts requesting one H200 GPU, six CPU cores, the `inferno` QOS, and `gts-jcgumbart-coda20`, and runs the local CUDA NAMD 3 binary with plain `+p`. No `--partition` is emitted by default: sites that route GPU jobs by `--gres` reject an explicit one, and on the cluster this targets `embers` is a QOS, not a partition. It only writes `.slurm` files and prints submission commands—it does not call `sbatch`. Override this setup with `--namd`, `--cpus`, `--gpus`, `--gpu-type`, `--gpu-devices`, `--cpu-affinity`, `--account`, `--partition`, `--qos`, or `--dependency JOB_ID[:JOB_ID]` (using `afterok`).

Submit a generated script yourself using its printed `sbatch` command. Alternatively, repeat the generation command with `--submit` to submit explicitly. Identical scripts can be reused; existing scripts with different settings are not overwritten. A scheduler dependency orders jobs; it does not create restart inputs or configure continuation. Completed jobs and existing output prefixes must not be resubmitted as fresh starts.

## Validation and remaining work

```bash
python -m unittest discover -s tests -v
```

Tests use synthetic PDB records, mocked VMD preparation, Tcl helpers, and generated configuration/job scripts. They do not run VMD plugins, NAMD simulations, or real scheduler submissions. Test the complete pipeline with your own installation and inspect preparation and simulation logs before scaling up.

Trajectory imaging/alignment, RMSD/RMSF and motif analysis, convergence/replicate assessment, MD summary CSVs, calibrated selection rules, and a `--md-csv` aggregation interface remain future work. No real project structures, trajectories, restart files, or force-field files are included in this repository.
