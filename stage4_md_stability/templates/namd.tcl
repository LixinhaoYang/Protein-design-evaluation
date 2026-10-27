# Experimental NPT segment. Run length and input/output settings are supplied above.
# This is a configurable execution template, not a calibrated equilibration protocol.
set case_dir [file normalize [file join [file dirname [info script]] ..]]
set output_prefix [file join $case_dir out $segment]
file mkdir [file dirname $output_prefix]
if {[llength [glob -nocomplain -directory [file dirname $output_prefix] "${segment}.*"]]} {
    error "output prefix already has files; choose a new segment name"
}
foreach path [concat [list [file join $case_dir system.psf] [file join $case_dir system.pdb]] $parameter_files] {
    if {![file isfile $path] || [file size $path] == 0} {error "missing or empty input file: $path"}
}
structure [file join $case_dir system.psf]
coordinates [file join $case_dir system.pdb]
paraTypeCharmm on
foreach path $parameter_files {parameters $path}

if {$restart_prefix ne ""} {
    if {$restart_relative} {set restart_prefix [file join $case_dir $restart_prefix]}
    foreach suffix {.coor .vel .xsc} {
        if {![file isfile ${restart_prefix}${suffix}]} {error "missing restart file: ${restart_prefix}${suffix}"}
    }
    set handle [open ${restart_prefix}.xsc r]
    set restart_ts -1
    while {[gets $handle line] >= 0} {
        if {[string trim $line] eq "" || [string match "#*" [string trimleft $line]]} {continue}
        set restart_ts [lindex $line 0]
    }
    close $handle
    if {![string is integer -strict $restart_ts] || $restart_ts < 0} {error "invalid restart timestep"}
    binCoordinates ${restart_prefix}.coor
    binVelocities ${restart_prefix}.vel
    extendedSystem ${restart_prefix}.xsc
    firsttimestep $restart_ts
} else {
    lassign $box_lengths a b c
    cellBasisVector1 $a 0.0 0.0
    cellBasisVector2 0.0 $b 0.0
    cellBasisVector3 0.0 0.0 $c
    cellOrigin {*}$box_origin
    firsttimestep 0
    temperature $temperature_k
}

# Output wrapping is off: NAMD 3.0.2 segfaults while wrapping a solvated box once
# many molecules have crossed the boundary since the previous frame, which happens
# at any realistic dcdfreq. Wrapping only rewrites output coordinates and does not
# affect the dynamics, so wrap at analysis time instead (VMD "pbc wrap",
# MDAnalysis transformations).
wrapAll off
wrapNearest off
exclude scaled1-4
oneFourScaling 1.0
cutoff 12.0
switching on
switchdist 10.0
pairlistdist 14.0
vdwForceSwitching on
timestep $timestep_fs
rigidBonds all
useSettle on
nonbondedFreq 1
fullElectFrequency 2
stepspercycle 10
PME yes
PMEGridSpacing 1.0
langevin on
langevinDamping 1.0
langevinTemp $temperature_k
# Group-based pressure is required with rigidBonds; NAMD forces it on regardless,
# so state it explicitly rather than relying on that override.
useGroupPressure yes
useFlexibleCell no
useConstantArea no
langevinPiston on
langevinPistonTarget $pressure_bar
langevinPistonPeriod 100.0
langevinPistonDecay 50.0
langevinPistonTemp $temperature_k
outputName $output_prefix
binaryoutput yes
binaryrestart yes
restartfreq $output_every
dcdfreq $output_every
xstFreq $output_every
outputEnergies $energy_every
outputPressure $energy_every

# GPU-resident integration keeps the trajectory on the GPU between steps. It needs a
# CUDA build and a single GPU, and silently disables features it cannot support, so
# it stays opt-in. "CUDASOAintegrate" is the spelling NAMD 3.0.x accepts; 3.1 renames
# it to "GPUresident".
if {$gpu_resident} {
    CUDASOAintegrate on
}

if {$restart_prefix eq "" && $minimize_steps > 0} {
    minimize $minimize_steps
    reinitvels $temperature_k
}
run $run_steps
