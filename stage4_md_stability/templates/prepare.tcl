# Variables are supplied by the per-case build.tcl driver.
package require psfgen
package require solvate
package require autoionize
resetpsf
foreach path $topology_files {topology $path}
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# One normalized protein chain, with the original workflow's terminal patches.
segment PRO {
    first NTER
    pdb protein.pdb
    last CTER
}
coordpdb protein.pdb PRO
guesscoord
regenerate angles dihedrals
writepsf protein_H.psf
writepdb protein_H.pdb

set protein_mol [mol new protein_H.psf waitfor all]
mol addfile protein_H.pdb molid $protein_mol waitfor all
set protein [atomselect $protein_mol all]
if {[$protein num] == 0} {error "PSFGEN produced no protein atoms"}
lassign [measure minmax $protein] minimum maximum
lassign [md_water_box $minimum $maximum $water_padding] center lengths lower upper
# Center the bounding box, not the atom centroid, so padding is correct on all faces.
$protein moveby [vecscale -1 $center]
$protein writepdb protein_centered.pdb
$protein delete
mol delete $protein_mol

solvate protein_H.psf protein_centered.pdb -minmax [list $lower $upper] -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -sc $salt_molar -cation POT -anion CLA -o ionized
set system_mol [mol new ionized.psf waitfor all]
mol addfile ionized.pdb molid $system_mol waitfor all
if {$use_hmr} {md_repartition_hydrogens $system_mol 3.024}
set all [atomselect $system_mol all]
$all writepsf system.psf
$all writepdb system.pdb
$all delete
mol delete $system_mol
set metadata [open box.dat w]
puts $metadata "[join $lengths { }] 0.0 0.0 0.0"
close $metadata
