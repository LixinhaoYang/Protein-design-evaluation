# Solute-only HMR; water masses are unchanged. Applied only when requested.
proc md_hmr_pair {hydrogen_mass heavy_mass target_mass} {
    set delta [expr {$target_mass - $hydrogen_mass}]
    set remaining [expr {$heavy_mass - $delta}]
    if {$delta <= 0 || $remaining <= 1.0} {error "invalid HMR mass transfer"}
    return [list $target_mass $remaining]
}

proc md_single_partner {bondlists} {
    if {[llength $bondlists] != 1 || [llength [lindex $bondlists 0]] != 1} {
        error "each solute hydrogen must have exactly one bonded partner"
    }
    return [lindex $bondlists 0 0]
}

proc md_total_mass {selection} {
    set total 0.0
    foreach mass [$selection get mass] {set total [expr {$total + $mass}]}
    return $total
}

proc md_repartition_hydrogens {molid target_mass} {
    set all [atomselect $molid all]
    set before [md_total_mass $all]
    set hydrogens [atomselect $molid "mass > 0.99 and mass < 1.1 and not water"]
    if {[$hydrogens num] == 0} {error "no solute hydrogens found for HMR"}
    foreach index [$hydrogens get index] {
        set hydrogen [atomselect $molid "index $index"]
        set partner [md_single_partner [$hydrogen getbonds]]
        set heavy [atomselect $molid "index $partner"]
        if {[$heavy num] != 1 || [lindex [$heavy get mass] 0] < 2.0} {
            error "hydrogen $index is not bonded to one heavy atom"
        }
        lassign [md_hmr_pair [lindex [$hydrogen get mass] 0] \
            [lindex [$heavy get mass] 0] $target_mass] hmass heavy_mass
        $hydrogen set mass $hmass
        $heavy set mass $heavy_mass
        $hydrogen delete
        $heavy delete
    }
    if {abs([md_total_mass $all] - $before) > 0.0001} {
        error "HMR did not conserve total system mass"
    }
    $hydrogens delete
    $all delete
}
