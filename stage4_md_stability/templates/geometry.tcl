# Pure geometry helpers; these can also be tested with tclsh, without VMD.
proc md_water_box {minimum maximum padding} {
    if {[llength $minimum] != 3 || [llength $maximum] != 3 || $padding <= 0} {
        error "invalid box bounds or padding"
    }
    set center {}
    set lengths {}
    set lower {}
    set upper {}
    foreach low $minimum high $maximum {
        if {$high < $low} {error "reversed box bounds"}
        set length [expr {$high - $low + 2.0 * $padding}]
        lappend center [expr {0.5 * ($high + $low)}]
        lappend lengths $length
        lappend lower [expr {-0.5 * $length}]
        lappend upper [expr {0.5 * $length}]
    }
    return [list $center $lengths $lower $upper]
}
