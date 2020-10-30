proc move_to_origin {selection} {

        set gc [veczero] ;# make a vector
        set sel $selection ;# make a selection
        foreach coord [$sel get {x y z}] {set gc [vecadd $gc $coord]} ;# add all coordinates to vector
        set gc [vecscale [expr 1.0 /[$sel num]] $gc] ;# normalize vector
        $sel moveby [vecscale -1 $gc] ;# move to origin
}
