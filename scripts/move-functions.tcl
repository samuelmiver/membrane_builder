
# rotate an atom selection about its center of mass
proc rotate_on_com {selection axis angle} {
	set trans_back [center_of_mass $selection]
	set trans [vecinvert $trans_back]
	$selection moveby $trans
	$selection move [transaxis $axis $angle]
	$selection moveby $trans_back
}

# Move an atom selection along a specified axis
proc translation { sel axis trans } {
  if {$axis == "x"} {set a $trans; set b 0; set c 0}
  if {$axis == "y"} {set a 0; set b $trans; set c 0}
  if {$axis == "z"} {set a 0; set b 0; set c $trans}
  $sel move [transoffset "$a $b $c" ]
} 

# center of mass calculation needed for the rotation procedure
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

