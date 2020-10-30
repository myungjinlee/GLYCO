set inpdb [lindex $argv 0]
set outpdb [lindex $argv 1]
mol load pdb $inpdb
move_to_pro_origin [atomselect top "protein"] [atomselect top "all"]
set sel [atomselect top "all"]
$sel writepdb $outpdb
mol delete all
exit
