#################### Inputs ####################
# 0 --> pdb file
# 1 --> output pdb file

# Load file 
mol new [lindex $argv 0] type pdb

# Get resids of chain A
set A [atomselect top "protein and chain A"]
set ResidLast [lindex [$A get resid] end]

set B [atomselect top "protein and chain B"]
set ResidStart [lindex [$B get resid] 0]
foreach i [$B get index] {
	set temp [atomselect top "index $i and chain B" ]
	$temp set chain "A"
	$temp set resid [expr [$temp get resid] + $ResidLast +1]
}

set total [atomselect top "protein and chain A"]
$total writepdb [lindex $argv 1] 

quit
