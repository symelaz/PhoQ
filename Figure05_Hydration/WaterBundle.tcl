
#################### Inputs ####################
# 0 --> psf file of the molecule
# 1 --> dcd file of the trajectory
# 2 --> Output file to be written

# Prepare Output file
set f_r_out [lindex $argv 2]
puts "Output File to be Written $f_r_out"
set cutoff 4; #3.5 or 4

# Load Trajectory
mol new [lindex $argv 0] type psf
mol addfile [lindex $argv 1] type dcd first 0 last -1 step 1 waitfor all

set nf [molinfo top get numframes]

set outfile [open $f_r_out w]
for {set i 0} {$i < $nf} {incr i} {
	puts "frame $i of $nf"
	set sel1 [atomselect top "water within $cutoff of (protein and resid 32 202 and sidechain)" frame $i]
	set counter [llength [lsearch -all -exact [$sel1 get type] OT]]
	set ids [lsort -unique [$sel1 get resid]]
	puts $outfile "$i,$counter,$ids"
	$sel1 delete
}
close $outfile
quit
