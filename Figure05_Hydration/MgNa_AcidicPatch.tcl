
#################### Inputs ####################
# 0 --> psf file of the molecule
# 1 --> dcd file of the trajectory
# 2 --> Output file to be written

# Get segnames of the protein
set segnames [lsort -uniq [[atomselect top "protein"] get segname]]

# Load Trajectory
mol new [lindex $argv 0] type psf
mol addfile [lindex $argv 1] type dcd first 0 last -1 step 1 waitfor all

# Get number of Frames of the simulation
set nf [molinfo top get numframes]

# Prepare file to be written
set f_r_out [lindex $argv 2]
set outfile [open $f_r_out w]
puts $outfile "Frame,Chain,Mg,Sod,Total" 

# Cut off distance to search for Mg or Na
set cutoff 3.2

# Iterate through the simulation
for {set i 0} {$i < $nf} {incr i} {
  # ITerate through the chains
	foreach seg $segnames {
		set seltext "all within $cutoff of (protein and segname $seg and ((resid 136 to 154) or (resid 123 to 128)))"
		set sel1 [atomselect top $seltext frame $i]
		set counter_Mg [llength [lsearch -all -exact [$sel1 get type] MG]]
		set counter_Sod [llength [lsearch -all -exact [$sel1 get type] SOD]]
		set counter [expr $counter_Mg + $counter_Sod]
		puts $outfile "$i,$seg,$counter_Mg,$counter_Sod"
		$sel1 delete
	}
}
close $outfile
quit
