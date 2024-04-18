#################### Inputs ####################
# 0 --> psf file of the molecule
# 1 --> dcd file of the trajectory
# 2 --> Output file to be written

# Ranges of regions to compute the distances (20-63, 185-227, 332-270)
set ranges {20 63 185 227 232 270}

proc point_to_point_distance {pointA pointB} {
	set Ax [lindex $pointA 0 0]
	set Ay [lindex $pointA 0 1]
	set Az [lindex $pointA 0 2]
	set Bx [lindex $pointB 0 0]
	set By [lindex $pointB 0 1]
	set Bz [lindex $pointB 0 2]
	set L [expr {sqrt(pow($Bx-$Ax,2) + pow($By-$Ay,2) + pow($Bz-$Az,2))}]
	return $L
}

# Load Trajectory
mol new [lindex $argv 0] type psf
mol addfile [lindex $argv 1] type dcd first 0 last -1 step 1 waitfor all

# Get segnames of the simulation
set segnames [lsort -uniq [[atomselect top "protein"] get segname]]

# Prepare Output File 
set f_r_out [lindex $argv 2]
set outfile [open $f_r_out w]

# get number of frames 
set nf [molinfo top get numframes];

# Iterate through all the ranges
foreach {start end} $ranges {
	# Iterate through all residues within the range
	for {set i $start} {$i < $end} {incr i} {
		set a [atomselect top "protein and segname [lindex $segnames 0] and resid $i and name CA"]
		set b [atomselect top "protein and segname [lindex $segnames 1] and resid $i and name CA"]
		if {[$a get resname] == "GLY"} {
			set c [atomselect top "protein and segname [lindex $segnames 0] and resid $i and name HA1"]
                	set d [atomselect top "protein and segname [lindex $segnames 1] and resid $i and name HA1"]
		} else {
			set c [atomselect top "protein and segname [lindex $segnames 0] and resid $i and name CB"]
			set d [atomselect top "protein and segname [lindex $segnames 1] and resid $i and name CB"]
		}

    # Iterate through the simulation frames
		for {set fr 0} {$fr < $nf} {incr fr} {
			$a frame $fr
			set coma [$a get {x y z}]		
			$b frame $fr
			set comb [$b get {x y z}]
			set disab [point_to_point_distance $coma $comb]

                        $c frame $fr
                        set comc [$c get {x y z}]
                        $d frame $fr
                        set comd [$d get {x y z}]
                        set discd [point_to_point_distance $comc $comd]

			puts $outfile "$i,$fr,$disab,$discd"
		}
	}
}
close $outfile
quit
