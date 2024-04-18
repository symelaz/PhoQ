#################### Inputs ####################
# 0 --> psf file of the molecule
# 1 --> dcd file of the trajectory
# 2 --> pdb file of the reference frame
# 3 --> Output file to be written
# 4 --> path for la1.0 library
# 5 --> path for orient library

# Set the points that define the slices for TM domain 
set points {28 202 21 210 36 194}

# Required libraries for the principle component analysis
lappend auto_path [lindex $argv 4]
lappend auto_path [lindex $argv 5]

package require Orient
namespace import Orient::orient

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

# Functions that checks if the angle between the surface vector and the z-axis is greater than pi/2
# and returns either the z or -z vector
proc check_absolute_value_greater_than_pi {variable_value} {
    # Define the value of pi
    set pi 3.14159265359

    # Get the absolute value of the variable_value
    if {$variable_value < 0} {
        set abs_value [expr {-$variable_value}]
    } else {
        set abs_value $variable_value
    }

    # Check if the absolute value of the variable is greater than pi/2
    if {$abs_value > [expr {$pi/2}]} {
        return [list 0 0 -1]
    } else {
        return [list 0 0 1]
    }
}

proc radians_to_degrees {angle_rad} {
    # Define the pi
    set pi 3.14159265359
    
    # Calculate the angle in degrees
    set angle_deg [expr {$angle_rad * 180 / $pi}]
    
    # Return the angle in degrees
    return $angle_deg
}

proc get_rotation_matrix_to_z {vector} {
	# normalize vector
	set vector_normalized [vecnorm $vector]
	
	# Get the anlge using the dot product of the vectors
	# Remeber: you don't have to devide by the product of the vectors length because they are both normalized
	set angle [expr {acos([vecdot $vector_normalized {0 0 1}])}]
	
	# Check if the angle bewteen the surface vector and the z-axis is greater than 90 degrees and
	# recalculate for the -z axis
	set z_axis_vector [check_absolute_value_greater_than_pi $angle]
	set angle [expr {acos([vecdot $vector_normalized $z_axis_vector])}]

	# Get sin and cosine of the angle
	set c [expr {cos($angle)}]
	set s [expr {sin($angle)}]

	# Define the rotation matrix around the z-axis
	set rotation_matrix [list [list $c -$s 0 0] [list $s $c 0 0] [list 0 0 1 0] [list 0 0 0 1]]
	return $rotation_matrix
}

proc write_coordinates_to_file {point_A point_B fr segnames outfile angle dz} {

	set a [atomselect top "protein and segname [lindex $segname 0] and resid $point_A and name CA" frame $fr]
        set b [atomselect top "protein and segname [lindex $segname 0] and resid $point_B and name CA" frame $fr]
        set c [atomselect top "protein and segname [lindex $segname 1] and resid $point_A and name CA" frame $fr]
        set d [atomselect top "protein and segname [lindex $segname 1] and resid $point_B and name CA" frame $fr]
	
	set coma [$a get {x y z}]
        set comb [$b get {x y z}]
        set comc [$c get {x y z}]
        set comd [$d get {x y z}]

        set disab [point_to_point_distance $coma $comb]
        set disbc [point_to_point_distance $comb $comc]
        set discd [point_to_point_distance $comc $comd]
        set disda [point_to_point_distance $comd $coma]

        puts $outfile "$point_A $point_B $fr [lindex $coma 0] [lindex $comb 0] [lindex $comc 0] [lindex $comd 0] $disab $disbc $discd $disda"
}

proc get_slice {point_A point_B nf segnames outfile} {
	
	# Get the reference slice
	set slice_ref [atomselect 0 "protein and resid $point_A $point_B and name CA"] 

	# Select the slice as defined by the 4 atoms
	set slice [atomselect top "protein and resid $point_A $point_B and name CA"]

	for {set fr 0} {$fr < $nf} {incr fr} {
		# Print current frame number
		puts "$fr out of [expr {$nf -1}]"
		# Define the selection frame
		$slice frame $fr
		
		# Align the slice to the zero frame
		set trans_mat [measure fit $slice $slice_ref]
		$slice move $trans_mat
		# $slice update
		
		# Draw the principal component to get the z axis component as a vector variable
		set slice_axis [draw principalaxes $slice]
		set vector [lindex $slice_axis 0]
		
		# Get rotation matrix of the vector and move selection
		set rotation_matrix [get_rotation_matrix_to_z $vector]
  		$slice move $rotation_matrix
		
		# Save the coordinates to a file
		write_coordinates_to_file $point_A $point_B $fr $segnames $outfile
	}
}

# Load reference molecule
mole new [lindex $argv 2]

# Load Trajectory
mol new [lindex $argv 0] type psf
mol addfile [lindex $argv 1] type dcd first 0 last -1 step 1 waitfor all

# Get segname
set segnames [lsort -uniq [[atomselect top "protein"] get segname]]

# Get number of Frames of the simulation
set nf [molinfo top get numframes]

# Prepare output file to be written
set output_file [lindex $argv 3]
set outfile [open $output_file w]
puts $outfile "PointA PointB Frame Ax Ay Az Bx By Bz Cx Cy Cz Dx Dy Dz ab bc cd da"

# Compute the positions of all slices 
foreach {i j} $points {
	puts "================================ Running for slice $i - $j ================================"
	get_slice $i $j $nf $segnames $outfile
}

# Close the output file
close $outfile

quit
