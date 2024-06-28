mkdir struct_in/multimers
files=$(ls struct_in/*.pdb)

for identifier in $files
do
	IFS="/" read -a prefix <<< "$identifier"
	pdb_id="${prefix[-1]}"
	IFS="." read -a prefix <<< "$pdb_id"
	vmd -dispdev text -e utils/dimer_to_monomer.tcl -args $identifier "struct_in/multimers/${prefix[0]}_multimer.pdb"
done
