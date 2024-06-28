#!/bin/bash

# Prepare structures in Multimer format
bash utils/prepare_structures.sh

# Folder with all the HAMP structures
structures=HAMPpred/data/input/struct_in/multimers

# Create a directory to save all the intermediate steps of FoldSeek search
mkdir -p Foldseek/Multimer/{Target,Query,Alignment,MSA}

# Prepare the databases of the pdbs
foldseek createdb "$structures/phoq_full_af_HAMP_multimer.pdb" Foldseek/Multimer/Query/queryDB > log
foldseek createdb $structures Foldseek/Multimer/Target/targetDB > log

# Search the HAMP domain of the alphafold predicted model against the HAMP database
foldseek search Foldseek/Multimer/Query/queryDB Foldseek/Multimer/Target/targetDB Foldseek/Multimer/Alignment/aln tmp --tmscore-threshold 0.5 -s 9.5 --min-seq-id 0.2 -c 0.8 --cov-mode 0 -a > log

# Prepare the Aligment Results into a readable tsv format
foldseek aln2tmscore Foldseek/Multimer/Query/queryDB Foldseek/Multimer/Target/targetDB Foldseek/Multimer/Alignment/aln Foldseek/Multimer/Alignment/aln_tmscore > log
foldseek createtsv Foldseek/Multimer/Query/queryDB Foldseek/Multimer/Target/targetDB Foldseek/Multimer/Alignment/aln_tmscore Foldseek/Multimer/Alignment/aln_tmscore.tsv > log

# Run the MSA between the structures
foldseek result2msa Foldseek/Multimer/Query/queryDB Foldseek/Multimer/Target/targetDB Foldseek/Multimer/Alignment/aln Foldseek/Multimer/MSA/msa --msa-format-mode 6 > log
foldseek unpackdb Foldseek/Multimer/MSA/msa Foldseek/Multimer/MSA/msa_output --unpack-suffix .a3m --unpack-name-mode 0 > log

rm log
