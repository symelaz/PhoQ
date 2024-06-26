This script runs a Principle Component Analysis between HAMP structures:

- Target: HAMP structures from [Aleksander Winski, Jan Ludwiczak, Malgorzata Orlowska, Rafal Madaj, Kamil Kaminski and Stanislaw Dunin-Horkawicz, AlphaFold2 captures the conformational landscape of the HAMP signaling domain, Protein Sci. 33, e4846 (2024), doi: 10.1002/pro.4846](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.4846)
- Query: HAMP structure of the AlphaFold predicted PhoQ structure (struct_in/phoq_full_af_HAMP.pdb)

## Prerequisites
Before running the script ensure that you have installed the [Foldseek](https://github.com/steineggerlab/foldseek) tool

## Steps
1. Convert the HAMP bundle into a single chain (to be considered by Foldseek as one chain)
2. Run Foldseek analysis to get the structural alignments and MSAs of the structures
3. Run pca analysis 
