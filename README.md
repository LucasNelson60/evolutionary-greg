## ACKNOWLEDGEMENTS: 
- GReg: Amos Zhang, Christopher Hennecker, and Anthony Mittermaier (https://github.com/Christopher-Hennecker/GReg)
- evoGReg, and modifications to GReg: Lucas Nelson
- Aligned sequence data: Dongjoon Lim and Mathieu Blanchette

Written as part of a joint research project between the McGill Computational Genomics Lab (under the supervision of Mathieu Blanchette) and Anthony Mittermaier's lab in the McGill Department of Chemistry.

Instructions:
1. Run copy_pickle.py to transfer alignment pickles to personal directory while deleting unused sequence alignments
2. Run evoGReg.py on a chrom-by-chrom basis to get initial file outputs
3. Run remove_redundant_files.py to remove files that are duplicates of existing files for a given gene and cause issues in the later simulation phase
4. Run fix_match_rate.py to update the match rates on each file to be the correct one

This data can then be used to generate G4-containing-region (G4CR) simulations to fit the real-life data into a distribution and calculate its significance.
