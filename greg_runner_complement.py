import pandas as pd
from greg import main as greg_main

def main():
    
    gene = input("please provide the name of the gene for which you wish to analyze the promoter region:\n")
    
    with open(f"./{gene}_comp.xlsx", "w") as output_file:
        
        output_file.write("Name of Sequence,G4CR Sequence,NT Polymorphism Values,G4CR Start Index,G4CR End Index,G4CR Length in NTs,G-content Percentage,N-total,N-tandem,Tm-max,Tm-median,Tm-min,Full Promoter Sequence\n")
        
        # NOTE:
        # Keys of 'alignment' are the sequence names, will want to loop through 'names'
        # Values of 'alignment' are lists of the sequence nucleotides

        # INPUT: the index of a TSS for the gene of interest, and the chromosome number where it is located
        # OUTPUT: the Gtand and Gtot values for every aligned sequence for the specified promoter region.

        max_loop = 7
        max_bulge = 3
        min_temp = 50
        
        tss = int(input("please provide the value for the index of the tss of the desired gene:\n"))
        chrom = int(input ("please provide the number of the chromosome on which this gene is located:\n"))
        alignment = pd.read_pickle(f'/home/mcb/users/dlim63/conservation/data/seqDictPad_chr{chrom}.pkl')

        names = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOO', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']

        for name in names:
            # Outer layer: run GReg on each sequence (human and ancestrals) successively
            promoter_list_seq = alignment[name][(tss-1999):(tss+2000)]
            promoter_list_seq.reverse()
            
            promoter_list_complement = [] # codons go here
            index_seq = [] # indices go here
            
            for index, codon in enumerate(promoter_list_seq):
                if codon == "A":
                    promoter_list_complement.append("T")
                elif codon == "T":
                    promoter_list_complement.append("A")
                elif codon == "C":
                    promoter_list_complement.append("G")
                elif codon == "G":
                    promoter_list_complement.append("C")
                if codon == "A" or codon == "C" or codon == "G" or codon == "T":
                    index_seq.append(index)
            
            promoter_seq_complement = "".join(promoter_list_complement)
            
            # Inner layer: for a given sequence, analyze every G4CR in the promoter region of the desired gene using GReg.
            this_output = greg_main(promoter_seq_complement, index_seq, max_loop, max_bulge, min_temp, name)
            output_file.write(this_output)
            
        output_file.write("TSS index,Chromosome Number,Gene Name/Number\n")
        output_file.write(f"{tss},{chrom},{gene}\n")

main()
