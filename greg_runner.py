import pandas as pd
from greg import main as greg_main

def main():
    
    gene = input("please provide the name of the gene for which you wish to analyze the promoter region:\n")
    
    with open(f"./{gene}.xlsx", "w") as output_file:
        
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
        chrom = input("please provide the number of the chromosome on which this gene is located:\n")
        alignment = pd.read_pickle(f'/home/mcb/users/dlim63/conservation/data/seqDictPad_chr{chrom}.pkl')

        names = ['hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOO', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD']

        for name in names:
            # Outer layer: run GReg on each sequence (human and ancestrals) successively
            sequence = alignment[name][(tss-50000):(tss+50000)]
            
            forward_index = 50000
            backward_index = 49999
            forward_counter = 0
            backward_counter = 0
            forward = []
            backward = []
            
            # Standardize every sequence to have length of exactly 4000 nucleotides
            while True:
                
                forward_nt = sequence[forward_index]
                backward_nt = sequence[backward_index]
                
                if forward_counter >= 2000 and backward_counter >= 2000: # break case: correct number of nucleotides
                    break
                
                if backward_counter < 2000 and (backward_nt == "A" or backward_nt == "C" or backward_nt == "G" or backward_nt == "T"): # increment the backward index
                    backward.append(backward_nt)
                    backward_index -= 1
                    backward_counter += 1
                elif backward_counter < 2000:
                    backward_index -= 1
                    
                if forward_counter < 2000 and (forward_nt == "A" or forward_nt == "C" or forward_nt == "G" or forward_nt == "T"): # increment the forward index
                    forward.append(forward_nt)
                    forward_index += 1
                    forward_counter += 1
                elif forward_counter < 2000:
                    forward_index += 1
            
            backward.reverse()
            promoter_list_seq = backward + forward
            promoter_seq = "".join(promoter_list_seq)
            
            
            # Inner layer: for a given sequence, analyze every G4CR in the promoter region of the desired gene using GReg.
            this_output = greg_main(promoter_seq, max_loop, max_bulge, min_temp, name)
            output_file.write(this_output)
            
        output_file.write("TSS index,Chromosome Number,Gene Name/Number\n")
        output_file.write(f"{tss},{chrom},{gene}\n")

main()
