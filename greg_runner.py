import pandas as pd
from greg import main as greg_main

def main():
    
    gene = input("please provide the name of the gene for which you wish to analyze the promoter region:\n")
    comp_input = input("Is this the complementary strand of the given gene's promoter region? Y/N\n")
    comp = True if comp_input.lower() == "y" else False
    filename = f"./{gene}.xlsx" if not comp else f"./{gene}_comp.xlsx"
    
    with open(filename, "w") as output_file:
        
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

        names = ['_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD', 'hg38', '_HP', '_HPG', '_HPGP', '_HPGPN', '_HPGPNRMPC', '_HPGPNRMPCCS', '_HPGPNRMPCCSO', '_HPGPNRMPCCSOT', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOO', 
                '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESC']
        
        human_promoter = ""
        ancestral_promoter = ""
        human_forward = 0
        human_backward = 0
        first_mammal_promoter_seq = ""
        
        for name in names:
            # Outer layer: run GReg on each sequence (human and ancestrals) successively
            sequence = alignment[name][(tss-50000):(tss+50000)]
            
            forward_index = 50000
            backward_index = 49999
            forward_counter = 0
            backward_counter = 0
            forward = []
            backward = []
            forward_max = human_forward if name == 'hg38' else 2000
            backward_max = human_backward if name == 'h38' else 2000
            # Standardize every sequence to have length of exactly 4000 nucleotides
            
            while True:
                forward_nt = sequence[forward_index]
                backward_nt = sequence[backward_index]
                
                if name == 'hg38' or forward_counter >= forward_max and backward_counter >= backward_max: # break case: correct number of nucleotides
                    break
                
                if backward_counter < backward_max and name == '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD':
                    human_backward += 1
                    backward.append(backward_nt)
                elif backward_counter < backward_max and (backward_nt == "A" or backward_nt == "C" or backward_nt == "G" or backward_nt == "T"): # increment the backward index
                    backward.append(backward_nt)
                    backward_index -= 1
                    backward_counter += 1
                elif backward_counter < backward_max:
                    backward_index -= 1
                 
                if forward_counter < forward_max and name == '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD':
                    human_forward += 1
                    forward.append(forward_nt)   
                elif forward_counter < forward_max and (forward_nt == "A" or forward_nt == "C" or forward_nt == "G" or forward_nt == "T"): # increment the forward index
                    forward.append(forward_nt)
                    forward_index += 1
                    forward_counter += 1
                elif forward_counter < forward_max:
                    forward_index += 1

            if name == 'hg38':
                human_promoter = sequence[(tss-human_backward):(tss+human_forward)]
                promoter_list_seq = list(sequence[(tss-2000):(tss+2000)])
            else:
                backward.reverse()
                promoter_list_seq = backward + forward
            
            if not comp:
                promoter_seq = "".join(promoter_list_seq)
            else:
                promoter_list_complement = []
                
                for nt in promoter_list_seq:
                    if nt == "A":
                        promoter_list_complement.append("T")
                    elif nt == "T":
                        promoter_list_complement.append("A")
                    elif nt == "C":
                        promoter_list_complement.append("G")
                    elif nt == "G":
                        promoter_list_complement.append("C")
                    else:
                        promoter_list_complement.append(nt)
                
                promoter_list_complement.reverse()
                promoter_seq = "".join(promoter_list_complement)
            
            if name == '_HPGPNRMPCCSOTSJMCMMRHCCOOOSVCTOPBOCECFCMAOLPPEMMESCLETCEOD':
                ancestral_promoter = promoter_seq
                new_promoter_seq = ""
                for i in range(len(promoter_seq)):
                    nt = promoter_seq[i]
                    if (nt == "A" or nt == "C" or nt == "G" or nt == "T"):
                        new_promoter_seq += nt
                first_mammal_promoter_seq = new_promoter_seq
            else:
                # Inner layer: for a given sequence, analyze every G4CR in the promoter region of the desired gene using GReg.
                this_output = greg_main(promoter_seq, max_loop, max_bulge, min_temp, name)
                output_file.write(this_output)


        this_output = greg_main(first_mammal_promoter_seq, max_loop, max_bulge, min_temp, name)
        output_file.write(this_output)
                
        output_file.write("TSS index,Chromosome Number,Gene Name/Number\n")
        output_file.write(f"{tss},{chrom},{gene}\n\n\n")
        print(f"Human forward counter: {human_forward}")
        print(f"Human backward counter: {human_backward}")
        print(ancestral_promoter)
        print()
        print(human_promoter)
        match_num = 0
        for i in range(len(human_promoter)):
            human_nt = human_promoter[i]
            ancestral_nt = ancestral_promoter[i]
            if human_nt == ancestral_nt:
                match_num += 1

        sim_gapless_percentage = match_num/4000
        output_file.write("Match rate in NTs (ungapped)\n")
        output_file.write(f"{sim_gapless_percentage}\n")
        

main()
