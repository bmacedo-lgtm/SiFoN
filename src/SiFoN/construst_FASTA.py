from selene_sdk.sequences import Genome
import pandas as pd
import itertools
import random

def retrieve(ref_genome_fasta):
    genome = Genome(input_path = ref_genome_fasta)
    
def write_haps_of_pairs(data, ref_filename, alt_filename, shift_filename, centered_filename,
         chrm, cutoff):
    alt_file = open(alt_filename, "w")
    ref_file = open(ref_filename, "w")
    shift_file = open(shift_filename, "w") 
    centered_file = open(centered_filename, "w")
    chrm, offset, cutoff = 2048, chrm, cutoff
    for pairs, ID in zip(itertools.combinations(data["POS"].values, 2), 
                            itertools.combinations(data["ID"].values, 2)):

        if abs(pairs[0] - pairs[1] > cutoff): continue # skip if snps are too far apart
        if pairs[0] == pairs[1]: continue # skip if snps are at the same position
        desc.append(write_to_file(pairs[0], pairs[1], 
                     ref_file, alt_file, shift_file, centered_file, 
                                 ID[0], ID[1]))
        desc.append(write_to_file(pairs[1], pairs[0], 
                     ref_file, alt_file, shift_file, centered_file, 
                                 ID[1], ID[0]))  
    ref_file.close()
    alt_file.close()
    shift_file.close()
    centered_file.close()
    
    
def writes_haps_around_SNP(data, ref_filename, alt_filename, shift_filename, centered_filename,
         chrm, rng, center_SNP):
    alt_file = open(alt_filename, "w")
    ref_file = open(ref_filename, "w")
    shift_file = open(shift_filename, "w") 
    centered_file = open(centered_filename, "w")
    offset, chrm, rng, center_SNP = 2048, chrm, rng, center_SNP
    for pairs in itertools.combinations(range(-rng, rng), 2):
        write_to_file(center + pairs[0], center + pairs[1], 
                     ref_file, alt_file, shift_file, centered_file)

        write_to_file(center + pairs[1], center + pairs[0], 
                     ref_file, alt_file, shift_file, centered_file, 
                      center + pairs[1], center + pairs[0])  
    ref_file.close()
    alt_file.close()
    shift_file.close()
    centered_file.close()


def write_to_file(snp1, snp2, 
                 ref_file, alt_file, shift_file, centered_file, 
                 id1, id2):
    seq = genome.get_sequence_from_coords(chrom = chrm,
                                 start = snp1 - offset,
                                 end = snp1 + offset)
    # add in reference haplotypes
    preface = str(str(snp1) +  "_" + str(snp2) + "_" + id1 + "_" + id2)
    
    ref_file.write(">" + str(preface) + "\n" + seq + "\n")
    seq_list = list(seq)
    
    # add shifted reference 
    seq_list[offset + snp1 - snp2] = ID_to_alt_dict[id2]
    mut_seq = "".join(seq_list)
    shift_file.write(">" + str(preface) + "\n" + mut_seq + "\n")
    
    # add centered reference
    seq_list2 = list(seq)
    seq_list2[offset - 1] = ID_to_alt_dict[id1]
    mut_seq = "".join(seq_list2)
    centered_file.write(">" + str(preface) + "\n" + mut_seq + "\n")
    
    # add haplotype
    seq_list[offset - 1] = ID_to_alt_dict[id1]
    mut_seq = "".join(seq_list)
    alt_file.write(">" + str(preface) + "\n" + mut_seq + "\n")
    return preface

def write_to_file_smooth(snp1, snp2, ref_file, alt_file, window=16):
    seq = genome.get_sequence_from_coords(chrom = chrm,
                                 start = snp1 - offset,
                                 end = snp1 + offset)

    for i in range(-window, window):
        ref_file.write(">" + str(snp1) +  "_" + str(snp2) + "_" + str(i + window) + "\n" + seq + "\n")
        seq_list2 = list(seq)
        seq_list2[offset - 1 + i] = alt_allele[seq_list2[offset - 1 + i].upper().upper()]
        mut_seq = "".join(seq_list2)
        alt_file.write(">" + str(snp1) +  "_" + str(snp2) + "_" + str(i + window) + "\n" + mut_seq + "\n")
