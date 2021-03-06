""" This module includes functions that create FASTA files containing haplotypes as input into Sei"""

from selene_sdk.sequences import Genome
import pandas as pd
import itertools
import random

alt_allele = {"A":"G", "G":"A", "T":"C", "C":"T"}
    
def write_haps_of_pairs(data, ref_filename, alt_filename, shift_filename, centered_filename,
         chrm, cutoff, gn):
    """Give a VCF file of SNPs, this module will create a FASTA file containing haplotype sequences that consist of all pairs of SNPs within a `cutoff` distance from one another. These FASTA files are in a format that can be read by Sei.

    Parameters
    ----------
    data : Pandas DataFrame 
        NumPy array of Sei chromatin profile scores. Should have shape (X, 21907), where 21,907 corresponds to the number of 
    ref_filename : string
        Filename for referance file
    alt_filename : string
        Filename for alternative file (containing all haplotype SNPs)
    shift_filename : string
        Filename for shifted file (only contains the off center SNP)
    centered_filename : string
        Filename for centered file (only contains the centered SNP)
    chrm : string
        Chromsome of haplotype. Should be in the format "chr" + number (e.g. chrm10).
    cutoff : int
        Specifies how far away two SNPs need to be to calculate a haplotype. The default is 2000. Since Sei works with 4kb input strings, this value should be <= 2000. 
    gn : selene_sdk.sequences.Genome
        Reference genome
    """
    alt_file = open(alt_filename, "w")
    ref_file = open(ref_filename, "w")
    shift_file = open(shift_filename, "w") 
    centered_file = open(centered_filename, "w")
    offset = 2048
    ID_to_alt_dict = dict(zip(data["ID"].values, data["ALT"].values))
    for pairs, ID in zip(itertools.combinations(data["POS"].values, 2), 
                            itertools.combinations(data["ID"].values, 2)):
        
        if abs(pairs[0] - pairs[1] > cutoff): continue # skip if snps are too far apart
        if pairs[0] == pairs[1]: continue # skip if snps are at the same position
        write_to_file(pairs[0], pairs[1], chrm,
                     ref_file, alt_file, shift_file, centered_file, 
                                 ID[0], ID[1], ID_to_alt_dict, gn)
        write_to_file(pairs[1], pairs[0], chrm,
                     ref_file, alt_file, shift_file, centered_file, 
                                 ID[1], ID[0], ID_to_alt_dict, gn)
    ref_file.close()
    alt_file.close()
    shift_file.close()
    centered_file.close()   

def write_to_file(snp1, snp2, chrm,
                 ref_file, alt_file, shift_file, centered_file, 
                 id1, id2, ID_to_alt_dict, gn):
    offset = 2048
    seq = gn.get_sequence_from_coords(chrom = chrm,
                                 start = snp1 - offset,
                                 end = snp1 + offset)
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