import pandas as pd
import re

def cdot_to_VCF_format(data, chrm, start=0):
    data["POS"] = [start + int(re.search('-?[0-9]+', alt)[0]) for alt in data["Alteration"]]
    data["REF"] = [re.search('(?<=[0-9])(.*)(?=>)', alt).group(0)[-1] for alt in data["Alteration"]]
    data["ALT"] = [re.search('>.*', alt).group(0)[-1] for alt in data["Alteration"]]
    data["#CHROM"] = [chrm for a in data["Alteration"]]
    data.rename(columns={"Alteration":"ID"}, inplace=True)
    return data[["#CHROM", "POS", "ID", "REF", "ALT"]]

def to_VCF(data):
    ccvs_sm["#CHROM"] = [ re.search('[A-z]*[a-z]*[0-9]*_', reg).group()[:-1]
                                for reg in data["Alteration"].to_list()]
    ccvs_sm["POS"] = [ int(re.search('_[0-9]*_', reg).group()[1:-1])
                                    for reg in data["Alteration"].to_list()]

