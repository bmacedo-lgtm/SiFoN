from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import seaborn as sns
import pyranges as pr
import pandas as pd
import numpy as np
%matplotlib inline

def _sort_intersection(intersect_ref):
    intersect_ref.sort_values(by="End", inplace=True, kind='mergesort')
    intersect_ref.sort_values(by="Start", inplace=True, kind='mergesort' )
    intersect_ref.sort_values(by="Chromosome", inplace=True, kind='mergesort')
    return intersect_ref
    
def find_enrichment(ccv, ref_bed):
    intersect_bed = pr.PyRanges(ccv)
    intersect_ref1 = ref_bed.intersect(intersect_bed).as_df()
    intersect_ref2 = _sort_intersection(intersect_ref1)
    intersect_ref2 = intersect_bed.intersect(ref_bed).as_df()
    intersect_ref2 = _sort_intersection(intersect_ref2)
    seqclass_names = pd.read_csv("../model_data/seqclass-names.txt", header=None, sep="\n").to_numpy()
    seqclass_names = np.append(seqclass_names, 21*["NA"])
    intersect_ref2["Sequence Identity"] = intersect_ref1["Name"]
    intersect_ref2['Sei Sequence Annotation'] = [seqclass_names[name] for name in intersect_ref2["Sequence Identity"]]
    
    ref = ref_bed.as_df() # proportions of counts across the genome
    ref['Sei Sequence Annotation'] = [seqclass_names[name] for name in ref["Name"]]
    total_counts = ref['Sei Sequence Annotation'].value_counts().rename_axis('Sei Sequence Annotation').reset_index(name='counts')
    total_counts["counts"] = total_counts["counts"]/sum(total_counts["counts"])
    total_counts.sort_values(by='Sei Sequence Annotation', inplace=True)
    
    counts_in_loci = intersect_ref2.groupby(by=["SNP",'Sei Sequence Annotation']).size().unstack(fill_value=0)
    counts_in_loci.sort_index(axis=1, inplace=True)
    counts_in_loci = counts_in_loci.div(counts_in_loci.sum(axis=1), axis=0)  # calculate proportion of each label per region
    counts_in_loci = counts_in_loci.div(total_counts["counts"].values, axis=1) # normalize by proportion across the full genome
    counts_in_loci = counts_in_loci + 0.00000001
    return counts_in_loci

def plot_enrichment(counts_in_loci, figname, vmin=0.1, vmax=10^1):
    sns.set(rc = {'figure.figsize':(24,16)})
    cmap = sns.diverging_palette(255, 0, sep=8, n=256)
    g = sns.clustermap(counts_in_loci, cmap=cmap, 
                   vmin=vmin, vmax=vmax, xticklabels=1,  norm=LogNorm())
    g.ax_col_dendrogram.set_xlim([0,0])
    g.ax_row_dendrogram.set_visible(False)
    g.savefig(figname)
    plt.show()