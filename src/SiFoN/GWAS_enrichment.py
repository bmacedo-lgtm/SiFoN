"""
This module helps prioritize risk loci based on sequence class enrichment.

Functions:
    find_and_plot_enrichment
"""
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import seaborn as sns
import pyranges as pr
import pandas as pd
import numpy as np

def _sort_intersection(intersect_ref):
    intersect_ref.sort_values(by="End", inplace=True, kind='mergesort')
    intersect_ref.sort_values(by="Start", inplace=True, kind='mergesort' )
    intersect_ref.sort_values(by="Chromosome", inplace=True, kind='mergesort')
    return intersect_ref
    
def find_and_plot_enrichment(ccv, ref_bed, figname, vmin=0.1, vmax=10^1):
    """Calculates sequence class enrichment in genomic sequences compared to the entire genome and then creates a heatmap to visualize results.

    Parameters
    ----------
    ccv : Pandas DataFrame
        Pandas DataFrame in BED format containing disease risk loci. Must include at least the following columns: Chrom, Start, End, SNP (identifier/name of sequence).
    ref_bed : pyranges BED object
        Should contain sequence class labels across (nearly) the entire genome to allow for background sequence class proportion calculations. Should use the same genomic coordinate system (e.g. hg19 or hg38) as `ccv`.
    figname : string
        Name of file that figure will be saved as.
    vmin : float, optional
        Value corresponding to minimum color in heatmap, default is 0.1.
    vmax : float, optional
        Value corresponding to minimum color in heatmap, default is 10^1.

    Returns
    -------
    Pandas DataFrame
        Sequence class enrichment dataframe. Rows correspond to risk loci and the indices are the "SNP" column of `ccv` (the name/ID for each sequence). Columns correspond to the 40 sequence classes.

    """
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
    sns.set(rc = {'figure.figsize':(24,16)})
    cmap = sns.diverging_palette(255, 0, sep=8, n=256)
    g = sns.clustermap(counts_in_loci, cmap=cmap, 
                   vmin=vmin, vmax=vmax, xticklabels=1,  norm=LogNorm())
    g.ax_col_dendrogram.set_xlim([0,0.0001])
    g.ax_row_dendrogram.set_visible(False)
    g.savefig(figname)
    plt.show()
    
    return counts_in_loci
