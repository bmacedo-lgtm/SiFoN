"""
This module makes it easy to combine clinical data with Sei sequence class predictions. In particular, this module includes functions to convert common clinical file formats to VCF files, a function to calculate odds ratios and assign SNPs to case or control labels based on these odds ratios, and a function to compare sequence class scores between cases and controls. 

Functions:
    cdot_to_VC
    
    odds_ratio
    
    add_case_control_label
    
    seq_class_t_tests
"""

import plotly.express as px
from scipy import stats
from SiFoN import viz
import pandas as pd
import numpy as np
import re

seqclass_names = pd.read_csv("../model_data/seqclass-names.txt", header=None, sep="\n").to_numpy()

def cdot_to_VCF(data, chrm, start=0):
    """Converts alterations from the c-dot (c. pos ref > alt e.g. c.300A>G) format to a VCF file.

    Parameters
    ----------
    data : Pandas DataFrame
        Clinical data. Must contain a column "Alteration" in c. notation that describes SNP's genomic coordinates.
    chrm : string
        Chromosome of alterations. Should be in the form "chr" + the chromosome number (e.g. chr10).
    start : int, optional
        Reference genomic position for alterations (e.g if all SNPs were in reference to position 1000 and you had a SNP c.300A>G than the true position would be 1300), default is 0.

    Returns
    -------
    Pandas DataFrame
        Input data stored in a VCF format dataframe.

    """
    output = pd.DataFrame(columns = ["#CHROM", "POS", "ID", "REF", "ALT"])
    output["#CHROM"] = [chrm for a in data["Alteration"]]
    output["POS"] = [start + int(re.search('-?[0-9]+', alt)[0]) for alt in data["Alteration"]]
    output["ID"] = data["Alteration"]
    output["REF"] = [re.search('(?<=[0-9])(.*)(?=>)', alt).group(0)[-1] for alt in data["Alteration"]]
    output["ALT"] = [re.search('>.*', alt).group(0)[-1] for alt in data["Alteration"]]
    return output



def odds_ratio(df, control_col, case_col, num_cases, num_controls, correction=0.01):
    """Calculates the odds ratio of case and control counts. Adds a new column to your dataframe `df` called "Odds Ratio".

    Parameters
    ----------
    df : Pandas DataFrame
        Contains case and control counts and any other metadata.
    control_col : string
        Name of DataFrame column with control counts.
    case_col : string
        Name of DataFrame column with case counts.
    num_cases : int
        Number of total cases.
    num_controls : int
        Number of total controls.
    correction : float, optional
        Added to all counts to prevent division by zero in the case of zero counts, default is 0.01.

    """
    df["Odds Ratio"] = [((case + correction)/num_cases)/((control + correction)/num_controls)
                         for case, control in zip(df[control_col], df[case_col])]
      
def add_case_control_label(df, case_cuttoff, control_cutoff):    
    """Assigns each SNP as a "Case", "Control", or "Equal" based on odds ratio scores. Adds dataframe column "Case/Control" to `df`.

    Parameters
    ----------
    df : Pandas DataFrame
        Contains case and control counts and any other metadata. Must have an "Odds Ratio" column.
    case_cuttoff : int
        Cutoff odds ratio to consider something a "Case" SNP.
    control_cutoff : int
        Cutoff odds ratio to consider something a "Control" SNP.
    """
    df["Case/Control"] = ["Case" if OR > case_cuttoff else "Control" if OR < control_cutoff else "Equal"
                              for OR in df["Odds Ratio"]]

def seq_class_t_tests(df, figname, fontsize=18, markersize=7):
    """Runs t-tests to compare sequence class scores in case and control populations. Returns p-values and then plots ranked p-values.

    Parameters
    ----------
    df : Pandas DataFrame
        Contains case and control counts and any other metadata. Must have scores for all 40 Sei sequence classes. Column names must correspond to sequence class names.
    figname : string
        Name of file that figure will be saved as.
    fontsize : int, optional.
        Fontsize for graph, defaults to 18,
    markersize : int, optional.
        Markersize for graph, defaults to 7,

    Returns
    -------
    Pandas DataFrame
        P-values for all 40 sequence class, ordered by ascending p-value.

    """
    pvals = np.zeros(len(seqclass_names))
    for ind, seq_class in enumerate(seqclass_names):
        pval = stats.ttest_ind(df[df["Case/Control"] == "Case"][seq_class], 
                    df[df["Case/Control"] == "Control"][seq_class]).pvalue
        pvals[ind] = pval   
    pvals_df = pd.DataFrame({"Sequence Class" : seqclass_names.flatten(), "P-value" : pvals})
    pvals_df.sort_values(by="P-value", ascending=True, inplace=True)
    pvals_df["Class"] =  [viz.get_category(index) for index in pvals_df["Sequence Class"]]
    pvals_df["Color"] =  [viz.color_map[index] for index in pvals_df["Class"]]
    pvals_df["Rank"] = range(0, len(pvals_df))
    pvals_df.head(5)
    
    fig = px.scatter(pvals_df, x="Rank", y="P-value", color="Class", hover_data=pvals_df.columns,
              color_discrete_sequence = pvals_df["Color"].unique())
    fig.update_layout(font=dict(size=fontsize), yaxis_title = "T-test P-values")
    fig.update_traces(marker={'size': markersize})
    viz.white_bg(fig)
    fig.write_html(figname)
    fig.show("png")
    return pvals_df
