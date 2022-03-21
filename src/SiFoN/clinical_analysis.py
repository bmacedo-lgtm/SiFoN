import plotly.express as px
from scipy import stats
from SiFoN import viz
import pandas as pd
import numpy as np
import re

seqclass_names = pd.read_csv("../model_data/seqclass-names.txt", header=None, sep="\n").to_numpy()

def cdot_to_VCF(data, chrm, start=0):
    output = pd.DataFrame(columns = ["#CHROM", "POS", "ID", "REF", "ALT"])
    output["#CHROM"] = [chrm for a in data["Alteration"]]
    output["POS"] = [start + int(re.search('-?[0-9]+', alt)[0]) for alt in data["Alteration"]]
    output["ID"] = data["Alteration"]
    output["REF"] = [re.search('(?<=[0-9])(.*)(?=>)', alt).group(0)[-1] for alt in data["Alteration"]]
    output["ALT"] = [re.search('>.*', alt).group(0)[-1] for alt in data["Alteration"]]
    return output

"""Converts alterations from the c-dot (c. pos ref > alt) format to a VCF file.
:param data: Clinical data. Must contain a column "Alteration" in c. notation that describes SNP's genomic coordinates.
:type data: Pandas DataFrame
:param chrm: Chromosome of alterations. Should be in the form "chr" + the chromosome number (e.g. chr10).
:type chrm: string
:param start: Reference genomic position for alterations, default is 0
:type start: int, optional

:return: Input data stored in a VCF format dataframe.
:rtype: Pandas DataFrame
"""

def odds_ratio(df, control_col, case_col, num_cases, num_controls, correction):
    df["Odds Ratio"] = [((case + correction)/num_cases)/((control + correction)/num_controls)
                         for case, control in zip(df[control_col], df[case_col])]
    
"""Calculates the odds ratio of case and control counts. Adds a new column to your dataframe `df` called "Odds Ratio"
:param df: Contains case and control counts and any other metadata.
:type df: Pandas DataFrame
:param control_col: name of DataFrame column with control counts
:type control_col: string
:param case_col: name of DataFrame column with case counts
:type case_col: string
:param num_cases: number of total cases
:type num_cases: int
:param num_controls: number of total controls
:type num_controls: int
:param correction: added to all counts to prevent division by zero in the case of zero counts.
:type correction: float
"""
    
def add_case_control_label(df, case_cuttoff, control_cutoff):    
    df["Case/Control"] = ["Case" if OR > case_cuttoff else "Control" if OR < control_cutoff else "Equal"
                              for OR in df["Odds Ratio"]]

"""Assigns each SNP as a "Case", "Control", or "Equal" based on odds ratio scores. Adds dataframe column "Case/Control" to `df`.
:param df: Contains case and control counts and any other metadata. Must have an "Odds Ratio" column.
:type df: Pandas DataFrame
:param case_cuttoff: cutoff odds ratio to consider something a "Case" SNP
:type case_cuttoff: int
:param control_cutoff: cutoff odds ratio to consider something a "Control" SNP
:type control_cutoff: int
"""   


def seq_class_t_tests(df, figname, fontsize=18, markersize=7):
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

"""Runs t-tests to compare sequence class scores in case and control populations. Returns p-vals and plots ranked p-vals.

:param df: Contains case and control counts and any other metadata. Must have scores for all 40 Sei sequence classes. Column names must correspond to sequence class names. 
:type df: Pandas DataFrame
:param figname: name of file that figure will be saves as.
:type figname: string
:param fontsize: fontsize for graph, defaults to 18
:type fontsize: int, optional.
:param markersize: markersize for graph, defaults to 7
:type markersize: int, optional.

:return: P-values for all 40 sequence class, ordered by ascending p-value. 
:rtype: Pandas DataFrame
"""