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

def odds_ratio(df, control_col, case_col, num_cases, num_controls, correction):
    df["Odds Ratio"] = [((case + correction)/num_cases)/((control + correction)/num_controls)
                         for case, control in zip(df[control_col], df[case_col])]
    
def add_case_control_label(df, case_cuttoff, control_cutoff):    
    df["Case/Control"] = ["Case" if OR > case_cuttoff else "Control" if OR < control_cutoff else "Equal"
                              for OR in df["Odds Ratio"]]

def seq_class_t_tests(df, filename, fontsize=18, markersize=7):
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
    fig.write_html("filename")
    fig.show()
    return pvals_df