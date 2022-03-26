"""
This module includes functions for visualizing sequence class scores along genomic regions. This module require Sei sequence class scores as input, and therefore assumes that predictions have already been calculated. 

Functions:

    white_bg
    
    preprocess
    
    postprocess
    
    find_max
    
    find_max_by_category
    
    plot_sequence_class
    
    plot_max
    
    plot_max_from_scores

Misc variables:

    seqclass
    color_map
    class_map

"""

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objs import *
import plotly.graph_objects as go

# TODO: add seqclass_names to github and directly link
seqclass = pd.read_csv("../model_data/seqclass-names.txt", header=None, sep="\n").to_numpy()
color_map = {
    'E': "#984ef3",
    'CTCF': "#fdb462",
    'P': "#ef3b2c",
    'PC': "#8abad4",
    'HET': "#662506",
    'TN': "#fb9a99",
    'TF': "#386cb0",
    'L': "#C0C0C0",
    'NA': "#000000"}

class_map = {'E': [12, 16, 36, 38, 5, 30, 7, 26, 6, 9, 17, 13, 
                  "E1 Stem cell", "E2 Multi-tissue", "E3 Brain / Melanocyte", "E4 Multi-tissue",
                   "E5 B-cell-like", "E6 Weak epithelial", "E7 Monocyte / Macrophage", "E8 Weak multi-tissue",
                   "E9 Liver / Intestine", "E10 Brain", "E11 T-cell", "E12 Erythroblast-like"],
             'CTCF': [27, "CTCF CTCF-Cohesin"],
             'P': [25, "P Promoter"],
             'PC': [34, 0, 15, 20, "PC1 Polycomb / Heterochromatin", 
                    "PC2 Weak Polycomb", "PC3 Polycomb", "PC4 Polycomb / Bivalent stem cell Enh"],
             'HET': [11, 23, 29, 32, 35, 39, 
                     "HET1 Heterochromatin", "HET2 Heterochromatin", "HET3 Heterochromatin", "HET4 Heterochromatin",
                    "HET5 Centromere", "HET6 Centromere"],
             'TN': [2, 3, 21, 28, 
                    "TN1 Transcription","TN2 Transcription","TN3 Transcription","TN4 Transcription"],
             'TF': [37, 19, 14, 31, 10,
                   "TF1 NANOG / FOXA1", "TF2 CEBPB", "TF3 FOXA1 / AR / ESR1", "TF4 OTX2", "TF5 AR"],
             'L': [8, 18, 22, 24, 33, 1, 4,
                   "L1 Low signal", "L2 Low signal", "L3 Low signal", "L4 Low signal",
                   "L5 Low signal", "L6 Low signal", "L7 Low signal"],
             'NA':[40]}

def white_bg(fig):
    fig.update_layout({"plot_bgcolor": "rgba(0, 0, 0, 0)",
                    "paper_bgcolor": "rgba(0, 0, 0, 0)"})
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', gridcolor="#DCDCDC")
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', gridcolor="#DCDCDC")

def _get_category(index):
    for category in class_map.keys():
        if index in class_map[category]: return category

def postprocess(scores_max, seqclass_names=[]):
    """Processes the `scores_max` array to have information necessary for graphing, including sequence class category and color.

    Parameters
    ----------
    scores_max :
        type scores_max: NumPy Array
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"

    Returns
    -------
    type
        an updated version of `scores_max` that adds the columns “Sequence Index”, “Class” (e..g E, P, TF), “Color” (corresponding to the “Class” column), “Function” (whether a sequence class is associated with repressed or active chromatin). SNPs with a max score in the Low Impact category are also removed.

    """
    if seqclass_names == []: seqclass_names=seqclass
    scores_max.columns = ["Max Score", "Sequence Name"]
    scores_max["Sequence Index"] = [np.where(seqclass_names == name)[0][0] for name in scores_max["Sequence Name"]]
    scores_max["Class"] =  [_get_category(index) for index in scores_max["Sequence Index"]]
    scores_max["Color"] =  [color_map[index] for index in scores_max["Class"]]
    scores_max["Function"] = ["Repressed" if cl in ["PC", "HET"] else "Active" for cl in scores_max["Class"]]
    scores_max = scores_max[scores_max["Class"] != "L"]
    scores_max = scores_max.reset_index()
    return scores_max

def preprocess(file, vcf, seqclass_names=[], signed=True):
    """Processes data so that there is a single row per position representing the alteration with the maximum average score at that position.

    Parameters
    ----------
    file : string
        filename of NumPy sequence class scores
    vcf : Pandas DataFrame
        VCF of SNPs of interest, including “Position” column.
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"
    signed : Boolean, optional
        signifies whether parameters are evaluated by absolute value (False) or not (True), default is True

    Returns
    -------
    NumPy Array
        array of scores, including a single row per position representing the alteration with the maximum average score

    """
    if seqclass_names == []: seqclass_names=seqclass
    scores = pd.DataFrame(np.load(file), columns=seqclass_names)
    scores.index = vcf["Position"]
    if signed:
        scores.sort_index(kind="mergesort", inplace=True) # sort by position
        scores_min, scores_max = scores.apply(lambda x: x.rolling(window=3).min()), scores.apply(lambda x: x.rolling(window=3).max()) 
        scores_min, scores_max = scores_min.iloc[::3, :], scores_max.iloc[::3, :]
        scores_min, scores_max = scores_min.iloc[1:,:], scores_max.iloc[1:,:]
        scores = scores_max.where(abs(scores_max) > abs(scores_min), scores_min)
    else:
        scores = scores.abs()
        scores.sort_index(kind="mergesort", inplace=True) # sort by position
        scores = scores.apply(lambda x: x.rolling(window=3).max()) # find maximum prediction from in silico mut
        scores = scores.iloc[::3, :]
        scores = scores.iloc[1:,:]
    scores.columns = [str(col)[2:-3] for col in scores.columns]
    return scores

def find_max(scores, signed=True, seqclass_names=[]):
    """Takes the (abs) max sequence class score for each position/alteration.

    Parameters
    ----------
    scores : NumPy array
        Takes the (abs) max sequence class score for each position/alteration. Should run the `preprocessing` function first. The output array has columns [“Max Score", "Sequence Name"] corresponding the (abs) max score and the corresponding sequence class. This array is then past through the `postprocess` function.
    signed : Boolean, optional
        signifies whether parameters are evaluated by absolute value (False) or not (True). This should be set to the same value as in the `preprocessing` function, default is True
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"

    Returns
    -------
    Pandas DataFrame
        Dataframe with one prediction per alteration, corresponding the maximum sequence class score and its corresponding sequence class name. Metadata used for graphing is also included in the dataframe.

    """
    if seqclass_names == []: seqclass_names=seqclass
    if signed:
        scores_max, scores_min = scores.max(axis = 1).to_frame(), scores.min(axis = 1).to_frame()
        max_idx, min_idx = scores.idxmax(axis = 1).to_frame(), scores.idxmin(axis = 1).to_frame()
        scores_max = pd.concat((scores_max.where(abs(scores_max) > abs(scores_min), scores_min), 
                        max_idx.where(abs(scores_max) > abs(scores_min), min_idx)), axis=1)
    else:
        scores_max = pd.concat((scores.max(axis = 1).to_frame(), scores.idxmax(axis = 1).to_frame()), axis=1)
        scores_max.columns = ["Max Score", "Sequence Name"]
    scores_max = postprocess(scores_max, seqclass_names=seqclass_names)
    return scores_max

def find_max_by_category(scores, signed=True, seqclass_names=[]):
    """Takes the (abs) max sequence class score within a category for each position/alteration.

    Parameters
    ----------
    scores : NumPy array
        Takes the (abs) max sequence class score within a category for each position/alteration. Should run the `preprocessing` function first. The output array has columns [“Max Score", "Sequence Name"] corresponding the (abs) max score and the corresponding sequence class. This array is then past through the `postprocess` function.
    signed : Boolean, optional
        signifies whether parameters are evaluated by absolute value (False) or not (True). This should be set to the same value as in the `preprocessing` function, default is True
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"

    Returns
    -------
    Pandas DataFrame
        Dataframe with seven predictions per alteration, corresponding to the maximum score within each of the seven sequence classes (excluding L and NA). Metadata used for graphing is also included in the dataframe.

    """
    if seqclass_names == []: seqclass_names=seqclass
    for key in color_map.keys():
        if key == "L" or key == "NA" : continue
        good_col = [col in class_map[key] for col in scores.columns]
        col_indices = [i for i, x in enumerate(good_col) if x]
        df = scores.iloc[:, col_indices]
        if key == "E":
            max_by_cat = find_max(df, signed)
        else:
            max_by_cat = pd.concat([max_by_cat, find_max(df, signed)])
    return max_by_cat

def plot_sequence_class(filename, file, vcf, TSS={}, category="All", signed=True, pre=False, seqclass_names=[]):
    """Preprocesses Sei prediction data and then plots the results. 

    Parameters
    ----------
    filename : string
        name of file that figure will be saved as.
    file : string
        filename of NumPy sequence class scores
    vcf : Pandas DataFrame
        VCF of SNPs of interest, including “Position” column.
    TSS : Dict[str, list[int, int]], optional
        Specifies regions of the genome to annotate in the figure. The first element of the dictionary (str) will be the name or label of the annotation which will appear in the legend. The list of ints will denote the end points of the region to be annotated, defaults to an empty dictionary
    category : string, optional
        Determines which sets of predictions will be plotted. Should be one of the category names listed in the `color_map`. If set to "All" then the full set of predictions is plotted, default is "All"
    signed : Boolean, optional
        signifies whether parameters are evaluated by absolute value (False) or not (True). This should be set to the same value as in the `preprocessing` function, default is True
    pre : Boolean, optional
        determines whether to run preprocessing or not. If there is already only one alteration per position, then this should be set to False, default is False.
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"
    """
    if seqclass_names == []: seqclass_names=seqclass
    if pre: 
        scores = preprocess(file, vcf, seqclass_names, signed=signed)
    else: 
        scores = pd.DataFrame(np.load(file), columns=seqclass_names)
        scores.index = vcf["Position"]
    scores = scores.melt(ignore_index=False)
    scores.rename(columns={"variable":"Sequence Name", "value":"Score"}, inplace=True)
    scores = scores[["Score", "Sequence Name"]]
    if category != "All":
        scores = scores[[_get_category(seq) == category for seq in scores["Sequence Name"]]]
    scores = postprocess(scores, seqclass_names=seqclass_names)
    plot_max(filename, scores, TSS, yaxis_title=category + " Scores")
    return scores
                
def plot_max(filename, scores_prune, TSS={}, yaxis_title="Max Score"):
    """Plots the maximum Sei classes along a genomic sequence.

    Parameters
    ----------
    filename : string
        name of file that figure will be saved as.
    scores_prune : Pandas DataFrame
        Dataframe that denotes the maximum Sei sequence classes along some region of the genome. Must contain the following columns: ["Max Score", "Class", "Position", "Sequence Name"]. This is generated as output from the `find_max` and `find_max_by_category` functions.
    TSS : Dict[str, list[int, int]]
        Specifies regions of the genome to annotate in the figure. The first element of the dictionary (str) will be the name or label of the annotation which will appear in the legend. The list of ints will denote the end points of the region to be annotated.
    """
    y, color = "Max Score", "Class"   
    fig = px.scatter(scores_prune, x="Position", y=y, color=color,
           hover_data=["Position", "Max Score", "Sequence Name", "Class"],
                    color_discrete_sequence = scores_prune["Color"].unique())
    TSS = {k: v for k, v in sorted(TSS.items(), key=lambda item: item[1])} 
    bot = min(scores_prune[y])
    for num, (tss_desc, tss_pos) in enumerate(TSS.items()):
        val = -num*2 - 2 + bot
        fig.add_trace(go.Scatter(y=[val, val], x=tss_pos, mode="lines+text", name=tss_desc, line=dict(width=7)))
    fig.update_layout(font=dict(size=18), yaxis_title=yaxis_title)
    fig.update_traces(marker={'size': 5, 'opacity':0.3})
    white_bg(fig)
    fig.write_html(filename)
    fig.show("png")
    
def plot_max_from_scores(filename, file, vcf, TSS={}, signed=True, plot_by="class", pre=False, seqclass_names=[]): 
    """Preprocesses Sei prediction data, finds the maximum, and then plots the results. Runs `preprocessing`, `find_max` or `find_max_by_category`, and then `plot_max`.

    Parameters
    ----------
    filename : string
        name of file that figure will be saved as.
    file : string
        filename of NumPy sequence class scores
    vcf : Pandas DataFrame
        VCF of SNPs of interest, including “Position” column.
    TSS : Dict[str, list[int, int]], optional
        Specifies regions of the genome to annotate in the figure. The first element of the dictionary (str) will be the name or label of the annotation which will appear in the legend. The list of ints will denote the end points of the region to be annotated, defaults to an empty dictionary
    signed : Boolean, optional
        signifies whether parameters are evaluated by absolute value (False) or not (True). This should be set to the same value as in the `preprocessing` function, default is True
    plot_by : string, optional
        Specifies whether the maximum of sequence class is plotted ("class") or the maximim within each sequence category ("category). Should be either "class" (calls `find_max`) or "category" (calls `find_max_by_category`), defaults to "class"
    pre : Boolean, optional
        determines whether to run preprocessing or not. If there is already only one alteration per position, then this should be set to False, default is False.
    seqclass_names : list of strings, optional
        List of names for profiles. Sei sequence classes names or chromatin profile names. Should have length equal to the number of rows in`scores_max`, defaults to Sei sequence class names found in "model_data/seqclass-names.txt"
    """
    if seqclass_names == []: seqclass_names=seqclass
    if pre: 
        scores = preprocess(file, vcf, seqclass_names, signed=signed)
    else: 
        scores = pd.DataFrame(np.load(file), columns=seqclass_names)
        scores.index = vcf["Position"]
    if "category":
        max_scores = find_max_by_category(scores, signed=True)
    elif "class":
        max_scores = find_max(scores, signed=signed)
    else: return("Error: `plot_by` should be `class` or `category`")
    plot_max(filename, max_scores, TSS=TSS)

        