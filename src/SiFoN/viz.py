import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objs import *
import plotly.graph_objects as go

# TODO: add seqclass_names to github and directly link
seqclass_names = pd.read_csv("../model_data/seqclass-names.txt", header=None, sep="\n").to_numpy()
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

def get_category(index):
    for category in class_map.keys():
        if index in class_map[category]: return category

def post_processing(scores_max, seqclass_names=seqclass_names):
    scores_max.columns = ["Max Score", "Sequence Name"]
    scores_max["Sequence Index"] = [np.where(seqclass_names == name)[0][0] for name in scores_max["Sequence Name"]]
    scores_max["Class"] =  [get_category(index) for index in scores_max["Sequence Index"]]
    scores_max["Color"] =  [color_map[index] for index in scores_max["Class"]]
    scores_max["Function"] = ["Repressed" if cl in ["PC", "HET"] else "Active" for cl in scores_max["Class"]]
    scores_max = scores_max[scores_max["Class"] != "L"]
    scores_max = scores_max.reset_index()
    return scores_max

# If you run `preprocess_unsigned`, you should run `find_max_unsigned`.
def preprocess_unsigned(file, vcf, seqclass_names=seqclass_names):
    scores = pd.DataFrame(np.load(file), columns=seqclass_names)
    scores.index = vcf["Position"]
    scores = scores.abs()
    scores.sort_index(kind="mergesort", inplace=True) # sort by position
    scores = scores.apply(lambda x: x.rolling(window=3).max()) # find maximum prediction from in silico mut
    scores = scores.iloc[::3, :]
    scores = scores.iloc[1:,:]
    scores.columns = [str(col)[2:-3] for col in scores.columns]
    return scores

# If you run `preprocess_signed`, you should run `find_max_signed`.
def preprocess_signed(file, vcf, seqclass_names=seqclass_names):
    scores = pd.DataFrame(np.load(file), columns=seqclass_names)
    scores.index = vcf["Position"]
    scores.sort_index(kind="mergesort", inplace=True) # sort by position
    scores_min, scores_max = scores.apply(lambda x: x.rolling(window=3).min()), scores.apply(lambda x: x.rolling(window=3).max()) 
    scores_min, scores_max = scores_min.iloc[::3, :], scores_max.iloc[::3, :]
    scores_min, scores_max = scores_min.iloc[1:,:], scores_max.iloc[1:,:]
    scores = scores_max.where(abs(scores_max) > abs(scores_min), scores_min)
    scores.columns = [str(col)[2:-3] for col in scores.columns]
    return scores

def preprocess(file, vcf, seqclass_names=seqclass_names, signed=True):
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

def find_max(scores, signed=True):
    if signed:
        scores_max, scores_min = scores.max(axis = 1).to_frame(), scores.min(axis = 1).to_frame()
        max_idx, min_idx = scores.idxmax(axis = 1).to_frame(), scores.idxmin(axis = 1).to_frame()
        scores_max = pd.concat((scores_max.where(abs(scores_max) > abs(scores_min), scores_min), 
                        max_idx.where(abs(scores_max) > abs(scores_min), min_idx)), axis=1)
    else:
        scores_max = pd.concat((scores.max(axis = 1).to_frame(), scores.idxmax(axis = 1).to_frame()), axis=1)
        scores_max.columns = ["Max Score", "Sequence Name"]
    scores_max = post_processing(scores_max, seqclass_names=seqclass_names)
    return scores_max

def find_max_by_category(scores, signed=True):
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

def plot_max(filename, scores_prune, TSS):
    y, color = "Max Score", "Class"   
    fig = px.scatter(scores_prune, x="Position", y=y, color=color,
           hover_data=["Position", "Max Score", "Sequence Name", "Class"],
                    color_discrete_sequence = scores_prune["Color"].unique())
    TSS = {k: v for k, v in sorted(TSS.items(), key=lambda item: item[1])} 
    bot = min(scores_prune[y])
    for num, (tss_desc, tss_pos) in enumerate(TSS.items()):
        val = -num*2 - 2 + bot
        fig.add_trace(go.Scatter(y=[val, val], x=tss_pos, mode="lines+text", name=tss_desc, line=dict(width=7)))
    fig.update_layout(font=dict(size=18))
    fig.update_traces(marker={'size': 5, 'opacity':0.3})
    white_bg(fig)
    fig.write_html(filename)
    fig.show()
    
def plot_from_scores(filename, file, vcf, TSS, signed=True, plot_by="class", preprocess=False):
    if preprocess: scores = preprocess(file, vcf, signed=signed)
    if "category":
        max_scores = find_max_by_category(scores, signed=True)
    elif "class":
        max_scores = find_max(scores, signed=signed)
    else: return("Error: `plot_by` should be `class` or `category`")
    plot_max(filename, max_scores, TSS=TSS)

        