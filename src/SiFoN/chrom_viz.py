import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm

### Heatmap
def white_bg(fig):
    fig.update_layout({"plot_bgcolor": "rgba(0, 0, 0, 0)",
                    "paper_bgcolor": "rgba(0, 0, 0, 0)"})
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', gridcolor="#DCDCDC")
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', gridcolor="#DCDCDC")
    
def preprocess_diff(diff, row_labels):   
    diff_avg = np.max(diff, axis=1) # take the max score across all the profiles
    diff_avg = diff_avg.reshape(int(diff_avg.size/3), 3) # reshape -> columns represent different alts.
    indices = [np.argmax(mut) for mut in diff_avg] # select alt with highest average score. 
    starts = [i for i in range(0, diff.shape[0] - 2, 3)]
    indices = [a+b for a,b in zip(starts, indices)] 
    return diff[indices, :], row_labels.iloc[indices, :]

def rank_scatter_plot(diff, row_labels, loc, loc_index, figname, fontsize=18):
    chrom_names = pd.read_csv("../chromatin_plot_library/target_names.txt", header=None, sep="\n").to_numpy()
    chrom_names = [name[0].split("|") for name in chrom_names]
    name_df = pd.DataFrame(chrom_names, columns=["Tissue", "Class", "ID", "none"])
    df = pd.DataFrame(diff[loc_index, :], columns=["Diff"])
    df["Abs Diff"] = df["Diff"].abs()
    df = df.join(name_df)
    df.sort_values(by="Diff", inplace=True)
    df.reset_index(inplace=True)
    title_info = row_labels[row_labels["pos"] == loc] 
    title = title_info["chrom"].values[0] + ":" + str(title_info["pos"].values[0]) + " " + title_info["ref"].values[0] + ">" + title_info["alt"].values[0]
    data = df["Diff"]
    center = - min(data) / (max(data) - min(data)) ## these lines make zero values white
    colorscale = [[0, 'rgba(0, 0, 150, 0.85)'],   
                   [center, '#D3D3D3'],  
                   [1, 'rgba(150, 0, 0, 0.85)']]
    fig = px.scatter(df, x=df.index, y="Diff", color_continuous_scale=colorscale, 
           hover_data=[df.index, "Diff", "Class", "Tissue"], color="Diff", title=title)

    fig.update_layout(font=dict(size=fontsize), xaxis_title="Rank")
    white_bg(fig)
    fig.write_html(figname)
    fig.show()

"""TODO: define the expected chromosome names.
`data`: `row_labels`:
`loc`: the chromosomal position of the SNP of interest.
`loc_index`: the index of the SNP of interest in the `data` array.
`top_X`: the number of top scoring chromatin profiles to plot for each SNP. Increasing top_X increases the height of
plot. There will be `top_X` * `pos_window` rows and pos_window columns.
`window`: the number of to the right and left of `loc` that will be plotted. 
`fontsize`: fontsize for graph."""
def chromatin_profile_heatmap(data, row_labels, loc, loc_index, figname, top_X=5, pos_window=5, fontsize=8):
    chrom_names = pd.read_csv("../chromatin_plot_library/target_names.txt", header=None, sep="\n").to_numpy()
    chrom_names = [name[0].split("|") for name in chrom_names]
    heat_data = data[loc_index - pos_window : loc_index + pos_window, :]
    top_vals = np.apply_along_axis(lambda x: np.abs(x).argsort()[-top_X:], 1, heat_data)
    top_chromatin = top_vals.flatten()
    heat_data = heat_data[:, top_chromatin]

    xticklabels = [i for i in range(loc - pos_window, loc + pos_window)]
    title_info = row_labels[row_labels["pos"] == loc] # information about SNPs
    title = title_info["chrom"].values[0] + " " + str(min(xticklabels)) + ":" + str(max(xticklabels))
    yticklabels = [id[1] + " | " + id[0] for id in np.asarray(chrom_names)[top_chromatin]]
    yticklabels = [a + " | " + str(b)[-2:] for a, b in zip(yticklabels, np.repeat(xticklabels, top_X))]
    a = [row_labels[row_labels["pos"] == pos] for pos in xticklabels]
    xticklabels = [str(title_info["pos"].values[0]) + " " + title_info["ref"].values[0] + ">" + title_info["alt"].values[0]
                   for title_info in a]
    bounds = min(np.abs(min(heat_data.flatten())), max(heat_data.flatten()))
    cmap = sns.diverging_palette(255, 0, sep=2, n=256)

    g = sns.heatmap(data=heat_data.T, cmap=cmap, vmin=-bounds, vmax=bounds,
                yticklabels=yticklabels, xticklabels=xticklabels, cbar_kws={"shrink": 0.3})
    # sns.set(rc = {'figure.figsize':(6,17)}, font_scale=0.08)
    sns.set(rc = {'figure.figsize':(10,10)})
    g.set_xticklabels(g.get_xmajorticklabels(), fontsize = fontsize)
    g.set_yticklabels(g.get_ymajorticklabels(), fontsize = fontsize)
    g.set_title(title, fontsize = fontsize + 10)
    plt.xticks(rotation=90) 
    plt.tight_layout()
    plt.savefig(figname)
    plt.show()