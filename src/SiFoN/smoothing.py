from scipy.signal import savgol_filter
import scipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def smooth_class(df, sequence_class_name, scores_colname="Raw Scores", window=15):
    smooth_df = df[df["Sequence Class"] == sequence_class_name]
    smooth_df["Smoothed Scores"] = smooth_df[scores_colname].rolling(window, center=True).mean().values
    # removes elements at the edge of the window  that cannot be avg'd
    smooth_df = smooth_df[smooth_df["Smoothed Scores"].notna()] 
    return smooth_df

def plot_smooth_v_raw(df, figname, sequence_class_name, s=0.2):
    fig, (ax_raw, ax_smooth) = plt.subplots(1, 2, figsize=(20, 6), sharex=True)
    ax_raw.scatter(df["Position"].values, df["Raw Scores"].values, s=s)
    ax_smooth.scatter(df["Position"].values, df["Smoothed Scores"].values, s=s)
    
    ax_raw.set_ylabel("Raw Scores")
    ax_smooth.set_ylabel("Smoothed Scores")
    fig.text(0.5, 0.00, 'Position', ha='center', va='center')
        
    plt.subplots_adjust(wspace=0.2)
    plt.rcParams.update({'font.size': 30})
    plt.savefig(figname)
    plt.show()
    
def smooth_all_scores(scores, window_length = 801, polyorder = 4, sav = True, sm = True):
    if sav: smoothed_scores = scores.apply(lambda x: savgol_filter(x, window_length, polyorder))
    if sm:  smoothed_scores.apply(lambda y: scipy.special.softmax(y), axis=1) 
    return  smoothed_scores