from scipy.signal import savgol_filter
import scipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def smooth_class(df, sequence_class_name, window=15):
    """Applies a rolling average on sequence claas scores for a particular class. 

    Parameters
    ----------
    df : Pandas DataFrame
        dataframe with the columns "Raw Scores", "Sequence Class", and "Position". "Raw scores" are Sei sequence class scores and "Sequence Class" denotes the corresponding sequence class of the score. It is assumed that all SNPs are along the same chromosome, allowing for smoothing via adjacency. "Position" denotes the position of the SNP on the genome. 
    sequence_class_name : string
        the name of the sequence class to apply rolling average to.
    window : int, optional
        window of rolling average, default is 15.
    Returns
    -------
    Pandas DataFrame
       Dataframe containing scores for the `sequence_class_name` class after the application of a rolling average. Note that elements on the window edge are removed from the dataframe such that the output will have a smaller size than the input dataframe.
"""
    smooth_df = df[df["Sequence Class"] == sequence_class_name]
    smooth_df["Smoothed Scores"] = smooth_df["Raw Scores"].rolling(window, center=True).mean().values
    # removes elements at the edge of the window  that cannot be avg'd
    smooth_df = smooth_df[smooth_df["Smoothed Scores"].notna()] 
    return smooth_df

def plot_smooth_v_raw(df, figname, sequence_class_name, s=0.2):
    """Plots smooth and raw sequence class scores side by side to compare visualizations.

    Parameters
    ----------
    df : Pandas DataFrame
        dataframe with the columns "Raw Scores", "Smoothed Scores", and "Position". 
    figname : string
        name of file that figure will be saved as.
    sequence_class_name : string
        the name of the sequence class to apply rolling average to.
    s : float, optional
        size of data points, default is 0.2
"""
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
    
def smooth_all_scores(df, window_length = 801, polyorder = 4, sav = True, sm = True):
    """Plots smooth and raw sequence class scores side by side to compare visualizations.

    Parameters
    ----------
    df : Pandas DataFrame
        dataframe with the columns "Raw Scores", "Smoothed Scores", and "Position". 
    window_length : int, optional
        the length of the Savitzy-Golay filter window, default is 801
    polyorder : int, optional
        the order of the polynomial used to fit the samples in Savitzy-Golay filter, default is 4
    sav : Boolean
        whether or not to apply Savitzy-Golay filtering, default is True
    sm : Boolean
        whether or not to apply Softmax, default is True
"""
    scores = df
    if sav: scores = scores.apply(lambda x: savgol_filter(x, window_length, polyorder))
    if sm:  scores.apply(lambda y: scipy.special.softmax(y), axis=1) 
    return  scores