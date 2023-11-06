
# Import the required packages
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text

def preprocess_volcano(df_volc):

    """
    This function takes a pd.DataFrame object and returns it with an extra column and 
    without missing values
    """
    
    try:
        # Create a negative log10 transformed column from the p adjusted values
        df_volc['nlog10'] = -np.log10(df_volc.padj)

        # Eliminate the missing or infinite values
        df_volc.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_volc.dropna(how='any', axis=0, inplace=True)

        return df_volc

    except KeyError:
        print("THERE IS A KEY ERROR, YOUR DATAFRAME MUST HAVE THE FOLLOWING COLUMNS:")
        print("1- 'gene_name': a column containing the gene symbols")
        print("2- 'log2FoldChange': a column containing the log2 fold changes values")
        print("3- 'padj': a column containing the adjusted p values")

        raise   


def map_color(a):
    """
    Accepts a tuple with a value of log2 fold change and negative log10 transformed
    p value and returns different strings depending on whether it is up|down regulated 
    or if its change is not significant
    """
    
    log2FoldChange, nlog10 = a
    
    if log2FoldChange > 1 and nlog10 > 2:
        return 'upregulated'
    elif log2FoldChange < -1 and nlog10 > 2:
        return 'downregulated'
    return 'not significant'

def create_texts(df_volc, texto):

    """
    This function takes a dataframe and returns a list with plt.text objects that contains significant regulated genes from a list of genes
    provided if it is provided or genes with the lower p adjusted values and the higher log2 fold changes esle
    """
    if type(texto)==list:
         texts = []
         for i in range(len(df_volc)):
            if df_volc["color"][i] != 'not significant' and df_volc["gene_name"][i] in texto:# esto es para que solo tenga que ser up o down
                texts.append(plt.text(x = df_volc.iloc[i].log2FoldChange, y = df_volc.iloc[i].nlog10, s = df_volc.iloc[i].gene_name,
                                    fontsize = 12, weight = 'bold'))
        

    else:
        texts = []
        for i in range(len(df_volc)):
            if df_volc.iloc[i].nlog10 > 25 and abs(df_volc.iloc[i].log2FoldChange) > 2:
                texts.append(plt.text(x = df_volc.iloc[i].log2FoldChange, y = df_volc.iloc[i].nlog10, s = df_volc.iloc[i].gene_name,
                                    fontsize = 12, weight = 'bold'))
                
    return texts


def volcano_plot(df_volc,path_save = False, texto =None):

    """
    This is the function that actually makes the plot, it needs a dataframe from a differenctial expression 
    analysis and create a vocano plot with the data in it, if a path is provided it is saved there and it can
    recive a list of genes to be marked in the volcano plot if they appear in the dataframe and has a signidicative
    change (padj<0.01 | log2FC>1)
    """

    df_volc = preprocess_volcano(df_volc)

    df_volc['color'] = df_volc[['log2FoldChange', 'nlog10']].apply(map_color, axis = 1)
    
    # Set figure size
    plt.figure(figsize = (6,6))

    # Plot graph
    ax = sns.scatterplot(data = df_volc, x = 'log2FoldChange', y = 'nlog10', hue = 'color', hue_order=['upregulated', 'downregulated', 'not significant'],
                        palette=['red', 'blue', 'grey'],alpha=0.5)

    # axis lines in the center of the graph
    ax.axhline(2, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(1, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-1, zorder = 0, c = 'k', lw = 2, ls = '--')

    # Gene names in the graph
    texts = create_texts(df_volc, texto)

    if len(texts)>0:
        # Arrows to the dots
        adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'))


    # Legend position
    plt.legend(loc = 1, bbox_to_anchor = (1.4,1), frameon = False, prop = {'weight':'bold'})

    # Axis spines
    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Axis ticks
    ax.tick_params(width = 2)

    plt.xticks(size = 12, weight = 'bold')
    plt.yticks(size = 12, weight = 'bold')

    # Axis labels
    plt.xlabel("$log_{2}$ fold change", size = 15)
    plt.ylabel("-$log_{10}$ FDR", size = 15)

    if path_save:
        # Save figure
        plt.savefig(path_save, dpi = 300, bbox_inches = 'tight', facecolor = 'white')

    # Plot figure
    plt.show()
