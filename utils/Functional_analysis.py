
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import gseapy as gp
from gseapy import barplot, dotplot


def generate_cls(condition,cls_path="./Data/gsea.cls"):

    """
    Parameters: 
    
    - condition: sorted  np.array or pd.Series with the condition that will be tested for each sample. 
                 Only two different conditions allowed.

    - cls_path: (str) path to the new cls file that will be created.

    ================================================================================

    Returns: 

    - cls_path: (str) path to the cls file that has been created.
 
    """
    if np.all(np.sort(condition)[::-1]==condition) or np.all(np.sort(condition)[::-1]==condition):
            conditions = np.unique(condition)
            # Creo las dos condiciones
            cls_list = "\t".join([ "1" if cond == conditions[0] else "0" for cond in condition])
            text= [str(len(condition)),
                    " 2", " 1\n", 
                    f"# {conditions[0]} {conditions[1]}\n"]

            with open(cls_path,"w") as f:
                    f.writelines(text)
                    f.writelines(cls_list)
            print(f"cls file succesfully created in {cls_path}")
            return cls_path
    else:
        print("This isnt a sorted dataset, please sort it and re run the function.\nTake into account that the first condition will be the tested.")
        

def sig_finder(gs_res,threshold=0.05):

    """ 
    DEFINO UNA FUNCION QUE ES CAPAZ DE CREAR UN DICCIONARIO Y UN DATAFRAME CON LOS RESULTADOS SIGNIFICATVOS DE UN GSEA, 
    HAY QUE PASARLE EL THRESHOLD DE FDR PVALUE DESEADO (DEFAULT = 0.05) Y UN OBJETO CON EL RESULTADO DE CORRER EL GSEA 
    """

    sig_results = {}

    name,Term, es , nes ,pval , fdr,fwerp,tag ,gene ,lead_genes, matched_genes,hits,RES = [],[],[],[],[],[],[],[],[],[],[],[],[]
    atts = (name, es , nes ,pval , fdr,fwerp,tag ,gene ,lead_genes, matched_genes,hits,RES) 

    for i in gs_res.results.keys():
        result = gs_res.results[i]
        if result['fdr'] < threshold:
            key = i.split(".")[0]
            sig_results[key] = result
            n=0
            for j in result.keys():
                atts[n].append(result[j])
                n+=1
            Term.append(key)
            
    sigres_df = {'name': name,
                'Term':Term,
                'es':es,
                'nes':nes,
                'pval':pval,
                'fdr':fdr,
                'fwerp':fwerp,
                'tag %':tag,
                'gene %':gene,
                'lead_genes':lead_genes,
                'matched_genes':matched_genes,
                'hits':hits,
                'RES':RES}  
        
    sig_res2d =   pd.DataFrame(sigres_df)

    # keep the top 25 values
    if len(sig_res2d)>25:
        sig_res2d = sig_res2d.sort_values(by="fdr").reset_index().loc[:26]

    return (sig_results,sig_res2d)




def NES_plot(df_result, outdir = False , dataset=""):

    if os.path.isfile(dataset): # If it is a filepath keep just the name of the file, without the extension
        dataset = dataset.split("/")[-1].split(".")[0]


    df_result.sort_values(by=["nes","fdr"],inplace=True) # Sort the values for visualization

    # NES PLOTS
    high = 0.3*len(df_result) # Adjust the plot high by the number of terms

    plt.figure(figsize=(5,high))

    # ESTO PERMITE QUE LUEGO SE COLOREEN LAS BARRAS SEGÚN EL -LOG10 (FDR) LOS MÁS SIGNIFICACIVOS SON MÁS OSCUROS
    data = np.array(-np.log10(df_result.fdr+0.0001))
    pal = sns.color_palette("hot", len(data))
    rank = data.argsort().argsort() 

    # EL PLOT EN SÍ
    ax = sns.barplot(
        data=df_result, x="nes", y="Term", 
        linewidth=4, palette= np.array(pal[::-1])[rank])

    # ESTO PERMITE PLOTEAR UNA LEYENDA PARA EL RANGO DE COLORES
    norm = plt.Normalize(data.min(),data.max())
    sm= plt.cm.ScalarMappable(cmap="hot_r", norm = norm)
    sm.set_array([])

    if len(df_result)<15:
        cbar = plt.colorbar(sm, location= "right",shrink = 0.4,anchor=(0,1),ax=ax)
    else:
        cbar = plt.colorbar(sm, location= "right",shrink = 0.2,anchor=(0,1),ax=ax)
    cbar.set_label("- Log10 ( FDR )",fontsize=8)

    # TÍTULO Y DEMAS PARÁMETROS DEL BARPLOT
    title = f"NES PLOT OF THE {dataset.upper()} DATASET\n"
    plt.title(title, fontweight=700,fontsize=15)
    plt.ylabel("")
    plt.xlabel("NES", fontsize=15)

    if not outdir:
        plt.show()
    else:
        plt.savefig(f"{outdir}/{dataset}_NES_plot.png",bbox_inches='tight')
    plt.close()
    
    
def run_ora(genes,datasets,save_path=False,figsize=(3,7),organism='human'):
    
    """
    Parameters:
    
    - genes: list of genes to be tested in the ora.
    
    - datasets: list of datasets # look at geseapy.get_library_name()
    
    - save_path: path where de dotplot will be saved.
    
    
    ==================================================================
    
    Returns: None. If savepath is provided it generates a dotplot of the significative terms in these path.
    
    """
    enr = gp.enrichr(gene_list=genes,
                    gene_sets=datasets,
                    organism=organism, 
                    outdir=None, # don't write to disk
                    )


    try:
        # categorical scatterplot
        ax = dotplot(enr.results,
                    column="Adjusted P-value",
                    x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
                    size=6,
                    top_term=10,
                    figsize=figsize,
                    title = "ORA",
                    xticklabels_rot=45, # rotate xtick labels
                    show_ring=True, # set to False to revmove outer ring
                    marker='o',
                    )
        
        if save_path:
            plt.savefig(save_path,bbox_inches="tight")
        plt.close()
        
    except ValueError:
        print('No enrich terms when cutoff = 0.05')