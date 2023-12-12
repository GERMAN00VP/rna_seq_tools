"""
Module created for handling of gmt files by: Germán Vallejo Palma
github:  https://github.com/GERMAN00VP

"""

import glob
from datetime import datetime
import os
import shutil


def gmt_to_gene_list(signature_path):

    """
    Input: .gmt file path (str)
    Output: gene list (list) | dictionary of gene lists (dict) if the gmt file contains multiple terms.
    """

    if signature_path[-4:]!= ".gmt":
        print("Make sure you have a gmt file, this file doesnt have a .gmt extenssion")

    signatures = open(signature_path,"r").read().split("\n")
    if len(signatures)==1:
        return open(signature_path,"r").read().split("\t")[2:]
    else:
        print("gmt file with multiple signatures, returning a dictionary!")
        return {signature.split("\t")[0] : signature.split("\t")[2:]  for signature in signatures}
    

def agraggate_signatures(path_dir_signgatures,new_signature_path=None):

    """
    Input: 
    
    - path_dir_signgatures: path to a directory with gmt signatures.

    - new_signature_path: path for the output file.

    Output: Creates a gmt file with a row for each gmt file in the input dir.
    """

    if path_dir_signgatures[-1] != "\\" or path_dir_signgatures[-1] != "/":
        path_dir_signgatures= path_dir_signgatures+"/"

    if new_signature_path==None:
        new_signature_path= path_dir_signgatures[:-2] +".gmt"

    new_signature= "\n".join([open(signature,"r").read() for signature in glob.glob(f"{path_dir_signgatures}*.gmt")])

    if len(new_signature)!=0:

        with open(new_signature_path,"w") as f:
            f.write(new_signature)
        
        print(f"The signature has been succesfully agregated i the file: {new_signature_path}")

    else:
        print("The provided path doesnt contain gmt files!")


def gene_list_to_gmt(gene_list,name="Customized_gene_set",description="Customized gmt",filepath=None, vervose=True):

        """
        Input: 
                - gene_list: gene info (list or dict).

                - name: if gene_list is a list it will be the name of the unique term in the file. Default: Customized_gene_set.

                - description: the second field in the gmt file (str). Default: Customized gmt.
                
                - filepath: path to save the gmat file, if filepath == None (the default) the path will be: f"./Customized_gene_set_{timestamp}.gmt".

        Output: 
                - Creates a .gmt file with the gene_list genes.

        Description: 
        This function takes a list of genes and creates a gmt file, if a dictionary is provided (dict="Term_name":gene_list) it takes the keys as term names and 
        values as gene lists, then creates a gmt with multiple rows.

        """
        timestamp= "_".join(str((datetime.timestamp(datetime.now()))).split("."))

        if filepath ==None:
                # Assign the default filepah
                filepath = f"./Customized_gene_set_{timestamp}.gmt"

        if type(gene_list)==dict:
                
                directory = f"./temp_{timestamp}/"
                os.mkdir(directory)

                for i,gene_set in enumerate(gene_list.keys()):
                        gene_list_to_gmt(gene_list[gene_set],name=gene_set,filepath=f"{directory}gene_set_{i}.gmt",vervose=False) 
                
                agraggate_signatures(path_dir_signgatures=directory,new_signature_path=filepath)
                
                shutil.rmtree(dirPath)


        elif type(gene_list)==list:

                if vervose:
                        print(f"Writing gmt file to {filepath}")
                text = "\t".join([name,description]+gene_list)
                with open(filepath,"w") as f:
                        f.write(text)

        else:
                print(f"The type of gene_list argument must be list or dict, you passed {type(gene_list)}")   