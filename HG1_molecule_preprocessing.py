#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
#if os.name == 'nt' or platform == 'win32':
#    print("Running on Windows")
#    if 'ptvsd' in sys.modules:
#        print("Running in Visual Studio")
#        try:
#            os.chdir(os.path.dirname(__file__))
#            print('CurrentDir: ', os.getcwd())
#        except:
#            pass
##--------------------------------------------------#
#    else:
#        print("Running outside Visual Studio")
#        try:
#            if not 'workbookDir' in globals():
#                workbookDir = os.getcwd()
#                print('workbookDir: ' + workbookDir)
#                os.chdir(workbookDir)
#        except:
#            pass
#--------------------------------------------------#
if os.name != 'nt' and platform != 'win32':
    print("Not Running on Windows")
#--------------------------------------------------#
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Fingerprints import FingerprintMols

#--------------------------------------------------#
import sys
import time
import numpy
import pickle
import typing
import itertools
from tqdm import tqdm
from copy import deepcopy
from pprint import pprint
from typing import Optional, Union, Tuple, Type, Set
#--------------------------------------------------#
import numpy as np
import pandas as pd
#--------------------------------------------------#
from PIL import Image
from cairosvg import svg2png

#--------------------------------------------------#
from HG_figure import *
from HG_rdkit import *
#--------------------------------------------------#





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def get_all_smiles_from_file(smiles_file):
    smiles_list = []
    with open(smiles_file) as f:
        lines = f.readlines()
        for one_line in lines:
            smiles_list.append(one_line.replace("\n", ""))
    return smiles_list

#============================================================================================================================#
def get_rxn_portion_from_smiles_list(smiles_list        : list  , 
                                     rxn_portion_SMARTS : tuple , ) \
                                     -> Tuple[dict, list]: # { smiles_str : mapping_list}, [mapping_list_i]
    #--------------------------------------------------#
    # Pattern Matching
    smiles_unannotated = []
    substrate_mapping_dict = dict([])
    for one_smiles in smiles_list:
        substrate_mapping_list = patterns_list_retrieving_AP(one_smiles, substructure_smarts_list = rxn_portion_SMARTS)
    #--------------------------------------------------#
    # Hard Coding for identifying substructure.
        if len(substrate_mapping_list[0])!=0:
            substrate_mapping_dict[one_smiles] = substrate_mapping_list[0]
        elif len(substrate_mapping_list[1])!=0:
            substrate_mapping_dict[one_smiles] = substrate_mapping_list[1]
        elif len(substrate_mapping_list[2])!=0:
            substrate_mapping_dict[one_smiles] = substrate_mapping_list[2]
        else:
            smiles_unannotated.append(one_smiles)
            print("unannotated smiles found!")
    #[print(i, substrate_mapping_dict[i]) for i in substrate_mapping_dict ]
    #--------------------------------------------------#
    # plot unannotated smiles
    print("number of smiles failed to identify reacting portion: ", len(smiles_unannotated))
    '''
    plot_smiles_list(   smiles_list = smiles_unannotated,
                        fig_folder = output_folder, 
                        img_size = (500, 500), 
                        molsPerRow = 5, 
                        separate_img = False, 
                        fig_name = dataset_nme)
                        '''
    #--------------------------------------------------#
    # Get rxn_portion
    substrate_rxn_portion_dict = dict([])
    substrate_rxn_portion_list = []
    for one_smiles in substrate_mapping_dict:

        mapping_set = set()
        for mapping in (substrate_mapping_dict[one_smiles]):
            mapping_set = mapping_set.union(set(mapping))
        substrate_rxn_portion_dict[one_smiles] = list(mapping_set)
        substrate_rxn_portion_list.append(list(mapping_set))
    return substrate_rxn_portion_dict, substrate_rxn_portion_list




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def main_x(dataset_nme, rxn_portion_SMARTS_dict, input_folder, output_folder, output_nme ,output_temp_folder):
    #============================================================================================================================#
    # Get all SMILES.
    smiles_list = get_all_smiles_from_file(input_folder / substrates_file)
    plot_smiles_list(   smiles_list = smiles_list,
                        fig_folder = output_folder, 
                        img_size = (500, 500), 
                        molsPerRow = 10, 
                        separate_img = False, 
                        fig_name = dataset_nme)
    #--------------------------------------------------#
    # Get substrate_rxn_portion_dict
    rxn_portion_SMARTS = rxn_portion_SMARTS_dict[dataset_nme]
    substrate_rxn_portion_dict, _ = get_rxn_portion_from_smiles_list(smiles_list, rxn_portion_SMARTS)
    #[print(i, substrate_rxn_portion_dict[i]) for i in substrate_rxn_portion_dict ]

    

    #============================================================================================================================#
    # Obtain graph components for each smiles.
    # components list:
    # 1. Atom_Attributes             size: (n_node, n_attr)               type: numpy.array
    # 2. Atom_RXN_Portion            size: (n_rxn_portion)                type: List
    # 3. Bond_Adjacency_Matrix       size: (n_node, n_node)               type: numpy.array
    # 4. Bond_Attributes             size: (n_node, n_node, dim_attr)     type: numpy.array
    # 5. Node_info
    Atom_Attributes_list = []
    Atom_RXN_Portion_list = []
    Bond_Adjacency_Matrix_list = []
    Bond_Attributes_list = []

    #print("locals(): ", locals())

    for one_smiles in smiles_list:
        #--------------------------------------------------#
        # 1. Get Atom_Attributes
        Atom_Attributes_list.append(smiles_to_nodes_encodings(one_smiles, smiles_list, radius = 3))
        Atom_RXN_Portion_list.append(substrate_rxn_portion_dict[one_smiles])

        Bond_Adjacency_Matrix, Bond_Attributes = smiles_to_bond_matrices(one_smiles)
        Bond_Adjacency_Matrix_list.append(Bond_Adjacency_Matrix)
        Bond_Attributes_list.append(Bond_Attributes)



    Graph_Attributes = ["Atom_Attributes_list", 
                        "Atom_RXN_Portion_list", 
                        "Bond_Adjacency_Matrix_list", 
                        "Bond_Attributes_list", ]


    Substrates_Graph_Attributes = { "Atom_Attributes_list"       : Atom_Attributes_list,
                                    "Atom_RXN_Portion_list"      : Atom_RXN_Portion_list,
                                    "Bond_Adjacency_Matrix_list" : Bond_Adjacency_Matrix_list,
                                    "Bond_Attributes_list"       : Bond_Attributes_list,
                                    }

    pickle.dump(Substrates_Graph_Attributes, open(output_folder / output_nme,"wb"))

    return
#######################################################################################################################################
#######################################################################################################################################








#######################################################################################################################################
#######################################################################################################################################
if __name__ == '__main__':
    # Directory and Files
    input_folder = Path("HG_data/")
    output_folder = Path("HG_results/")
    output_temp_folder = Path("HG_results/temp/")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_temp_folder):
        os.makedirs(output_temp_folder)

    #--------------------------------------------------#
    dataset_nme = ["phosphatase", "kinase", "halogenase", "esterase"][0]
    substrates_file = "X00_" + dataset_nme + ".smiles"
    smiles_list = get_all_smiles_from_file(input_folder / substrates_file)

    output_nme = "X01_" + dataset_nme + "_Substrates_Graph_Attributes.p"
    #--------------------------------------------------#
    # Interested Substructures
    rxn_portion_SMARTS_dict = { "phosphatase": ("[P](=[O])([OH])([OH])", "[O]([P](=[O])([OH]))([P](=[O])([OH]))", "[P](=[O])([OH])"),
                                "kinase"     : "[#15]",
                                "halogenase" : "[#15]",
                                "esterase"   : "[#15]",
                                                            }

    main_x(dataset_nme, rxn_portion_SMARTS_dict, input_folder, output_folder, output_nme ,output_temp_folder)

 