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
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
from alfabet import model
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
from typing import Optional, Union, Tuple, Type, Set, List, Dict
#--------------------------------------------------#
import numpy as np
import pandas as pd
#--------------------------------------------------#
from PIL import Image
from cairosvg import svg2png

#--------------------------------------------------#
from HG_figure import *
from AP_convert import *
#--------------------------------------------------#





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def get_all_smiles_from_file(smiles_file):  #####
    smiles_list = []
    with open(smiles_file) as f:
        lines = f.readlines()
        for one_line in lines:
            smiles_list.append(one_line.replace("\n", ""))
    return smiles_list

#============================================================================================================================#
def patterns_list_retrieving_AP(cmpd_smiles, substructure_smarts_list, bad_ss_dict={}):
    #####----------return mappings of substructures
    mol_x = MolFromSmiles_AP(cmpd_smiles, bad_ss_dict)
    try:
        substructure_mapping_list = []
        for substructure_smarts in substructure_smarts_list:
            pattern_matching = pattern_matching_AP(cmpd_smiles, substructure_smarts, bad_ss_dict)
            if pattern_matching:
                substructure_mapping = mol_x.GetSubstructMatches(Chem.MolFromSmarts(substructure_smarts))
                substructure_mapping_list.append(substructure_mapping)
            else:
                substructure_mapping_list.append( tuple([]) )
    except Exception:
        substructure_mapping_list = [ tuple([]) for i in range(len(substructure_mapping_list)) ]
        print("patterns_list_retrieving_AP ERROR, ", cmpd_smiles, "\nsubstructure_mapping_list: ", substructure_mapping_list)
    return substructure_mapping_list

#============================================================================================================================#
def smiles_to_node_Morgan_dict(smiles_x, radius = 3):
    #--------------------------------------------------#
    # Inputs
    mol_x = Chem.MolFromSmiles(smiles_x)
    num_atom_x = len([atom for atom in mol_x.GetAtoms()])
    # Get Morgan dict
    Morgan_dict = {}
    Morgan_fingerprint = AllChem.GetMorganFingerprint(mol_x, radius = radius, bitInfo = Morgan_dict)
    #--------------------------------------------------#
    #[print(one_Morgan_info, Morgan_dict[one_Morgan_info]) for one_Morgan_info in Morgan_dict]
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    #--------------------------------------------------#
    Morgan_dict_reverse = {}
    for one_MorganFP in Morgan_dict:
        for one_substructure in Morgan_dict[one_MorganFP]:
            Morgan_dict_reverse[one_substructure] = one_MorganFP
    #print(Morgan_dict_reverse)
    #--------------------------------------------------#
    # Output
    node_Morgan_dict = {}
    for one_idx in atom_idx_list_x:
        node_Morgan_dict[one_idx] = []
    node_Morgan_list = []
    
    for one_idx in atom_idx_list_x:
        for one_radius in range(radius):
            one_substructure_mapping = (one_idx, one_radius)
            #print(one_substructure_mapping)
            if one_substructure_mapping not in Morgan_dict_reverse:
                continue
            node_Morgan_dict[one_idx].append(Morgan_dict_reverse[one_substructure_mapping])
    node_Morgan_list = [[one_node, node_Morgan_dict[one_node]] for one_node in node_Morgan_dict]
    #--------------------------------------------------#
    return node_Morgan_dict, node_Morgan_list

#============================================================================================================================#
def smiles_list_to_all_Morgan_list(smiles_list: List[str], radius: int = 3) -> List: 
    all_Morgan_set = set([])
    for one_smiles in smiles_list:
        _, node_Morgan_list = smiles_to_node_Morgan_dict(one_smiles, radius = radius)
        for one_node_info in node_Morgan_list:
            all_Morgan_set = all_Morgan_set.union((one_node_info[1]))
    all_Morgan_list = list(all_Morgan_set)
    return all_Morgan_list

#============================================================================================================================#
def smiles_to_all_nodes_Morgan_list_dict(smiles_x: str, radius: int = 3) -> Tuple[List, Dict]:
    #--------------------------------------------------#
    # Inputs
    mol_x = Chem.MolFromSmiles(smiles_x)
    num_atom_x = len([atom for atom in mol_x.GetAtoms()])
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    #--------------------------------------------------#   
    all_Morgan = set([])
    all_nodes_Morgan_dict = dict([])
    node_Morgan_dict, node_Morgan_list = smiles_to_node_Morgan_dict(smiles_x, radius = radius)
    count_x = 0

    for one_atom_id in atom_idx_list_x:
        #print(count_x," out of ", len(atom_idx_list_x) )
        count_x += 1
        all_Morgan = all_Morgan.union(set(node_Morgan_dict[one_atom_id]))
    all_nodes_Morgan_dict = node_Morgan_dict
    return list(all_Morgan), all_nodes_Morgan_dict

#============================================================================================================================#
def smiles_to_bond_bde_dict(smiles_x, drop_duplicates = False):
    #--------------------------------------------------#
    # Input SMILES
    mol_x = Chem.MolFromSmiles(smiles_x)
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    #--------------------------------------------------#
    # Get bond_index_bde_dict.
    # {idx : (bond_type, bond_bde) }
    df_alfabet_pred_BDE = model.predict([smiles_x], drop_duplicates = drop_duplicates) #BDE: bond desociation energy
    num_row = len(df_alfabet_pred_BDE)
    bond_index_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bond_index"] for one_row in range(num_row)}
    bond_type_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bond_type"] for one_row in range(num_row)}
    bond_bde_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bde_pred"] for one_row in range(num_row)}
    bond_index_bde_dict = {bond_index_dict[one_row] : (bond_type_dict[one_row], bond_bde_dict[one_row])  for one_row in range(num_row)}
    return bond_index_bde_dict

#============================================================================================================================#
def smiles_to_bond_info(smiles_x):
    #--------------------------------------------------#
    # Input SMILES.
    mol_x = Chem.MolFromSmiles(smiles_x)
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    #--------------------------------------------------#
    # Get atom and bond attributes.
    df_atom_attributes, df_bond_attributes, atom_dict, bond_dict = smiles_to_attributes(smiles_x)
    df_bond_info = df_bond_attributes[["bond_idx", "BondTypeAsDouble", "bond_idx_tuple"]]
    #--------------------------------------------------#
    # Get bond energy for bonds with order >= 2.
    bond_type_SMARTS_dict = {1.0 : "-", 2.0 : "=", 3.0 : "#", 1.5 : ":"}
    bde_kJ_dict = { "C=C": 602, "C:C": 518, "C#C": 835, 
                    "C=O": 749, "O=P": 544, "C:N": 508, 
                    "C=N": 615, "C#N": 887, "C=S": 573, 
                    "N=O": 607, "P=S": 335, "C:O": 621,
                    "C:S": 456,

                    
                    "C-C": 418, "C-O": 470, "C-S": 317, # on ring

                    "X": 420}
    bde_kcal_dict = {one_bond : bde_kJ_dict[one_bond] * 0.239006 for one_bond in bde_kJ_dict}

    bond_index_bde_dict = smiles_to_bond_bde_dict(smiles_x)
    bond_index_len_dict = smiles_to_bond_len_dict(smiles_x)

    #--------------------------------------------------#
    # Get bond string. (e.g., "C-O")
    bond_str_list = []
    bond_idx_str_dict = {}
    for bond_idx in bond_idx_list_x:
        bond_middle = bond_type_SMARTS_dict[bond_dict[bond_idx]["BondTypeAsDouble"]]
        bond_head = atom_dict[bond_dict[bond_idx]["bond_idx_tuple"][0]]["Symbol"]
        bond_tail = atom_dict[bond_dict[bond_idx]["bond_idx_tuple"][1]]["Symbol"]
        one_bond_str = bond_head + bond_middle + bond_tail
        one_bond_str_rev = one_bond_str[::-1]
        one_bond_str = min((one_bond_str, one_bond_str_rev))
        bond_str_list.append(one_bond_str)
        bond_idx_str_dict[bond_idx] = one_bond_str
    df_bond_info["bond_str"] = bond_str_list
    #--------------------------------------------------#
    # Get df_bond_info
    bond_bde_list = []
    bond_len_list = []
    for bond_idx in bond_idx_list_x:
        if bond_idx in bond_index_bde_dict.keys():
            bond_bde_list.append(bond_index_bde_dict[bond_idx][1])
        else:
            if bond_idx_str_dict[bond_idx] in bde_kcal_dict:
                bond_bde_list.append(bde_kcal_dict[bond_idx_str_dict[bond_idx]])
            else:
                print("\nCAUTION: Found Unknown Bond Type !!! " + bond_idx_str_dict[bond_idx])
                print(smiles_x, bond_index_bde_dict)
                bond_bde_list.append(bde_kcal_dict["X"])
        bond_len_list.append(bond_index_len_dict[bond_idx])
        
    df_bond_info["bond_bde"] = bond_bde_list
    df_bond_info["bond_len"] = bond_len_list

    return df_bond_info

#============================================================================================================================#
def smiles_to_attributes(smiles_x):
    #--------------------------------------------------#
    # Input SMILES
    mol_x = Chem.MolFromSmiles(smiles_x)
    #============================================================================================================================#
    # Go through all atoms, get all properties and write to a dict.
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    atom_dict = {}
    # atom_dict             <-  { atom_idx  : atom_attributes_dict }
    # atom_attributes_dict  <-  { atom_attr : atom_value }
    #--------------------------------------------------#
    # Atoms
    for atom in mol_x.GetAtoms():
        # Add mapping number to SMARTS, mapping number = idx + 1. (mapping number cant be 0).
        atom.SetAtomMapNum(atom.GetIdx()+1)
        # Get All Attributes of Each Atom.
        atom_attributes_dict = {}
        atom_attributes_dict["atom_idx"]            = atom.GetIdx()
        atom_attributes_dict["atom_map_num"]        = atom.GetAtomMapNum()
        atom_attributes_dict["atomic_num"]          = atom.GetAtomicNum()
        atom_attributes_dict["ExplicitValence"]     = atom.GetExplicitValence()
        atom_attributes_dict["FormalCharge"]        = atom.GetFormalCharge()
        atom_attributes_dict["ImplicitValence"]     = atom.GetImplicitValence()
        atom_attributes_dict["IsAromatic"]          = atom.GetIsAromatic()
        atom_attributes_dict["Isotope"]             = atom.GetIsotope()
        atom_attributes_dict["Mass"]                = atom.GetMass()
        atom_attributes_dict["MonomerInfo"]         = atom.GetMonomerInfo()
        atom_attributes_dict["NoImplicit"]          = atom.GetNoImplicit()
        atom_attributes_dict["NumExplicitHs"]       = atom.GetNumExplicitHs()
        atom_attributes_dict["NumImplicitHs"]       = atom.GetNumImplicitHs()
        atom_attributes_dict["NumRadicalElectrons"] = atom.GetNumRadicalElectrons()
        atom_attributes_dict["PDBResidueInfo"]      = atom.GetPDBResidueInfo()
        atom_attributes_dict["Smarts"]              = atom.GetSmarts()
        atom_attributes_dict["Symbol"]              = atom.GetSymbol()
        atom_attributes_dict["TotalDegree"]         = atom.GetTotalDegree() # Degree is defined to be its number of directly-bonded neighbors.
        atom_attributes_dict["TotalNumHs"]          = atom.GetTotalNumHs()
        atom_attributes_dict["TotalValence"]        = atom.GetTotalValence()
        atom_attributes_dict["ChiralTag"]           = atom.GetChiralTag()
        atom_attributes_dict["Neighbors"]           = [atom.GetIdx() for atom in atom.GetNeighbors()]
        atom_attributes_dict["Bonds"]               = [bond.GetIdx() for bond in atom.GetBonds()]
        atom_dict[atom.GetIdx()] = atom_attributes_dict
    #--------------------------------------------------#
    # Get dict and list to generate the dataframe.
    atom_attributes_key_list = list(atom_dict[0].keys())
    atom_attributes_value_dict = {attr : [atom_dict[one_atom][attr] for one_atom in atom_dict]  for attr in atom_attributes_key_list }
    #print(atom_attributes_value_dict)
    #print(atom_attributes_key_list)
    #--------------------------------------------------#
    # Get the Dataframe.
    df_atom_attributes = pd.DataFrame(data = atom_attributes_value_dict, columns = atom_attributes_key_list)
    #============================================================================================================================#
    # Go through all bonds, get all properties and write to a dict.
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    bond_dict = {}
    # bond_dict             <-  { bond_idx  : bond_attributes_dict }
    # bond_attributes_dict  <-  { bond_attr : bond_value }
    #--------------------------------------------------#
    # Bonds
    for bond in mol_x.GetBonds():
        bond_attributes_dict = {}
        # Get All Attributes of Each Bond.
        bond_attributes_dict["bond_idx"]           =     bond.GetIdx()
        bond_attributes_dict["BondDir"]            =     bond.GetBondDir()
        bond_attributes_dict["BondType"]           =     bond.GetBondType()
        bond_attributes_dict["BondTypeAsDouble"]   =     bond.GetBondTypeAsDouble()
        bond_attributes_dict["IsAromatic"]         =     bond.GetIsAromatic()
        bond_attributes_dict["IsConjugated"]       =     bond.GetIsConjugated()
        bond_attributes_dict["Smarts"]             =     bond.GetSmarts()
        bond_attributes_dict["Stereo"]             =     bond.GetStereo()
        bond_attributes_dict["bond_idx_tuple"]     =     (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        bond_dict[bond.GetIdx()] = bond_attributes_dict
    #--------------------------------------------------#
    # Get dict and list to generate the dataframe.
    bond_attributes_key_list = list(bond_dict[0].keys())
    bond_attributes_value_dict = {attr : [bond_dict[one_bond][attr] for one_bond in bond_dict]  for attr in bond_attributes_key_list }
    #print(bond_attributes_value_dict)
    #print(bond_attributes_key_list)
    #--------------------------------------------------#
    # Get the Dataframe.
    df_bond_attributes = pd.DataFrame(data = bond_attributes_value_dict, columns = bond_attributes_key_list)

    return df_atom_attributes, df_bond_attributes, atom_dict, bond_dict

#============================================================================================================================#
def smiles_to_bond_bde_dict(smiles_x, drop_duplicates = False):
    #--------------------------------------------------#
    # Input SMILES
    mol_x = Chem.MolFromSmiles(smiles_x)
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    #--------------------------------------------------#
    # Get bond_index_bde_dict.
    # {idx : (bond_type, bond_bde) }
    df_alfabet_pred_BDE = model.predict([smiles_x], drop_duplicates = drop_duplicates) #BDE: bond desociation energy
    num_row = len(df_alfabet_pred_BDE)
    bond_index_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bond_index"] for one_row in range(num_row)}
    bond_type_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bond_type"] for one_row in range(num_row)}
    bond_bde_dict = {one_row: df_alfabet_pred_BDE.loc[one_row, "bde_pred"] for one_row in range(num_row)}
    bond_index_bde_dict = {bond_index_dict[one_row] : (bond_type_dict[one_row], bond_bde_dict[one_row])  for one_row in range(num_row)}
    return bond_index_bde_dict

#============================================================================================================================#
def smiles_to_bond_len_dict(smiles_x):
    #--------------------------------------------------#
    # Input SMILES
    mol_x = Chem.MolFromSmiles(smiles_x)
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    #--------------------------------------------------#
    # Get bond_index_len_dict.

    AllChem.EmbedMolecule(mol_x, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol_x)

    bond_index_len_dict = {}
    for bond in mol_x.GetBonds():
        bond_length = rdMolTransforms.GetBondLength(mol_x.GetConformer(), bond.GetBeginAtomIdx(),bond.GetEndAtomIdx() )
        bond_index_len_dict[bond.GetIdx()] = bond_length
    #--------------------------------------------------#
    #print(bond_index_len_dict)
    return bond_index_len_dict




#######################################################################################################################################
#######################################################################################################################################
def get_rxn_portion_from_smiles_list(smiles_list        : list  , 
                                     rxn_portion_SMARTS : tuple , ) \
                                     -> Tuple[dict, list]: # { smiles_str : mapping_list}, [mapping_list_i] #####
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

#============================================================================================================================#
def smiles_to_nodes_encodings(smiles_x: str, smiles_list, radius: int = 3):
    #--------------------------------------------------#
    # Inputs
    mol_x = Chem.MolFromSmiles(smiles_x)
    num_atom_x = len([atom for atom in mol_x.GetAtoms()])
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    #--------------------------------------------------#
    # 
    all_Morgan_list = smiles_list_to_all_Morgan_list(smiles_list, radius = radius)
    _, all_nodes_Morgan_dict = smiles_to_all_nodes_Morgan_list_dict(smiles_x, radius = radius)
    node_attrs_dim = len(all_Morgan_list)
    nodes_encodings = []
    for one_node in atom_idx_list_x:
        Xi = [0] * node_attrs_dim
        Xi_Morgan_list = all_nodes_Morgan_dict[one_node]
        for one_Morgan in Xi_Morgan_list:
            Xi[all_Morgan_list.index(one_Morgan)] = Xi_Morgan_list.count(one_Morgan)
        nodes_encodings.append(Xi) # np.zeros((len(atom_idx_list_x), node_attrs_dim))

    nodes_encodings_np = np.array(nodes_encodings)
    #print(nodes_encodings_np.shape)
    return nodes_encodings_np

#============================================================================================================================#
def smiles_to_bond_matrices(smiles_x: str) -> numpy.ndarray: # size: (n_node, n_node)
    #--------------------------------------------------#
    # Input SMILES
    mol_x = Chem.MolFromSmiles(smiles_x)
    bond_idx_list_x = [bond.GetIdx() for bond in mol_x.GetBonds()]
    atom_idx_list_x = [atom.GetIdx() for atom in mol_x.GetAtoms()]
    #--------------------------------------------------#
    #
    bond_index_bde_dict = smiles_to_bond_bde_dict(smiles_x)
    df_bond_info = smiles_to_bond_info(smiles_x)
    #--------------------------------------------------#
    # bond_adjacency_matrix
    bond_adjacency_matrix = np.zeros((len(atom_idx_list_x), len(atom_idx_list_x)))
    for one_bond_tuple in df_bond_info.loc[:, "bond_idx_tuple"]:
        bond_adjacency_matrix[one_bond_tuple[0]][one_bond_tuple[1]] = 1
        bond_adjacency_matrix[one_bond_tuple[1]][one_bond_tuple[0]] = 1
    #--------------------------------------------------#
    # bond_graph_attributes
    bond_graph_attributes = np.zeros((len(atom_idx_list_x), len(atom_idx_list_x), 2))
    for one_bond_tuple in df_bond_info.loc[:, "bond_idx_tuple"]:
        bond_bde = df_bond_info.loc[ df_bond_info["bond_idx_tuple"] == one_bond_tuple, ["bond_bde"]].values
        bond_len = df_bond_info.loc[ df_bond_info["bond_idx_tuple"] == one_bond_tuple, ["bond_len"]].values
        bond_graph_attributes[one_bond_tuple[0]][one_bond_tuple[1]][0] = bond_bde
        bond_graph_attributes[one_bond_tuple[1]][one_bond_tuple[0]][0] = bond_bde
        bond_graph_attributes[one_bond_tuple[0]][one_bond_tuple[1]][1] = bond_len
        bond_graph_attributes[one_bond_tuple[1]][one_bond_tuple[0]][1] = bond_len
    #print (bond_adjacency_matrix)
    #print (bond_graph_attributes)
    return bond_adjacency_matrix, bond_graph_attributes



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def main_x(dataset_nme, rxn_portion_SMARTS_dict, input_folder, output_folder, output_nme ,output_temp_folder):
    #============================================================================================================================#
    # Get all SMILES.
    smiles_list = get_all_smiles_from_file(input_folder / substrates_file)
    # plot_smiles_list(   smiles_list = smiles_list,
    #                     fig_folder = output_folder, 
    #                     img_size = (500, 500), 
    #                     molsPerRow = 10, 
    #                     separate_img = False, 
    #                     fig_name = dataset_nme)
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

 