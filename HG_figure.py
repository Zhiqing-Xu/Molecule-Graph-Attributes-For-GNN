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
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
import sys
import time
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
from AP_funcs import *
from AP_convert import *

#######################################################################################################################################
#######################################################################################################################################
def get_concat_h_blank(im1, im2, color=(1, 1, 1)):
    dst = Image.new('RGB', (im1.width + im2.width, max(im1.height, im2.height)), color)
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst
#============================================================================================================================#
def get_concat_v_blank(im1, im2, color=(1, 1, 1)):
    dst = Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height), color)
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst
#============================================================================================================================#
def get_concat_multi_h_blank(image_list, color=(1, 1, 1)):
    current_image = image_list[0]
    for i in range(len(image_list)-1):
        current_image = get_concat_h_blank(current_image, image_list[i+1], color=(1, 1, 1))
    return current_image
#============================================================================================================================#
def get_concat_multi_v_blank(image_list, color=(1, 1, 1)):
    current_image = image_list[0]
    for i in range(len(image_list)-1):
        current_image = get_concat_v_blank(current_image, image_list[i+1], color=(1, 1, 1))
    return current_image
#============================================================================================================================#
def get_concat_h_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    if im1.height == im2.height:
        _im1 = im1
        _im2 = im2
    elif (((im1.height > im2.height) and resize_big_image) or
          ((im1.height < im2.height) and not resize_big_image)):
        _im1 = im1.resize((int(im1.width * im2.height / im1.height), im2.height), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((int(im2.width * im1.height / im2.height), im1.height), resample=resample)
    dst = Image.new('RGB', (_im1.width + _im2.width, _im1.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (_im1.width, 0))
    return dst
#============================================================================================================================#
def get_concat_v_resize(im1, im2, resample=Image.BICUBIC, resize_big_image=True):
    if im1.width == im2.width:
        _im1 = im1
        _im2 = im2
    elif (((im1.width > im2.width) and resize_big_image) or
          ((im1.width < im2.width) and not resize_big_image)):
        _im1 = im1.resize((im2.width, int(im1.height * im2.width / im1.width)), resample=resample)
        _im2 = im2
    else:
        _im1 = im1
        _im2 = im2.resize((im1.width, int(im2.height * im1.width / im2.width)), resample=resample)
    dst = Image.new('RGB', (_im1.width, _im1.height + _im2.height))
    dst.paste(_im1, (0, 0))
    dst.paste(_im2, (0, _im1.height))
    return dst


#######################################################################################################################################
#######################################################################################################################################
def test():

    #--------------------------------------------------#
    # directory
    input_folder = Path("HG_test_in/")
    output_folder = Path("HG_test_out/")
    output_temp_folder = Path("HG_test_out/temp/")


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_temp_folder):
        os.makedirs(output_temp_folder)

    #--------------------------------------------------#
    im1 = Image.open(output_temp_folder / 'MorganFP_testlevel_0.png')
    im2 = Image.open(output_temp_folder / 'MorganFP_testlevel_1.png')
    im3 = Image.open(output_temp_folder / 'MorganFP_testlevel_2.png')

    get_concat_multi_h_blank([im1, im2, im3]).save('./test1.png')
    get_concat_multi_v_blank([im1, im2, im3]).save('./test2.png')

    get_concat_h_resize(im1, im2).save('./test3.png')
    get_concat_v_resize(im1, im2, resize_big_image=False).save('./test4.png')

#######################################################################################################################################
#######################################################################################################################################
if __name__ == '__main__':
    test()

 