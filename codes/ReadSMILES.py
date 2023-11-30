# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:04:22 2022

@author: PI
"""

# Import required libraries
import rdkit
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit import Chem
from rdkit.Chem import rdMolHash
import pandas as pd

# Function to draw a molecule from its SMILES representation and save the image
def draw_smiles(smiles, name):
    opts = DrawingOptions()
    m = Chem.MolFromSmiles(smiles)
    opts.includeAtomNumbers = True
    opts.bondLineWidth = 2.8
    draw = Draw.MolToImage(m, options=opts)
    draw.save(f'..\data\CutoffData.xlsx\structure\{name}.jpg')

# Function to generate Murcko scaffold hash for a list of molecules
def gujia(mList):
    mMols = [Chem.MolFromSmiles(m) for m in mList]
    murckoHashList = [rdMolHash.MolHash(mMol, rdkit.Chem.rdMolHash.HashFunction.MurckoScaffold) for mMol in mMols]
    return murckoHashList

# Function to find the most frequent element in a list
def mostFreq(lst):
    return max(set(lst), key=lst.count)

# Function to read SMILES representations from an Excel file
def read_smiles(path):
    all_smis = pd.read_excel(path)
    all_smis = list(all_smis['SMILES'])
    return all_smis

# Main block
if __name__ == '__main__':
    # Define the path to the Excel file containing SMILES representations
    excel_path = r"..\data\CutoffData.xlsx"
    
    # Read SMILES representations from the Excel file
    mList = read_smiles(excel_path)
    
    # Generate Murcko scaffold hash for each molecule
    m = gujia(mList)
    
    # Find the most frequent Murcko scaffold hash
    n = mostFreq(m)
    mostFreq_murckoHash_mol = Chem.MolFromSmiles(n)

    # Loop through the first 10 molecules, for example
    for i in range(10):
        # Create a list of molecules from SMILES representations
        mMols = [Chem.MolFromSmiles(smiles) for smiles in mList[i:i + 1]]
        
        # Find substructure matches for the most frequent Murcko scaffold in each molecule
