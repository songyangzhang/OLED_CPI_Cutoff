# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:07:44 2022

@author: PI
"""


from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import xlrd
import xlwt


# Load molecular structures from an Excel file in mol file format
excel_file = xlrd.open_workbook(r'..\data\CutoffData.xlsx')

# Retrieve the first sheet from the Excel workbook
table = excel_file.sheets()[0]

# Get the number of rows in the sheet
rows = table.nrows

# Retrieve the second column's content, excluding the header
col = table.col_values(1)
col.pop(0)

# Create Descriptors.xlsx to store descriptors and their meanings
descriptors_xlsx = xlwt.Workbook(encoding='utf-8')
desc_sheet = descriptors_xlsx.add_sheet("descriptors")

# Create Des_data.xlsx to store computed descriptor data
des_data_xlsx = xlwt.Workbook(encoding='utf-8')
data_sheet = des_data_xlsx.add_sheet('descriptors_data')

# Import the SMILES structure of the first molecule
mol = Chem.MolFromSmiles(col[0])

# Generate a list of all calculable descriptors
d_list = [x[0] for x in Descriptors._descList]

# Calculate all descriptors for the molecule
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
desc_list_note = calculator.GetDescriptorSummaries()

# Remove descriptors without meanings
x_list = [x for x in desc_list_note if x != 'N/A']

# Generate an excel sheet with descriptors and their meanings
for i in range(0, len(d_list)):
    desc_sheet.write(i, 0, d_list[i])
    desc_sheet.write(i, 1, desc_list_note[i])

# Calculate descriptors for the first molecule and store data in Des_data.xlsx
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
calculator = list(calculator.CalcDescriptors(mol))

for i in range(0, len(d_list)):
    data_sheet.write(0, i + 1, d_list[i])
    data_sheet.write(1, i + 1, calculator[i])

data_sheet.write(0, 0, "SMILES")

for i in range(len(col)):
    data_sheet.write(i + 1, 0, col[i])

# Loop through each molecule's SMILES structure
for i in range(1, len(col)):
    mol = Chem.MolFromSmiles(col[i])

    # Calculate descriptors for the molecule
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
    calculator = list(calculator.CalcDescriptors(mol))

    # Write descriptor data to Des_data.xlsx
    for j in range(0, len(d_list)):
        data_sheet.write(i + 1, j + 1, calculator[j])

# Save the created Excel files
descriptors_xlsx.save('..\data\Molecular_Descriptor_Meaning.xlsx')
des_data_xlsx.save('..\data\Molecular_Descriptor.xlsx')


