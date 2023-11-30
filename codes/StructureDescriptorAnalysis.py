# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:04:21 2022

@author: PI
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import xlrd
import xlwt


SDA = ['a', 'b', 'c', 'd']

for i in range(len(SDA)):
    # Open Excel file using xlrd
    excel_path = r'..\data\{}.xlsx'.format(SDA[i])
    excel_file = xlrd.open_workbook(excel_path)

    # Get the first sheet
    sheet = excel_file.sheets()[0]


    # Get the number of rows
    rows = sheet.nrows

    # Get the second column as a list
    col = sheet.col_values(1)
    col.pop(0)

    # Create Descriptors.xlsx to store descriptors and their meanings
    Descriptors_xlsx = xlwt.Workbook(encoding='utf-8')
    DX_sheet1 = Descriptors_xlsx.add_sheet("{}_validation_descriptors".format(SDA[i]))

    # Create Des_data.xlsx to store computed descriptor data
    Des_data_xlsx = xlwt.Workbook(encoding='utf-8')
    DDX_sheet1 = Des_data_xlsx.add_sheet('{}_validation_descriptors_data'.format(SDA[i]))

    # Import the SMILES structure of the first molecule
    mol = Chem.MolFromSmiles(col[0])

    # Get all calculable descriptors and generate a list
    d_list = [x[0] for x in Descriptors._descList]

    # Calculate descriptors for the first molecule
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
    descList_note = calculator.GetDescriptorSummaries()
    x_list = [x for x in descList_note if x != 'N/A']

    # Generate excel with descriptors and their meanings
    for i in range(0, len(d_list)):
        DX_sheet1.write(i, 0, d_list[i])
        DX_sheet1.write(i, 1, descList_note[i])

    # Calculate descriptors for the first molecule and store data
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
    calculator = list(calculator.CalcDescriptors(mol))

    for i in range(0, len(d_list)):
        DDX_sheet1.write(0, i + 1, d_list[i])
        DDX_sheet1.write(1, i + 1, calculator[i])
    DDX_sheet1.write(0, 0, "SMILES")
    for i in range(len(col)):
        DDX_sheet1.write(i + 1, 0, col[i])

    # Iterate over the remaining molecules
    for i in range(1, len(col)):
        mol = Chem.MolFromSmiles(col[i])
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(d_list)
        calculator = list(calculator.CalcDescriptors(mol))

        for j in range(0, len(d_list)):
            DDX_sheet1.write(i + 1, j + 1, calculator[j])

    # Save the Des_data.xlsx workbook
    Des_data_xlsx.save(r'..\data\{}_validation_Molecular_Descriptor.xlsx'.format(SDA[i]))


