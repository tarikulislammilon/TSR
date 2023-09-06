#Program to search specfic keys consist of thrianges CA, CB and CG atoms.
#Author: Tarikul Islam Milon
#Created on: 06/30/2023

import os
import csv #Changed


# folder path
dir_path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/output_CA'

OutputFile=open('CA_CB_CG_key_extraction.txt','w')

File=open("sample_details_TYR.csv",'r')
reader=csv.reader(File)
next(reader)

atom_list=['CA','CB','CG']
for row in reader:
    FileName = row[0]
    f = open(f'{dir_path}/{FileName}.triplets_theta29_dist17', 'r')
    next(f)
    for row in f:
        data = row.split()
        Atom1 = data[0].split("_")[3]
        Atom2 = data[1].split("_")[3]
        Atom3 = data[2].split("_")[3]
        if Atom1 in atom_list and Atom2 in atom_list and Atom3 in atom_list:
            OutputFile.write(f'{FileName}\t{row}\n')

print('Completed Successfully')

