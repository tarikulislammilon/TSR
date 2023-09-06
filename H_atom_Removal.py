#Program to remove Hydrogen atom from PDB file
#Author:Tarikul Islam Milon

import pandas as pd
import os
import csv

list_file = open("sample_details_17.csv",'r')
reader=csv.reader(list_file)
next(reader)
for row in reader:
    pdb=row[0]
    if pdb=='1T3R':
        PDB_File = open('/Users/c00479477/Desktop/PycharmProjects/pythonProject/PDB_datadir_drug/{}.pdb'.format(pdb),   'r')
        PDB_File_Without_H_Atom = open(  '/Users/c00479477/Desktop/PycharmProjects/pythonProject/PDB_datadir_H_removed/{}.pdb'.format(pdb), 'w')
        lines = PDB_File.readlines()
        for line in lines:
            Atom_type = line[0:6].strip()
            if Atom_type != 'HETATM':
                PDB_File_Without_H_Atom.writelines(line)
            elif Atom_type == 'HETATM':
                data=line.split()
                if data[3]=='017':
                    Line=line.replace("017",'A17')
                    PDB_File_Without_H_Atom.writelines(Line)

        PDB_File.close()
        PDB_File_Without_H_Atom.close()



 











