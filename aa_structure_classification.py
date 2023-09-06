# python code to get the different structure group information of amino acids in a protein
# This program groups the amino acids based on two different structure 1) alpha-Helix 3)Beta sheet 2) Linear chain
# The x coordinate of each residue as added to each aa residue data file to take into account  100A, 100B, 100C are differnt residue
#Input files :1) Sample details 2) PDB File containing folder path
#Author:Tarikul Islam Milon
#Created on:03/04/2023

import csv
import math
import Bio.PDB
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np

AminoAcidName='TYR' # change here for different aa
df=pd.read_csv('sample_details_s_h_l_mix5_original_pdb.csv')
PDB_list=df['protein'].to_list()
Chain=df['chain'].to_list()

structure_info_dict={"helix":[],"sheet":[],"linear_chain":[]}
outputFile=open('aa_structure_classification.txt','w')


def get_structure_info(PDB,**kwargs): # Gives different structures info (i.e. Helix,sheet)
    f=open(f"{PDB}.pdb")
    lines=f.readlines()
    h={}
    s={}
    counter_h=0
    counter_s=0
    for line in lines:
        if line[0:6].strip()=='HELIX':
            heLX_chain=line[19]
            heLX_seqStart=line[20:25].strip()
            heLX_seqEnd=line[33:37].strip()
            for i in np.arange(int(heLX_seqStart),int(heLX_seqEnd)+1):
                h[f'{heLX_chain}_{i}']=counter_h
                counter_h+=1

        elif line[0:6].strip()=='SHEET':
            Sheet_chain=line[21]
            Sheet_seqStart=line[22:26].strip()
            Sheet_seqEnd=line[33:37].strip()
            for i in np.arange(int(Sheet_seqStart),int(Sheet_seqEnd)+1):
                s[f'{Sheet_chain}_{i}']=counter_s
                counter_s+=1

    return h,s

for i in range(len(PDB_list)):
    PDB_ID = PDB_list[i]
    # Get the information of different structures
    helix, sheet = get_structure_info(PDB_ID)
    Chain_Name = Chain[i]
    PDB_File_Path = "{}.pdb".format(PDB_ID.lower())
    p = Bio.PDB.PDBParser()
    Structure = p.get_structure('PrimaryStructureChain', PDB_File_Path)  # Get the protein structure
    model = Structure[0]
    for chain in model:
        if chain.id == Chain_Name:
            for residue in chain:
                if str(residue)[9:12] == AminoAcidName.upper():
                    numeric_filter = filter(str.isdigit, str(residue.id)[6:10])
                    Res_Id = "".join(numeric_filter)
                    for atom0 in residue:
                        Coord_end = atom0.get_vector()
                        X1_end = Coord_end[0]
                    X1_end = '{:.2f}'.format(X1_end)  # last x coordinate of the aa residue

                    structure_info = f'{Chain_Name}_{Res_Id}'
                    data = f'{PDB_ID}_{Chain_Name}_{AminoAcidName}_{Res_Id}_{X1_end}'
                    if structure_info in helix:
                        structure_info_dict["helix"].append(data)
                    elif structure_info in sheet:
                        structure_info_dict["sheet"].append(data)
                    else:
                        structure_info_dict["linear_chain"].append(data)

# Writting output to file
for structure in structure_info_dict:
    outputFile.write(f'{structure}\n')
    for aa in structure_info_dict[structure]:
        outputFile.write(f'{aa}\n')

outputFile.close()
print('Completed Successfully')



