import os
import csv
from csv import writer
import os
import urllib.request
import Bio.PDB
from Bio.PDB import PDBParser

"""
New_File_Name_list=[]
for i in range(84):
    New_File_Name_list.append(f'rcsb_pdb_custom_report_{i}.csv')
counter=0
entries = os.listdir('/Users/c00479477/Downloads/Downloaded_PDB_files')
for entry in entries:
    os.rename(entry, New_File_Name_list[counter])
    counter+=1

"""


def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        return None

outputfile=open('PTR_SEP_TRO_aa_PDB_list.txt','w')
Path='/Users/c00479477/Desktop/PycharmProjects/pythonProject/Downloaded_PDB_PDBdatabase'
Downloaded_PDB_path='/Users/c00479477/Desktop/PycharmProjects/pythonProject/PDB_Files'
entries = os.listdir(Path)
for entry in entries:
    file = open(f'{Path}/{entry}', 'r')
    reader = csv.reader(file)
    next(reader)
    next(reader)
    for row in reader:
        download_pdb(row[1], Downloaded_PDB_path)
        Path_ = f'{Downloaded_PDB_path}/{row[1]}.pdb'
        if os.path.exists(Path_):
            p = Bio.PDB.PDBParser()
            Structure = p.get_structure('PrimaryStructureChain', Path_)
            model = Structure[0]
            for chain in model:
                for residue in chain:
                    resName = str(residue)[9:12].strip()
                    if resName == 'PTR' or resName == 'SEP' or resName == 'TPO':
                        numeric_filter = filter(str.isdigit, str(residue.id))
                        Res_Id = "".join(numeric_filter)
                        outputfile.write(f'{row[1]}\t{chain.id}\t{resName}\t{Res_Id}')
        else:
            print(f'{row[1]} could not be downloaded')

outputfile.close()
print('completed Successfully')



