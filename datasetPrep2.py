# Python scripts for drug datasets preparation
# Author: Tarikul Islam Milon

import csv
from csv import writer
import urllib.request
import pandas as pd
import os

class DatasetPrep:
    def __init__(self,*args):
        self.mainFile=args[0]
        self.File_to_read=args[1]
        self.Exist_PDB_list=args[2]
        self.drug_name=args[3]

    def readInputCsvFile(self):
        PDB_list=[]
        reader=csv.reader(File_to_read)
        for row in reader:
            if len(row[0])!=0:
                current_PDB = row[0]

            if row[5]=='PNN' or row[5].strip()=='PNM' :
                if row[0] not in self.Exist_PDB_list:
                    List=[current_PDB,row[7],row[5]]
                    with open(f'{self.mainFile}', 'a') as f_object:
                        writer_object = writer(f_object)
                        writer_object.writerow(List)
                        f_object.close()

                    PDB_list.append(row[0])
        return PDB_list

# CSV file in which the data will be prepared
Drug_Name='PNN'
mainFile=f'rcsb_pdb_{Drug_Name}_data.csv'
#Make list using already existing PDB
Exist_PDB_list=[]
f=open(mainFile,'r')
reader=csv.reader(f)
for data in reader:
    print(data)
    Exist_PDB_list.append(data[0])
f.close()

# Downloaded  csv file that we need to read for data
File_to_read=open('rcsb_pdb_PNN_data1.csv','r')
#Enter the drug name to look for
a = DatasetPrep(mainFile,File_to_read,Exist_PDB_list,Drug_Name)

#Returns new list of PDB files that were not present in previous file
New_PDB_list=a.readInputCsvFile()


print('Completed successfully')


