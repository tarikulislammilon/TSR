# Python scripts for preparing the dataset for drugs
# Author: Tarikul Islam Milon
# Created On:
import csv
from csv import writer
import os
import urllib.request
import Bio.PDB
from Bio.PDB import PDBParser


# 3 change required

class checkCorrectChain_Gene:
    def __init__(self):
        pass

    @staticmethod
    def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
        pdbfn = pdbcode + ".pdb"
        url = downloadurl + pdbfn
        outfnm = os.path.join(datadir, pdbfn)
        try:
            urllib.request.urlretrieve(url, outfnm)
            return outfnm
        except Exception as err:
            return None

    @staticmethod
    def chainCheck(pdbcode, datadir, drug_):
        '''

        :param pdbcode: PDB file
        :param datadir: directory that contains PDB file
        :drug: drug name
        :return: Correct chain and other information of drug(resID,drugName etc.)
        change the condition for different drug
        '''

        drug_info = []
        path = f'{datadir}/{pdbcode}.pdb'
        p = Bio.PDB.PDBParser()
        Structure = p.get_structure('PrimaryStructureChain', path)
        model = Structure[0]
        for chain in model:
            for residue in chain:
                drug_Identifier = str(residue)[17:18].strip()
                resName = str(residue)[9:12].strip()
                if drug_Identifier == 'H' and resName == drug_:  # Enter the drugs name #changed
                    numeric_filter = filter(str.isdigit, str(residue.id))
                    Res_Id = "".join(numeric_filter)
                    data = f'{pdbcode}_{resName}_{Res_Id}_{chain.id}'
                    drug_info.append(data)
        return drug_info

    @staticmethod
    def geneCheck(Input_data_):
        gene_info = {}
        for row in Input_data_:
            if len(row[0]) != 0:
                Current_PDB = row[0]
            if len(row[2]) != 0:
                Current_Chain = row[2]
            if len(row[0]) != 0 and len(row[2]) != 0 and len(row[3]) == 0:
                Current_gene = 'None'
                # data=f'{Current_PDB}_{Current_Chain}_None'
            elif len(row[3]) != 0:
                Current_gene = row[3]

            data = f'{Current_PDB}_{Current_Chain}'

            gene_info[data] = Current_gene
        return gene_info


class datasetPrep:
    """
    prepares the dataset for drugs with the chain information of protein to which the drugs binds
    Reads the downloaded .csv file from rcsb database
    adds new drugs (PDB file of that drug) to main output file
    """

    def __init__(self, *args):
        self.outputFile=args[0]
        self.Input_data = args[1]
        self.InputFile2=args[2]


    def readInputCsvFile(self):

        Gene_info = checkCorrectChain_Gene.geneCheck(Input_data)
        PDB_datadir= '/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/Drug/protein_drug/dataset_prep_rxt/PDB_datadir' #2.Change here
        f = open(InputFile2, 'r')
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            PDB_ID=row[0]
            Drug=row[6]
            checkCorrectChain_Gene.download_pdb(PDB_ID, PDB_datadir) #download the PDB file first
            Path=Path=f'{PDB_datadir}/{PDB_ID}.pdb'
            if os.path.exists(Path):
                drug_Info_List = checkCorrectChain_Gene.chainCheck(PDB_ID, PDB_datadir, Drug)

                for i in drug_Info_List:
                    List = i.split('_')
                    prot_Chain = f'{List[0]}_{List[3]}'
                    List.append(Gene_info[prot_Chain])
                    with open(f'{outputFile}', 'a') as f_object:
                        writer_object = writer(f_object)
                        writer_object.writerow(List)
                        f_object.close()
            else:
                print(f'{PDB_ID} can not be downloaded')

outputFile=f'Gene_information_NIH_datasets.csv'

# downloaded csv file from rcsb database
InputFile_ = f'Sample_details_NIH_datasets_1_2500.csv'
InputFile2=f'sample_details_pdbbind_mix2.csv'
File_to_read = open(InputFile_, 'r')
reader2 = csv.reader(File_to_read)
Input_data = []
for row in reader2:
    Input_data.append(row)

a = datasetPrep(outputFile,Input_data,InputFile2)
a.readInputCsvFile()


print('Completed Successfully')
