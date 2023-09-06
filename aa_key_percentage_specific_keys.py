#Python scripts to find the specific keys from aa structure and the key percentage in each group
#Author: Tarikul Islam Milon
#Created on:
#aa=amino acid
#Creating the group information, each aa as a group

import csv
import pandas as pd
import os

def Create_group_info(aa,Folder_path,*args): #Creates group information and calculates key (%) in each group
    key_dict_list=[]
    key_percentage = {}
    Specfic_key_data={}


    df_sample_details= pd.read_csv(f'{Folder_path}/sample_details_{aa}.csv')
    aa_list = df_sample_details['protein'].to_list()
    counter=0
    Total_key = 0
    for aa_data in aa_list:
        Key_dict={}
        Path=f'{Folder_path}/key_{aa_data}.keys_theta29_dist17'
        if os.path.exists(Path):
            print(aa_data)
            InputKeyFile = open(f'{Folder_path}/key_{aa_data}.keys_theta29_dist17', 'r')  # Enter the path of the Key file
            # InputKeyFile = open(f'key_{aa}.keys_theta29_dist17', 'r')
            lines = InputKeyFile.readlines()
            lines.remove(lines[0])  # Remove the title
            for line in lines:
                Key = line.split()[0]
                Key_freq = line.split()[1]
                Key_dict[Key] = Key_freq
                Total_key += int(Key_freq)
                if Key not in key_percentage:
                    key_percentage[Key] = 1
                else:
                    key_percentage[Key] += 1
            counter += 1
            key_dict_list.append(Key_dict)

    for i in key_percentage:
        key_percentage[i] = float("%.2f" % ((key_percentage[i] / counter) * 100))
        OutputFile1.write(f'{i}\t{key_percentage[i]}\t\t{aa}\t\t{len(key_percentage)}\t\t{Total_key}\n')
        if key_percentage[i]==100: # Change it here to calculate specfic keys for a particular percent
            Specfic_key_data[i]=key_percentage[i]

    return key_dict_list,key_percentage,Specfic_key_data


"""
def Specific_key_cal(aa_group_key_percentage_,aa_group_spec_key_data_,amino_acids,*args): #Calculates Specific keys
    Counter1=0
    for data1 in aa_group_spec_key_data_:
        Spec_key_count=1
        for Spec_key in aa_group_spec_key_data_[data1].keys(): #Each 100% keys in the group
            Counter2 = 0
            Specific_key_check=0
            for data2 in aa_group_key_percentage_: #Search for other groups
                if Counter1!=Counter2:
                    if Spec_key in aa_group_key_percentage_[data2]:
                        if aa_group_key_percentage_[data2][Spec_key] >= 0:
                            Specific_key_check = 1
                Counter2+=1

            if Specific_key_check==0:
                Counter3 = 0
                for data2 in aa_group_key_percentage_:
                    if Counter3!=Counter1:
                        OutputFile2.write(f"{amino_acids[Counter1]}\t{Spec_key}\t{aa_group_spec_key_data_[data1][Spec_key]}\t{Spec_key_count}")
                        if Spec_key in aa_group_key_percentage_[data2]:
                            OutputFile2.write(f"\t\t{aa_group_key_percentage_[data2][Spec_key]}\n")
                        else:
                            OutputFile2.write("\t\tNot Present\n")
                    Counter3+=1
                Spec_key_count+=1
        Counter1+=1
"""
if __name__=="__main__":
    Amino_acid_list = ['GLN', 'TYR']
    aa_group_info = {}
    aa_group_key_percentage = {}
    aa_group_spec_key_data={}

    OutputFile1 = open('Key_percentage_each_group.txt', 'w')
    OutputFile1.write(f'Key\t\tkey (%) \tGroup_Info\tDistinct_key\tTotal_key\n')

    OutputFile2 = open('Specific_Key_Each_Group.txt', 'w')
    OutputFile2.write("Group\tKey\t\t(%)\tTotal_specific\tOther_groups(%)\n")

    for aa in Amino_acid_list:
        Input_files_folder_path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/'
        # Input_files_folder_path = f'/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/intra_aa_tsr_{aa.lower()}/sample_s_h_l_mix5_original_pdb/output'
        key_dict_list_, key_percentage_ ,Specfic_key_data_= Create_group_info(aa, Input_files_folder_path)
        aa_group_key_percentage[aa] = key_percentage_
        aa_group_info[aa] = key_dict_list_
        aa_group_spec_key_data[aa]=Specfic_key_data_

    Specific_key_cal(aa_group_key_percentage,aa_group_spec_key_data,Amino_acid_list)











