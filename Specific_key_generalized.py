#Python scripts to find the specific keys from aa structure and the key percentage in each group
#Author: Tarikul Islam Milon
#Modified on:04/13/2023
#aa=amino acid
#Creates group information, each aa as a group
#It can calculates one specific keys at a time and a set of keys

import csv
import pandas as pd
import time
import itertools

def Create_group_info(aa,Folder_path,Group_,*args): #Creates group information and calculates key (%) in each group
    key_dict_list=[]
    key_percentage = {}
    Specfic_key_data={}

    #df_sample_details= pd.read_csv(f'{Folder_path}/sample_details_{aa}.csv')
    #aa_list = df_sample_details['protein'].to_list()
    aa_list=Group_[aa]
    counter=0
    Total_key = 0
    for aa_data in aa_list:
        Key_dict={}
        InputKeyFile = open(f'{Folder_path}/key_{aa_data}.keys_theta29_dist17', 'r')  # Enter the path of the Key file
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
        if key_percentage[i]==100: # Change it here to calculate specfic keys for a particular (%)
            Specfic_key_data[i]=key_percentage[i]

    return key_dict_list,key_percentage,Specfic_key_data

def Specific_key_cal(*args): #For one specific key at a time
    Counter1=0
    for data1 in args[1]:
        Spec_key_count=1
        for Spec_key in args[1][data1].keys():
            Counter2=0
            Specific_key_check = 0
            for data2 in args[0]:
                if Counter1!=Counter2:
                    if Spec_key in args[0][data2]:
                        if args[0][data2][Spec_key]>=0: #Change here for other group (%)
                            Specific_key_check=1
                Counter2+=1

            if Specific_key_check == 0:
                Counter3=0
                for data3 in args[0]:
                    if Counter1!=Counter3:
                        OutputFile2.write(f'{args[2][Counter1]}\t{Spec_key}\t{args[1][data1][Spec_key]}\t{Spec_key_count}')
                        if Spec_key in args[0][data3]:
                            OutputFile2.write(f'\t\t{args[0][data3][Spec_key]}\n')
                        else:
                            OutputFile2.write('\t\tNone\n')
                    Counter3+=1
                Spec_key_count+=1
        Counter1+=1


def Specific_key_cal_set(*args): #For a set of specific key at a time
    Counter1 = 0
    for data1 in args[1]: #Each 100% keys in the group
        Spec_key_count = 1
        # Makes set from the 100% keys
        Set=itertools.combinations(args[1][data1].keys(),args[3])
        for i in Set: #Checks which set is specific
            Counter2 = 0
            Specific_key_check = 0
            for group in args[0]: #Check other groups
                if Counter1!=Counter2:
                    for data2 in args[0][group]: #Each PDB in the group
                        if set(i).issubset(set(data2)):
                            Specific_key_check = 1
                Counter2+=1
            if Specific_key_check==0:
                OutputFile3.write(f"{args[2][Counter1]}\t{i}\t{Spec_key_count}\n")
                Spec_key_count+=1
        Counter1+=1

def Group_info_sample_details(*args):
    Group_dict={}
    df2=pd.read_csv(f'{args[0]}/tyr_group_generalised_s_h_l_mix5.csv')
    Protein_list = df2['protein'].to_list()
    Group_list=df2['group'].to_list()
    for i in range(len(Protein_list)):
        if Group_list[i] not in Group_dict:
            Group_dict[Group_list[i]]=[]
            Group_dict[Group_list[i]].append(Protein_list[i])
        else:
            Group_dict[Group_list[i]].append(Protein_list[i])
    return Group_dict



start=time.time()
if __name__=="__main__":
    #Group info
    #Group_information = ['GLN', 'TYR'] #amino acids
    Group_information = ['A', 'B'] #Divided into groups
    aa_group_info = {}
    aa_group_key_percentage = {}
    aa_group_spec_key_data={}
    Set_size=2 #change here for different set size

    OutputFile1 = open('Key(%)_each_group.txt', 'w')
    OutputFile1.write(f'Key\t\tkey (%) \tGroup_Info\tDistinct_key\tTotal_key\n')

    OutputFile2 = open('Specific_Key_Each_Group.txt', 'w')     #For one specific key at a time
    OutputFile2.write("Group\tKey\t\t(%)\tTotal_specific\tOther_groups(%)\n")


    OutputFile3 = open('Specific_Key_Each_Group_set.txt', 'w') #For a set of specific key at a time
    OutputFile3.write("Group\tKey_set\t\t\t\tTotal_specific\n")

    Input_files_folder_path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/'
    Group_dict_=Group_info_sample_details(Input_files_folder_path)


    for aa in Group_information:
        #Input_files_folder_path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/'
        #Input_files_folder_path = f'/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/intra_aa_tsr_{aa.lower()}/sample_s_h_l_mix5_original_pdb/output'
        key_dict_list_,key_percentage_,Specfic_key_data_= Create_group_info(aa, Input_files_folder_path,Group_dict_)
        aa_group_key_percentage[aa] = key_percentage_
        aa_group_info[aa] = key_dict_list_
        aa_group_spec_key_data[aa]=Specfic_key_data_

    Specific_key_cal(aa_group_key_percentage, aa_group_spec_key_data, Group_information) #For one specific key at a time
    Specific_key_cal_set(aa_group_info,aa_group_spec_key_data,Group_information,Set_size) #For a set of specific keys at a time


End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))











