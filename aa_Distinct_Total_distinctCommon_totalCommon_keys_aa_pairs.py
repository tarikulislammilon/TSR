# python scripts to calculate the distinct keys,total keys, distinct common keys ,total common keys from different pairs aa structure
# Author:Tarikul Islam Milon
# Created on: 04/06/2023
# Two input required: 1) sample details and 2) key File 3) Output Folder path
# aa = Amino acids

import time
import pandas as pd
import os
import itertools

def Common_key_find(key_set, *args):  # Finds common keys
    counter = 0
    for i in key_set:
        if counter == 0:
            Com_key_set = i
        else:
            Com_key_set.intersection_update(i)
        counter += 1
    return Com_key_set

def Find_key_info_all_aa(aa_set_data,*args):  # Finds different parameters for each aa (i.e Total_keys,Distinct_keys etc)
    Distict_key_set_ = []
    All_key_dict_ = []
    aa_dis_tot_keys_dict_ = {}
    Total_keys = 0
    Distinct_keys = 0
    Key_set = set({})
    Key_dict = {}

    for i in aa_set_data:
        AminoAcidName = i
        #Contains output file from previous code (Input for this code)
        Output_Folder_Path = f'/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/intra_aa_tsr_{AminoAcidName.lower()}/sample_s_h_l_mix5_original_pdb/output'
        df = pd.read_csv(f'{Output_Folder_Path}/sample_details_{AminoAcidName}.csv')  # Input File (1)
        aa_list = df['protein'].to_list()
        for aa in aa_list:
            Total_keys = 0
            Distinct_keys = 0
            Key_set = set({})
            Key_dict = {}

            InputKeyFile = open(f'{Output_Folder_Path}/key_{aa}.keys_theta29_dist17', 'r')  # Enter the path of the Key file
            lines = InputKeyFile.readlines()
            lines.remove(lines[0])  # Remove the title
            for line in lines:
                Total_keys += int(line.split()[1])
                Key_set.add(line.split()[0])
                Key_dict[line.split()[0]] = line.split()[1]
                Distinct_keys += 1

            Distict_key_set_.append(Key_set)
            All_key_dict_.append(Key_dict)
            aa_dis_tot_keys_dict_[aa] = [Distinct_keys, Total_keys]
            InputKeyFile.close()

    return Distict_key_set_, All_key_dict_, aa_dis_tot_keys_dict_


start=time.time()
if __name__=="__main__":

    Amino_acid_list = ['TYR','GLN'] #Change the number of aa here
    aa_set=itertools.combinations(Amino_acid_list,2) # For 2 pairs of aa

    for set_data in aa_set:
        OutputFile1 = open(f'aa_{set_data[0]}_{set_data[1]}_Distinct_Total_distinctCommon_totalCommonKeys.txt', 'w')
        OutputFile1.write(f'aa\t\t\tDistinct_keys\tTotal_keys\tDistinct_common\tTotal_common\n')
        OutputFile2 = open(f'aa_{set_data[0]}_{set_data[1]}_Common_key_Names.txt', 'w')
        OutputFile2.write(f'Common_keys\n')

        Distict_key_set, All_key_dict, aa_dis_tot_keys_dict = Find_key_info_all_aa(set_data) # Finds all key information for a particular aa
        Common_key_set = Common_key_find(Distict_key_set)  # Finding the common keys

        # Writting to output file
        for com_key in Common_key_set:
            OutputFile2.write(f'{com_key}\n')

        aa_counter = 0
        for aa in aa_dis_tot_keys_dict:
            Total_com_key = 0
            if len(Common_key_set) != 0:
                for Com_key in Common_key_set:
                    Total_com_key += int(All_key_dict[aa_counter][Com_key])
            OutputFile1.write( f'{aa}\t{aa_dis_tot_keys_dict[aa][0]}\t\t{aa_dis_tot_keys_dict[aa][1]}\t\t{len(Common_key_set)}\t\t{Total_com_key}\n')
            aa_counter += 1
        OutputFile1.close()
        OutputFile2.close()

End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))

