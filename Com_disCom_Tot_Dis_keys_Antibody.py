# python scripts to calculate the distinct keys,total keys, distinct common keys ,total common keys from different aa structure
# Author:Tarikul Islam Milon
# Created on: 22/06/2023

import csv
import time
import pandas as pd
import os

#4 changes needed

def Common_key_find(key_set, *args):  # Finds common keys
    counter = 0
    for i in key_set:
        if counter == 0:
            Com_key_set = i
        else:
            Com_key_set.intersection_update(i)
        counter += 1
    return Com_key_set


def Find_key_info_all_aa(Folder_path,antibody_):  # Finds different parameters for each aa (i.e Total_keys,Distinct_keys etc)
    Distict_key_set_ = []
    All_key_dict_ = []
    aa_dis_tot_keys_dict_ = {}
    
    #reads sample details file
    File=open(f'{Folder_path}/sample_details_Corrected.csv','r') #3.Change here
    reader=csv.reader(File)
    next(reader)
    for row in reader:
        antiB_region=row[1] #4.CDRH3=row[1], CDRL3=row[2]

        Path = f'{Folder_path}/key_{row[3]}_{antiB_region}_{antibody_}.keys_theta29_dist17'
        if os.path.exists(Path):
            File_ID = f'{row[3]}_{antiB_region}_{antibody_}'
            Total_keys = 0
            Distinct_keys = 0
            Key_set = set({})
            Key_dict = {}
            # Change extension
            InputKeyFile = open(f'{Folder_path}/key_{row[3]}_{antiB_region}_{antibody_}.keys_theta29_dist17', 'r')

            lines = InputKeyFile.readlines()
            lines.remove(lines[0])  # Remove the title

            for line in lines:
                Total_keys += int(line.split()[1])
                Key_set.add(line.split()[0])
                Key_dict[line.split()[0]] = line.split()[1]
                Distinct_keys += 1

            Distict_key_set_.append(Key_set)
            All_key_dict_.append(Key_dict)
            aa_dis_tot_keys_dict_[File_ID] = [Distinct_keys, Total_keys]
            InputKeyFile.close()

    return Distict_key_set_, All_key_dict_, aa_dis_tot_keys_dict_


start=time.time()
if __name__=="__main__":

    # Contains output file from previous code (Input for this code)
    Output_Folder_Path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/' #1.Change here
    antibody='CDRH3'  #2.Change here
    OutputFile1 = open(f'{Output_Folder_Path}/{antibody}_Distinct_Total_distinctCommon_totalCommonKeys.txt', 'w')
    # OutputFile=open(f'aa_{AminoAcidName}_Distinct_Total_distinctCommon_totalCommon_Keys.txt','w')
    OutputFile1.write(f'aa\t\t\tDistinct_keys\tTotal_keys\tDistinct_common\tTotal_common\n')
    OutputFile2 = open(f'{Output_Folder_Path}/{antibody}_Common_key_Names.txt', 'w')
    OutputFile2.write(f'Common_keys\n')

    Distict_key_set, All_key_dict, aa_dis_tot_keys_dict = Find_key_info_all_aa(Output_Folder_Path,antibody)  # Finds all key information for a particular aa

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
        OutputFile1.write(f'{aa}\t{aa_dis_tot_keys_dict[aa][0]}\t\t{aa_dis_tot_keys_dict[aa][1]}\t\t{len(Common_key_set)}\t\t{Total_com_key}\n')
        aa_counter += 1
    OutputFile1.close()
    OutputFile2.close()

    for i in All_key_dict:
        print(i)

End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))

