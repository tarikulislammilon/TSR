#Program to calculte the specfic keys from Intra-TSR in amino acid in Lilly Datasets
#Author:Tarikul Islam Milon
#Two inputs: (i) sample_details.csv with group information; (ii) key files

import csv
import os
import itertools


Output_Folder_Path="/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/sample_lilly_mix5_s_h_2mol/output"
Group_Info={} # Creating the Group
Group_Name=[]

File1=open(f'{Output_Folder_Path}/sample_details_lilly_s_h_mix5_2mol_5A.csv','r') # Reads the Group info from a txt File and column number for group information
reader=csv.reader(File1)
next(reader)



for row in reader:
    if row[5] not in Group_Info:
        Group_Info[row[5]]=[]
        Group_Name.append(row[5])
        Group_Info[row[5]].append(f"{row[0]}_{row[1]}_{row[2]}")
    else:
        Group_Info[row[5]].append(f"{row[0]}_{row[1]}_{row[2]}")




#Creating the key percentage in Each Group
Key_Per_All_Group=[]
Spec_Key_Per_All_Group=[]
key_Dict_Total=[]
#OutputFile1=open(f'{Output_Folder_Path}/Key_Percentage_Each_Group.txt','w')
#OutputFile1.write(f'Key\t\tkey_per\t\tGroup_Info\tDistict_key\tTotal_key\n')


for Group in Group_Info:
    key_percentage={}
    Specfic_key_data={}
    counter=0
    Total_key=0
    key_Dict_Each_Group=()
    for data in Group_Info[Group]:
        File2 = open(f'{Output_Folder_Path}/key_{data}_H.keys_theta29_dist17', 'r')  # Change to L for light chain
        lines = File2.readlines()
        lines.remove(lines[0])
        Key_Dict = {}
        for line in lines:
            Key = line.split()[0]
            key_Freq=line.split()[1]
            Total_key += int(key_Freq)
            Key_Dict[Key]= key_Freq

            if Key not in key_percentage:
                key_percentage[Key] = 1
            else:
                key_percentage[Key] += 1
        counter += 1
        key_Dict_Each_Group+=(Key_Dict,)
        File2.close()
    key_Dict_Total.append(key_Dict_Each_Group)

    for i in key_percentage:
        key_percentage[i] = float("%.2f" % ((key_percentage[i] / counter) * 100))
        #OutputFile1.write(f'{i}\t{key_percentage[i]}\t\t{Group}\t{len(key_percentage)}\t{Total_key}\n')
        if key_percentage[i] == 100:  # Change it here to calculate specfic keys for a particular percent
            Specfic_key_data[i] = key_percentage[i]
    Key_Per_All_Group.append(key_percentage)
    Spec_Key_Per_All_Group.append(Specfic_key_data)



#OutputFile1.close()
File1.close()

#Fining the Specific keys
i=0 #Change it here for the group index
OutputFile2=open(f'{Output_Folder_Path}/Specific_Key_Each_Group_{i}.txt','w')
OutputFile2.write("Group_Index\tKey\t\tTotal_specific\n")

counter = 1
for data in list(itertools.combinations(Spec_Key_Per_All_Group[i], 2)): #Change it here for the set number
    Specific_key_check = 0
    for j in range(len(key_Dict_Total)):
        if i != j:
            for dict_keys in key_Dict_Total[j]:
                if dict_keys.keys() >= set(data):
                    Specific_key_check = 1

    if Specific_key_check == 0:
        OutputFile2.write(f"{Group_Name[i]}\t{data}\t{counter}\n")  # Change in here more more than 2 sets
        counter += 1

OutputFile2.close()
print("Completed_Successfully")

#print(key_Dict_Total)


