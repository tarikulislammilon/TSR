#Program to calculate the specific keys and common keys from Intra-TSR Amino acid keys_Spike_Ace2 (modified)
#Author: Tarikul Islam milon
#Requires two input : 1)sample_details.csv and 2) key file 

import csv
Output_Folder_Path='/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/sample_lilly_mix5_s_h_2mol/output'
#Temp=open(f'{Output_Folder_Path}/Remove_PDB.txt','w')
Summary_f=open(f"{Output_Folder_Path}/Summary_DisTot_Keys_lilly_mix5_s_h_2mol_AllAtom.txt",'w')
Summary_f.write('Residue-ID\t   Distinct\t\tTotal\n')

Group_Info={}
Group_Name=[]
Distict_Key_sets=[]
Sample_List=[]
Key_Dict=[]

File1=open(f'{Output_Folder_Path}/sample_details_lilly_mix5_s_h_2mol.csv','r')
reader=csv.reader(File1)
next(reader)
for row in reader:
    Sample_Id=row[0]
    Chain1=row[1]
    Chain2=row[2]
    if row[5] not in Group_Info: #start with column 0
        Group_Info[row[5]] = []
        Group_Name.append(row[5])
        Group_Info[row[5]].append(row[0])
    else:
        Group_Info[row[5]].append(row[0])

    File2 = open(f'{Output_Folder_Path}/key_{Sample_Id}_{Chain1}_{Chain2}_H.keys_theta29_dist17','r')# Change to L for Light Chain
    data1 = File2.readlines()
    data1.remove(data1[0])
    Total_Keys = 0
    distinct_keys = 0
    Key_set = set({})
    key_dict = {}

    for value in data1:
        Total_Keys += int(value.split()[1])
        Key_set.add(value.split()[0])
        key_dict[value.split()[0]] = value.split()[1]
        distinct_keys += 1

    Distict_Key_sets.append(Key_set)
    # Temp.write(line.strip()+'\n')
    Sample_List.append(Sample_Id)
    Key_Dict.append(key_dict)
    Summary_f.write(f"{Sample_Id}\t\t{distinct_keys}\t\t{Total_Keys}\n")
    File2.close()

File1.close()
#Finding the common keys

counter=0
for i in Distict_Key_sets:
    if counter==0:
        Common_key_set=i
    else:
        Common_key_set.intersection_update(i)

    counter+=1

#Writting the Commom keys output to a Output File
Output_Com_Keys=open(f"{Output_Folder_Path}/Common_keys_lilly_mix5_s_h_2mol_AllAtom.txt",'w')
Output_Com_Summary_Keys=open(f"{Output_Folder_Path}/Summary_Common_keys_lilly_mix5_s_h_2mol_AllAtom.txt",'w')

Output_Com_Keys.write('File_Name\tcommon_Keys\t  keysFreq\n')
Output_Com_Summary_Keys.write('File_Name\t Total_common_Keys\t  Sum_keysFreq\n')
for i in range(len(Sample_List)):
    sum=0
    for j in Common_key_set:
        sum+=int(Key_Dict[i][j])
        Output_Com_Keys.write(f'{Sample_List[i]}\t\t{j}\t\t{Key_Dict[i][j]}\n')
    Output_Com_Summary_Keys.write(f'{Sample_List[i]}\t\t{len(Common_key_set)}\t\t{sum}\n')

print('Completed_Successfully')


