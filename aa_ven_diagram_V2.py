#Python scripts to make inputs for venn diagram
#Author: Tarikul Islam Milon
#Created on: 04/11/2023
#Inputs Files:aa_TYR_key_distribution_{i}.txt

import itertools
import time
import pandas as pd

def Find_key_info_all_aa(*args):
    key_percentage = {}
    counter=0
    for aa in args[0]:
        InputKeyFile = open(f'{args[1]}/key_{aa}.keys_theta29_dist17', 'r')  # Enter the path of the Key file
        lines = InputKeyFile.readlines()
        lines.remove(lines[0])  # Remove the title
        for line in lines:
            Key = line.split()[0]
            if Key not in key_percentage:
                key_percentage[Key] = 1
            else:
                key_percentage[Key] += 1
        counter += 1

    for i in key_percentage:
        key_percentage[i] = float("%.2f" % ((key_percentage[i] / counter) * 100))
        OutputFile3.write(f'{args[2]}\t{i}\t{key_percentage[i]}\n')


    return key_percentage

def All_intersection(**kwargs):
    Counter=0
    k=0
    Total_distinct_dict={}
    data=list(kwargs.keys())
    for i in kwargs:
        for values in kwargs[i]:
            if values not in Total_distinct_dict:
                Total_distinct_dict[values]=k
                k+=1
        if Counter == 0:
            Com_keys=set(kwargs[i])
        else:
            Com_keys.intersection_update(set(kwargs[i]))
        Counter += 1

    OutputFile1.write(f'{data[0]}-{data[1]}-{data[2]}\t{len(Com_keys)}\t\n')
    OutputFile1.write(f'Total\t{len(Total_distinct_dict)}\n')


def Two_intersection(**kwargs):
    pairs=itertools.combinations(kwargs,2)
    for j in pairs:
        Com_key_Bet_pairs=set(kwargs[j[0]]).intersection(set(kwargs[j[1]]))
        OutputFile1.write(f'{j[0]}-{j[1]}\t{len(Com_key_Bet_pairs)}\n')

def Single_only(Key_percent_dict_,**kwargs): #Calculates the keys that are present only in one region not in others
    Count1=0
    for i in kwargs:
        Count2=0
        temp_dict=set(kwargs[i])
        for j in kwargs:
            if Count1 != Count2:
                Com_keys=set(kwargs[i]).intersection(set(kwargs[j]))
                temp_dict=temp_dict-Com_keys

            Count2+=1
        for keys in temp_dict:
            OutputFile2.write(f'{i}\t{keys}\t{len(temp_dict)}\t{Key_percent_dict_[i][keys]}\n')
        Count1+=1

start=time.time()
if __name__=='__main__':
    # Dictionary for three regions of Spike and their monoclonal antibodies
    AminoAcidName='TYR'
    s = {}
    h = {}
    l = {}

    List = ['s', 'h', 'l'] #Enter the list of dict for venn diagram
    Key_percent_dict={}
    for i in List:
        f = open(f"aa_TYR_key_distribution_{i}.txt", 'r')
        lines = f.readlines()
        lines.remove(lines[0])
        for line in lines:
            Key = line.split()[0]
            Key_Freq = line.split()[1]
            if i == 's':
                s[Key] = Key_Freq
            elif i == 'h':
                h[Key] = Key_Freq
            else:
                l[Key] = Key_Freq

    OutputFile1=open('aa_Ven_Diagram.txt','w')
    OutputFile2=open('aa_Ven_Diagram_single_value.txt','w')
    OutputFile2.write('Region\tkeys\t\tTotal_keys\tkey(%)\n')
    OutputFile3 = open('aa_Ven_Diagram_key_(%).txt', 'w')
    OutputFile3.write('Group\tkey\t\tkey(%)\n')

    for i in List:
        Output_Folder_Path = f'/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/intra_aa_tsr_{AminoAcidName.lower()}/sample_{i}_mix5_original_pdb/output'
        #Output_Folder_Path = '/Users/c00479477/Desktop/PycharmProjects/pythonProject/'
        df = pd.read_csv(f'{Output_Folder_Path}/sample_details_{AminoAcidName}.csv')  # Input File (1)
        aa_list = df['protein'].to_list()
        key_percentage_ = Find_key_info_all_aa(aa_list, Output_Folder_Path, i)
        Key_percent_dict[i] = key_percentage_

    All_intersection(s=s,h=h,l=l) #Change here for different dict
    Two_intersection(s=s,h=h,l=l) #Change here
    Single_only(Key_percent_dict,s=s,h=h,l=l)#Change here
    OutputFile1.write(f's\t{len(s)}\nh\t{len(h)}\nl\t{len(l)}') #Change here

    OutputFile1.close()
    OutputFile2.close()
    OutputFile3.close()


End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))













