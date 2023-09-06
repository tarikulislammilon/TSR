#Python scripts to make inputs for venn diagram
#Author: Tarikul Islam Milon
#Created on: 04/11/2023
#Inputs Files:aa_TYR_key_distribution_{i}.txt

import itertools
import time

def All_intersection(**kwargs):
    Counter=0
    data=list(kwargs.keys())
    for i in kwargs:
        if Counter == 0:
            Com_keys=set(kwargs[i])
        else:
            Com_keys.intersection_update(set(kwargs[i]))
        Counter += 1
    OutputFile.write(f'{data[0]}-{data[1]}-{data[2]}\t{len(Com_keys)}\n')


def Two_intersection(**kwargs):
    pairs=itertools.combinations(kwargs,2)
    for j in pairs:
        Com_key_Bet_pairs=set(kwargs[j[0]]).intersection(set(kwargs[j[1]]))
        OutputFile.write(f'{j[0]}-{j[1]}\t{len(Com_key_Bet_pairs)}\n')


start=time.time()
if __name__=='__main__':
    # Dictionary for three regions of Spike and their monoclonal antibodies
    s = {}
    h = {}
    l = {}

    List = ['s', 'h', 'l'] #Enter the list of dict for venn diagram
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

    OutputFile=open('aa_Ven_Diagram.txt','w')

    All_inter=All_intersection(s=s,h=h,l=l) #Change here for different dict
    Two_inter=Two_intersection(s=s,h=h,l=l) #Change here
    OutputFile.write(f's\t{len(s)}\nh\t{len(h)}\nl\t{len(l)}') #Change here
End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))













