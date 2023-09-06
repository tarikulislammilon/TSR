#Program to calculate venn diagram using 2,3,4,5,6 sets
#Author: Tarikul Islam Milom
#Created on:06/07/2023

import matplotlib.pyplot as plt
from venn import venn
import csv

f=open("venn_4group_raw_data.csv",'r')
reader=csv.reader(f)
next(reader)
ChlA=[]
ChlD=[]
ChlF=[]
Pheophytin=[]
for row in  reader:
    if len(row[0])!=0:
        ChlA.append(row[0])
    if len(row[1])!=0:
        ChlD.append(row[1])
    if len(row[2])!=0:
        ChlF.append(row[2])
    if len(row[3])!=0:
        Pheophytin.append(row[3])

dataset_dict = {'ChlA':set(ChlA),'ChlD':set(ChlD),'ChlF':set(ChlF),'Pheophytin':set(Pheophytin)}
venn(dataset_dict, fmt="{percentage:.1f}%", fontsize=8, legend_loc="upper left")
plt.savefig('Venn_diagram_plot')
plt.show()