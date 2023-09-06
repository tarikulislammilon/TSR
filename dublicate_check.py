import sys
import openeye
import csv

f=open("sample_details_nih_dataset_mix1.csv",'r')
reader=csv.reader(f)
next(reader)
data=[]
for row in reader:
    data1=f'{row[0]}_{row[1]}_{row[2]}_{row[3]}_{row[4]}_{row[5]}_{row[6]}_{row[7]}_{row[8]}_{row[9]}_{row[10]}'
    if data1 not in data:
        data.append(data1)
    else:
        print(data1)


