#Python scripts for key searching for triplets file
#Author:Tarikul Islam Milon
#Created on:04/19/2023
#Two inputs file required, 1)key_search_list.csv and 2) Triplets file folder path
import os
import time

start=time.time()
Key_dict={}
Space =' '
Outputfile1=open('Key_information.txt','w')
Outputfile1.write(f'PDB_ID\t\t\tResidue1\tResidue2\tResidue3\tEdge1\t{Space*3}Edge2{Space*6}Edge3\t{Space}'
                  f'CoorR1\t\t\tCoorR2\t\t\tCoorR3\t\t\tTheta\t{Space*2}Max_dist{Space*2}d_3{Space*7}key3D\n')
Outputfile2=open('Summary_Key_information.txt','w')
Outputfile2.write(f'PDB_ID\t\t\tKey\t\tFreq\n')

f1=open('key_search_list.csv','r') # File contains key information
lines=f1.readlines()
lines.remove(lines[0]) #Removing title
Counter=0
for line in lines:
    Key=line.split(',')[1].strip()
    Key_dict[Key]=Counter
    Counter+=1


for file in os.listdir("/Users/c00479477/Desktop/PycharmProjects/pythonProject/"): #Enter the folder path here
    if file.endswith(".triplets_theta29_dist17"):
        f2=open(f'/Users/c00479477/Desktop/PycharmProjects/pythonProject/{file}','r')
        data=f2.readlines() #Removing title
        data.remove(data[0])
        PDB_ID = file.replace('.triplets_theta29_dist17', '')
        Key_freq_dict={}
        for line in data:
            Info=line.split()
            Only_Key_info=Info[len(Info)-1]
            for Keys in Key_dict:
                if Keys==Only_Key_info:
                    Outputfile1.write(f'{PDB_ID:<10}\t{Info[0]:<10}\t{Info[1]:<10}\t{Info[2]:<10}\t{Info[3]:<10}{Space}{Info[4]:<10}'
                                      f'{Space}{Info[5]:<10}{Space}{Info[6]:<10}\t{Info[7]:<10}\t{Info[8]:<10}\t{Info[9]:<10}{Info[10]:<10}'
                                      f'{Info[11]:<10}{Info[12]:<10}\n')
                    if Keys not in Key_freq_dict:
                        Key_freq_dict[Keys] = 1
                    else:
                        Key_freq_dict[Keys] += 1

        for key in Key_freq_dict:
            Outputfile2.write(f'{PDB_ID}\t{key}\t{Key_freq_dict[key]}\n')
Outputfile1.close()
Outputfile2.close()
End = time.time()
Total_time = End - start
print("Completed successfully\nExecution Time: {:.2f} s".format(Total_time))


