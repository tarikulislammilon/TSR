#Program to calculte the specfic keys from Intra-TSR in amino acid in Lilly Datasets
#Author:Tarikul Islam Milon

import csv
Output_Folder_Path="C:/Users/C00479477/PycharmProjects/pythonProject/Output_Folder_Check"
Group_Info={} # Creating the Group
Group_Name=[]

File1=open(f'{Output_Folder_Path}/sample_details.csv','r') # Reads the Group info from a txt File
reader=csv.reader(File1)
next(reader)

for row in reader:
    if row[4] not in Group_Info:
        Group_Info[row[4]]=[]
        Group_Name.append(row[4])
        Group_Info[row[4]].append(row[0])
    else:
        Group_Info[row[4]].append(row[0])




#Creating the key percentage in Each Group
Key_Per_All_Group=[]
Spec_Key_Per_All_Group=[]
OutputFile1=open(f'{Output_Folder_Path}/Key_Percentage_Each_Group.txt','w')
OutputFile1.write(f'Key\t\tkey_per\t\tGroup_Info\tDistict_key\tTotal_key\n')
for Group in Group_Info:
    key_percentage={}
    Specfic_key_data={}
    counter=0
    Total_key=0
    for data in Group_Info[Group]:
        File2=open(f'{Output_Folder_Path}/key_{data}.keys_theta29_dist17','r')
        lines=File2.readlines()
        lines.remove(lines[0])
        for line in lines:
            Key=line.split()[0]
            Total_key+=int(line.split()[1])
            if Key not in key_percentage:
                key_percentage[Key]=1
            else:
                key_percentage[Key] += 1
        counter+=1
        File2.close()
    for i in key_percentage:
        key_percentage[i]= float("%.2f" % ((key_percentage[i]/counter)*100))
        OutputFile1.write(f'{i}\t{key_percentage[i]}\t\t{Group}\t{len(key_percentage)}\t{Total_key}\n')
        if key_percentage[i]==100: # Change it here to calculate specfic keys for a particular percent
            Specfic_key_data[i]=key_percentage[i]
    Key_Per_All_Group.append(key_percentage)
    Spec_Key_Per_All_Group.append(Specfic_key_data)

OutputFile1.close()
File1.close()



#Fining the Specific keys

OutputFile2=open(f'{Output_Folder_Path}/Specific_Key_Each_Group.txt','w')
OutputFile2.write("Group_Index\tKey\tPercentage\tTotal_specific\tPercent_Other_Groups\n")
for i in range(len(Spec_Key_Per_All_Group)):
    counter = 1
    for data in Spec_Key_Per_All_Group[i]:
        Specific_key_check = 0
        for j in range(len(Key_Per_All_Group)):
            if i != j:
                #if data in Key_Per_All_Group[j]:
                if data in Key_Per_All_Group[j]:
                    if Key_Per_All_Group[j][data] >= 0:  # Change it here to calculate specfic keys for a particular percent
                        Specific_key_check = 1

        if Specific_key_check ==0:

            for j in range(len(Key_Per_All_Group)):
                if i != j:
                    OutputFile2.write(f"{Group_Name[i]}\t{data}\t{Spec_Key_Per_All_Group[i][data]}\t{counter}")
                    if data in Key_Per_All_Group[j]:
                        OutputFile2.write(f"{Key_Per_All_Group[j][data]}\n")
                    else:
                        OutputFile2.write("\tNot Present\n")
            counter+=1




OutputFile2.close()
print("Completed_Successfully")









