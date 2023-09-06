import csv
from csv import writer

"""
file1=open('sample_details_nih_dataset_mix1_data_deleted.csv','r')
outputFile=f'sample_details_478_017_STI_nih_datasets.csv'
reader=csv.reader(file1)
next(reader)
drug_list=['478','STI','17']
for row in file1:
    drug=row.split(',')[2].strip()
    if drug in drug_list:
        All_data = row.split(',')
        with open(f'{outputFile}', 'a') as f_object:
            writer_object = writer(f_object)
            writer_object.writerow(All_data)
            f_object.close()
print('Completed Successfully')

"""

dict1={'100':2,'300':3,'400':4,'-100':1,'-300':5,'-500':6,'400':1,'-400':5}
dict1_arr=list(dict1.keys())

sum=0
Total_sum=0

for i in range(0,len(dict1_arr)):
    for j in range(i+1,len(dict1_arr)):
        if abs(int(dict1_arr[i]))==abs(int(dict1_arr[j])):
            sum+=2
            Total_sum+=dict1[dict1_arr[i]]+dict1[dict1_arr[j]]

print('{:.2f}'.format((sum/len(dict1)*100)))
print(Total_sum)



