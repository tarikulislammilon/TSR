# Program to calculate the similarity between amino acid from simillarity matrix of clustering
#Authror:Tarikul Islam Milon

#aa=amino acids
# Requires two inputs 1) aa_list 2) simillarity matrix as .csv file(No header, aa index format:"7KMG_C_LYS_356_-5.45")

import itertools
#Enter the amino acids for which simillarity should be calculated
aa_list=['VAL','TYR']
aa_pairs=itertools.combinations(aa_list,2) #Creating aa pairs for simillarity calculations
aa_pairs=list(aa_pairs)
for l in aa_list:
    aa_pairs.append((l,l)) #adding the aa itself pairs

#Step1
#Keeps the position index and names of all aa
pos_dict={} #position index
aa_name_dict={} # names
f1=open('generalised.csv','r')
lines=f1.readlines()
for i in aa_list:
    Pos_Index=0
    for line in lines:
        row = line.split(',')
        aa = row[0].split('_')[2] #Reads only the aa (Change here for different aa index format)
        if aa == i:
            if aa not in pos_dict:
                pos_dict[aa]=[]
                aa_name_dict[aa]=[]
                pos_dict[aa].append(Pos_Index)
                aa_name_dict[aa].append(row[0])
            else:
                pos_dict[aa].append(Pos_Index)
                aa_name_dict[aa].append(row[0])
        Pos_Index+=1

#Step2
#Calculates simillarity between amino acids pairs
output_File1 = open('Average_Similarity_between_amino_acid_pairs.txt', 'w')
output_File1.write('aa_pair\t Similarity(%)\n')

for aa_pair_data in aa_pairs: #pass each pair
    aa_1 = aa_pair_data[0]
    aa_2 = aa_pair_data[1]
    output_File2 = open(f'Similarity_between_amino_acid_pairs_{aa_1}_{aa_2}.txt', 'w')
    output_File2.write('aa1\t\t\taa2\t\t\tSimilarity(%)\n')

    Simillarity_dict = {}
    Simill_sum = 0
    counter = 0

    for line_ in lines:
        row = line_.split(',')
        aa = row[0].split('_')[2] #Reads only the aa (Change here for different aa index format)

        if aa == aa_1: # 1st aa from the pair
            aa_1_name=row[0]
            aa_2_positions = pos_dict[aa_2] # positions of 2nd aa from the pair
            aa_2_names=aa_name_dict[aa_2]  # Names of 2nd aa from the pair
            row.remove(row[0])
            sum = 0
            aa_name_Index=0

            for j in aa_2_positions:
                sum += (1 - float(row[j]))
                if aa_1_name != aa_2_names[aa_name_Index]:
                    output_File2.write(f'{aa_1_name}\t{aa_2_names[aa_name_Index]}\t{"{:.2f}".format((1 - float(row[j])) * 100)}\n') #writes the similarity for each individual pair
                aa_name_Index += 1

            average_simillairy = sum / len(aa_2_positions)
            Simill_sum += average_simillairy

            counter += 1

    Total_Avg_Simill_each_pair = "{:.2f}".format((Simill_sum / counter)*100)
    output_File1.write(f'{aa_1}-{aa_2}\t {Total_Avg_Simill_each_pair}\n')
    output_File2.close()

output_File1.close()
print("completed Successfully")








