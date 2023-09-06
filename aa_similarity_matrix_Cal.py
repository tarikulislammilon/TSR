#Python scripts to obtain normalized similarity matrix
#Author: Tarikul Islam Milon

import numpy as np
import pandas as pd
from sinkhorn_knopp import sinkhorn_knopp as skp
sk = skp.SinkhornKnopp()

Non_diagonal_values={}
diagonal_values={}
aa_sequence=[]

f=open('Average_Similarity_between_amino_acid_pairs.txt','r')
lines=f.readlines()
lines.remove(lines[0])
Counter1=0
Counter2=0
for line in lines:
    aa_pair=line.split()[0].split('-')
    aa1=aa_pair[0]
    aa2=aa_pair[1]
    Similarity=line.split()[1]
    if aa1!=aa2:
        Non_diagonal_values[Counter1] = Similarity
        Counter1 += 1
    else:
        diagonal_values[Counter2]=Similarity
        aa_sequence.append(aa1)
        Counter2+=1

B=np.zeros((20,20))
l=0
k=0
for i in range(len(B)):
    for j in range(l,len(B)):
        if i!=j:
            B[i, j] = Non_diagonal_values[k]
            B[j, i] = Non_diagonal_values[k]
            k += 1
    B[i,i]=diagonal_values[l]
    l+=1

df1 = pd.DataFrame(B,index=aa_sequence,columns=aa_sequence)
df1.to_csv("Data_before_normalize.csv")
#after normalization
B_ds = sk.fit(B)
B_ds=B_ds/np.min(B_ds)

df2 = pd.DataFrame(B_ds,index=aa_sequence,columns=aa_sequence)
df2.to_csv("Data_after_normalize.csv")
print('Completed Successfully')
