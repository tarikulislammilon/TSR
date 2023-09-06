#Help code -aa count

f1=open('sample_details.csv','r')
lines=f1.readlines()
lines.remove(lines[0])
counter=0

aa_count={}
for line in lines:
    aa=line.split(',')[2].strip()
    if aa not in aa_count:
        aa_count[aa]=1
    else:
        aa_count[aa] += 1

f2=open("aa_count.txt",'w')
for i in aa_count:
    f2.write(f'{i}\t{aa_count[i]}\n')

print('completed Successfully')




