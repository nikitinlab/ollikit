
import pandas as pd


A = pd.read_csv('./mature.txt', sep='\t', header = None).to_numpy()
miRNA = {}
for i in range(len(A)//2):
    if 'Homo sapiens' in str(A[i]):
        key = str(A[i]).split(' ')[0][3:]
        seq = str(A[i+1])[2:-2]
        miRNA[key] =  seq
names = list(miRNA.keys())
print (miRNA)