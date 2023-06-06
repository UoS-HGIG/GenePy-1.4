import numpy as np
import sys

db = ([line.rstrip().split('\t') for line in open(sys.argv[1],'r')])
print ('imported')
db = np.array(db)
geni = np.unique(db[:,2])
samples = np.unique(db[:,0])

M = np.zeros((samples.shape[0],geni.shape[0]))
M[:] = np.nan

dic_geni = {}
dic_pz = {}

for i,n in enumerate(geni):
    dic_geni[n]=i
for i,n in enumerate(samples):
    dic_pz[n]=i
print ('dic generated')
for x in db:
	i = dic_pz[x[0]]
	j = dic_geni[x[2]]
	M[i,j] = float(x[1])
print('matrix built')
samples = samples.reshape((samples.shape[0],1))
geni = geni.reshape((1, geni.shape[0]))

geni = np.insert(geni, 0, 'Samid')

M = np.hstack((samples, M))
M = np.vstack((geni, M))

np.savetxt(sys.argv[2], M, delimiter='\t', fmt='%s')
print('matrix exported')
