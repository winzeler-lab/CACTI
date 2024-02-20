# Main
# Karla Godinez

# Declarations
import os
import argparse
import time
import sys

from config import *
import getSimilarCompounds as similComp
import getMoleculeDescriptor as molDescr

# Initial variables
data_out = {}
simil_threshold = 0.80

# Separate compound name from SMILES
cmp_name = []
cmp_smi = []
header = ['Compound','SMILES','MW','Cluster']
for line in open('big_test.txt','r').readlines()[1:]:
	tmp = line.strip().split('\t')
	cmp_name.append(tmp[0])
	cmp_smi.append(tmp[1])


# Get molecule descriptors
ms = molDescr.getMS(cmp_smi) # molecule description
mw = molDescr.getMW(ms) # molecular weight
#fp = molDescr.getMorganFingerprint4(ms) # fingerprint
#fp = molDescr.getMorganBit(ms)
#print(cmp_name[6400],cmp_name[6484])
# Cluster
#clust = similComp.ClusterFps(fp,simil_threshold)
#clust = ['']*len(fp)



# Generate report
f = open('clusV.txt','w')
f.write('\t'.join(header)+'\n')
for i in range(len(cmp_name)):
	f.write(cmp_name[i]+'\t'+cmp_smi[i]+'\t'+str(mw[i])+'\n')#+'\t'+str(clust[i])+'\n')
f.close()




print('End')




