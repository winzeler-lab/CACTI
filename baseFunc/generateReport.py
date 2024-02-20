# Write txt files
# Karla Godinez

# DECLARATIONS
import pandas as pd


# Create txt from dictionary
def fromDict(prefix,path_output,head,dic):
	
	f = open(path_output+'/'+prefix+'_litReport.txt','w')
	f.write('\t'.join(head)+'\n')
	for key,value in dic.items():
		f.write(key+'\t'+'\t'.join(value)+'\n')
	f.close()
	
# Ceate excel from network
def excelNetNode(prefix,path_output,dfnet,dfnode):
	with pd.ExcelWriter(path_output+'/'+prefix+'_similarityNetwork.xlsx') as writer:
		dfnode.to_excel(writer,sheet_name='Nodes Details',index = False)
		dfnet.to_excel(writer,sheet_name='Similarity Network',index = False)

# Create txt from network
def fromNodes(prefix,path_output,head,nodes):
	f = open(path_output+'/'+prefix+'_similarityNodeReport.txt','w')
	f.write('\t'.join(head)+'\n')
	for item in nodes:
		f.write('\t'.join(item)+'\n')
	f.close()


# Create txt from list compounds
def writeCmpds(prefix,path_output,cmpds):
	f = open(path_output+'/'+prefix+'_dbSimilarCmpd.txt','w')
	f.write('Smiles\n')
	for item in cmpds:
		f.write(item+'\n')
	f.close()











