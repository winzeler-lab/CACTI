# Update reference files
# Karla Godinez

import time
import pandas as pd
import getMoleculeDescriptor as mol

def updateAntimSet(ref,upd):
	
	dic_ref = {}	
	head_ref = ref.readline().strip().split('\t')
	print(head_ref)
	for line in open('ref','r').readlines()[1:]:
		tmp = line.strip().split('\t')
		dic_ref[tmp[0]] = tmp[1]

	dic_upd = {}
	head_upd = upd.readline().strip().split('\t')
	print(head_upd)
	for line in open('upd','r').readlines()[1:]:
		tmp = line.strip().split('\t')
		dic_upd[tmp[0]] = tmp[1]

	
def removeDuplicates(df,col,outdir,name):
	#dup = data_df[data_df.duplicated(['can'],keep=False)]
	#dup.to_excel('duplicate.xlsx',header=True)
	
	#for ind in df.index:
		#print(df[col][ind])
	#dup = df.duplicated(col,keep=False)
	dup = df.drop_duplicates(subset=col,keep='first')
	
	dup.to_excel(outdir+'/'+name+'.xlsx',header=True,index=False)


	
def updateRef(infile,updater,insep,upsep,outdir):
	inf = pd.read_csv(infile,sep=insep)
	upd = pd.read_csv(updater,sep=upsep)
	name = infile.strip().split('.')[0].split('/')[-1]

	col_list = upd.columns.tolist()
	
	result = inf
	#result = inf.loc[:, col_list].fillna('')
	result.loc[:, col_list].fillna('')
	for index, row in upd.iterrows():
		if(result['Compound Name'] == row[0]).any():
			for col in col_list[2:]:
				result.at[result.index[result['Compound Name'] == row[0]],col] = row[col]

		else:
			result.append(row,ignore_index=True)
			#break
	
	result.to_csv(outdir+'/'+name+'_updated.txt',header=True,index=False,sep=insep)
	print('File updated')
	
	


def removeDuplicateFile(infile,sep,outdir):
	inf = pd.read_csv(infile,sep=sep)
	name = infile.strip().split('.')[0].split('/')[-1]

	inf['can'] = mol.getCanonical(inf['SMILES'].tolist())
	dup = inf.drop_duplicates(subset=['can'],keep='first')
	del dup['can']
	dup.to_csv(outdir+'/'+name+'_updated.txt',header=True,index=False,sep=sep)
	print('Removed ' + str(len(inf)-len(dup)) + ' duplicates')





