# Module 1: Synonyms, literature, datasets and patent search

# DECLARATIONS
from config import *
import parseInputFunc as pi
import parseJsonFunc as pj
import getMoleculeDescriptor as mol
import generateReport as report
import pandas as pd
from qualityCheck import getStatusSMILES as getStatusSMILES
from time import sleep


# SET DEFAULTS
head_lit = ['MW','Synonym','CID','CID_PMID','CID_Patent','PMID','PMC','EBI PMID','Chembl ID','','','']

# Main structure
def getLit(type,prefix,infile,path_output,sep,api_ncbi,keywords):
	status = False
	msg = 'Processing Literature Search\n'

	# Check input file validity
	try:
		f=open(infile,'r')
		head = f.readline()
		reads = f.readlines()
		f.close()
	except:
		msg += 'Error: Reading '+infile+' \n'
		return(False,msg)

	# Obtain lists for smiles - names from input file
	if type == 2:
		smiset = []
		names = []
		for r in reads:
			tmp = r.strip().split(sep)
			names.append(tmp[0])
			smiset.append(tmp[1])
	else:
		smiset = [r.strip() for r in reads]
		names = ['unknown_'+str(i) for i in range(1,len(reads)+1)]		
	
	data = pd.read_csv(infile,header=0,sep=sep)

	# Gather information for each SMILES
	for i in range(len(smiset)):
		#print(str(i+1)+'\t'+names[i]) # To keep a counter on which SMILES have been compared, uncomment this line

		if (getStatusSMILES(smiset[i]) == True):
			
			data.loc[i, 'MW'] = mol.getMWfromSMI(smiset[i]) # Get molecular weight

			# Obtain synonyms
			if type == 1:
				data.loc[i, 'Synonym'] = pj.getSynonym('',smiset[i])
			else:
				data.loc[i, 'Synonym'] = pj.getSynonym(names[i],smiset[i])

			# Obtain PubChem Compound ID
			if type == 1:
				data.loc[i, 'CID'] = pj.getPubchemCID('',smiset[i])
			else:
				data.loc[i, 'CID'] = pj.getPubchemCID(names[i],smiset[i])

			# Obtain ChEMBL Compound ID
			data.loc[i, 'ChEMBL ID'] = pj.getChemblSmiID(smiset[i],100)


			# Search Literature
			### Pubchem CID search
			data.loc[i,'CID_PMID'] = pj.getPubmedID_CID(data.loc[i,'CID'])
			data.loc[i,'CID_Patent'] = pj.getPatent_CID(data.loc[i,'CID'])

			### NCBI database search
			if type == 1:
				data.loc[i,'PMID'] = pj.getPubmedID(data.loc[i, 'Synonym'],'pubmed',api_ncbi,keywords)
				data.loc[i,'PMC'] = pj.getPubmedID(data.loc[i, 'Synonym'],'pmc',api_ncbi,keywords)
			else:
				data.loc[i,'PMID'] = pj.getPubmedID(names[i]+'\n'+data.loc[i, 'Synonym'],'pubmed',api_ncbi,keywords)
				data.loc[i,'PMC'] = pj.getPubmedID(names[i]+'\n'+data.loc[i, 'Synonym'],'pmc',api_ncbi,keywords)

			### EMBL-EBI database search
			if type == 1:
				data.loc[i,'Europe PMC'] = pj.getEBIID(data.loc[i, 'Synonym'])
			else:
				data.loc[i,'Europe PMC'] = pj.getEBIID(names[i]+'\n'+data.loc[i, 'Synonym'])

			### ChEMBL search, executed only if an ID was identified
			if (data.loc[i, 'ChEMBL ID'] != ''):
				data.loc[i,'ChEMBL Lit'] = pj.getChemblLit(data.loc[i, 'ChEMBL ID'])
				data.loc[i,'ChEMBL Mechanism'] = pj.getChemblMech(data.loc[i, 'ChEMBL ID'])

			### BindingDB
			data.loc[i,'BindingDB'] = pj.getBindingDB(smiset[i],1.0)
			print(names[i])
			sleep(4) # Request time lag controller
			#data.to_csv(path_output+'/'+prefix+'_lit.txt',sep='\t',index=False,header=True,na_rep='') # To generate an intermediate CSV between calls, uncomment this line
		else:
			msg += smiset[i] + ' is invalid SMILES\n'



	# Generate report
	try:
		data.replace('',' ',inplace=True)
		data.to_excel(path_output+'/'+prefix+'_lit.xlsx',sheet_name='Literature search',header=True,index=False)
		#data.to_csv(path_output+'/'+prefix+'_lit.txt',sep='\t',index=False,header=True,na_rep='') # To save output as csv, uncomment this line
	except:
		msg += 'Error: Could not generate output file\n'

	msg += 'Success: Complete literature search\n'
	return(True,msg)










