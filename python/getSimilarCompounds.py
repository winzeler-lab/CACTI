# Module 2: Close analog search, and similarity between found set

# DECLARATIONS
from config import *
import getMoleculeDescriptor as mol
import parseInputFunc as pi
import parseJsonFunc as pj
import generateReport as report
import generateFigure as figure

import pandas as pd
import networkx as nx
import numpy as np
from rdkit import DataStructs

api_ncbi = ''
keywords = ''

def getSimPair(ids,fps,cutoff):
	pairs = []
	for i in range(len(fps)):
		sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[i+1:])
		for j in range(len(sims)):
			if sims[j] >= cutoff:
				pairs.append([ids[i],ids[i+j+1],sims[j]])
	return(pairs)

# Main structure
def getNetwork(type,prefix,infile,path_output,sep,treshold_simil,sim_conf):
	msg = 'Processing similarity network\n'

	# Check input file validity
	try:
		f=open(infile,'r')
		header = f.readline()
		reads = f.readlines()
		f.close()
	except:
		msg += 'Error: Reading ' + infile + ' file\n'
		return(False,msg)


	# Obtain lists for smiles - names from input file
	if type == 2:
		cmp_name = pi.getCmpName(reads,sep)
		cmp_smi = pi.getCmpSmi(reads,sep)
	else:
		cmp_name = ['unknown_cmp'+str(i+1) for i in range(len(reads))]
		cmp_smi = [r.strip() for r in reads]

	# Transform SMILES for clustering among them
	can_smi = mol.getCanonical(cmp_smi)
	msset = mol.getMS(can_smi)
	fpset = mol.getFingerprint(msset)

	# Structure to save pairs
	# Source = compound 1, Target = compound 2, Simil = similarity percentage between source and target molecules
	df_head = ['Source', 'Target', 'Simil']
	df = pd.DataFrame(columns = df_head)
	db_cmp_smi = {}

	# Similarity based on dataset given
	pairs = getSimPair(cmp_name,fpset,treshold_simil)

	if pairs:
		for p in pairs:
			s = pd.Series([p[0],p[1],str(float(p[2])*100)],index=df_head)
			df = df.append(s,ignore_index=True)
	pairs_all = pairs #structure to save pairs of similar compounds, in case database search is selected

	# Find database close analogs, and calculate/obtain similarities
	if ('simnetdb' in sim_conf):
		db_cmp = {}
		for i in range(len(can_smi)):

			# Search in CHEMBL
			pairs = pj.getChemblSimilPerc(can_smi[i],int((treshold_simil*100)))

			if pairs:
				for p in pairs:
					s = pd.Series([cmp_name[i],p[0],p[2]],index=df_head)
					df = df.append(s,ignore_index=True)
					db_cmp_smi[p[0]] = p[1] #Save found compound's smiles
					pairs_all.append([cmp_name[i],p[0],p[2]])

			# Search in PUBCHEM
			similar_struc = pj.getPubChemSubstructure(can_smi[i],int((treshold_simil*100)))
			if similar_struc:
				for s in similar_struc:
					simil_percentage = mol.similPercentSmiles(cmp_smi[i],s[1]) #Fingerprint similarity checkup
					#simil_percentage = treshold_simil*100  # Fingerprint similarity checkup
					if (simil_percentage >= treshold_simil*100 and simil_percentage < 100):
						pair = pd.Series([cmp_name[i],'CID '+str(s[0]),simil_percentage],index=df_head)
						df = df.append(pair,ignore_index=True)
						db_cmp_smi['CID '+str(s[0])] = pj.getCIDSmiles(str(s[0]))
						pairs_all.append([cmp_name[i],'CID '+str(s[0]),simil_percentage])

			# Search in BINDINGDB
			similar_bindingdb = pj.getBindingDBSimil(can_smi[i], treshold_simil)
			if similar_bindingdb:
				for s in similar_bindingdb:
					simil_percentage = mol.similPercentSmiles(cmp_smi[i], s[1])  # Fingerprint similarity checkup
					if (simil_percentage >= treshold_simil * 100 and simil_percentage < 100):
						pair = pd.Series([cmp_name[i], str(s[0]), simil_percentage], index=df_head)
						df = df.append(pair, ignore_index=True)
						db_cmp_smi[str(s[0])] = s[1]
						pairs_all.append([cmp_name[i], str(s[0]), simil_percentage])
			sleep(4)  # Request time lag controller

		# Find references (Module 1) for close analogs, if selected
		if ('simnetlit' in sim_conf and pairs_all):
			import getLiterature as lit
			try:
				report.writeCmpds(prefix,path_output,db_cmp_smi.values())
				msg += 'Txt with similar compound SMILES from search created\n'
				msg += 'Procesing literature of similar compounds\n'
				lit.getLit(1,prefix+'_dbSimilarCmpd',path_output+'/'+prefix+'_dbSimilarCmpd.txt',path_output,sep,api_ncbi,keywords)
				msg += 'Literature of similar compounds completed\n'
			except:
				msg += 'Unable to get literature for similar compounds\n'

	# If similarity pairs were found, generate report
	if pairs_all:
		# Obtain nodes details
		dfn_head = ['Compound','Smiles','Canonical smiles','Type']
		df_nodes = pd.DataFrame(columns = dfn_head)

		## Compounds from dataset
		for i in range(len(cmp_name)):
			s = pd.Series([cmp_name[i], cmp_smi[i], cmp_smi[i], 'dataset'], index=dfn_head)
			df_nodes = df_nodes.append(s,ignore_index=True)
		## Compounds from database
		for key,val in db_cmp_smi.items():
			s = pd.Series([key,'',val,'database'],index=dfn_head)
			df_nodes = df_nodes.append(s,ignore_index=True)

		# Generate report
		df['Simil'] = df['Simil'].astype(float)
		df.sort_values(by=['Simil'])

		report.excelNetNode(prefix,path_output,df,df_nodes)
	else:
		msg += 'No network to be generated\n'

	msg += 'Success: Complete similarity network\n'
	return (True,msg)











