# Module 3: Chemical clustering analysis

# DECLARATIONS
from config import *
import getMoleculeDescriptor as mol
from rdkit import DataStructs
import parseInputFunc as pi
import parseJsonFunc as pj
import generateReport as report
import generateFigure as figure

import pandas as pd
import networkx as nx


# Convert cluster tuples into list of compound and cluster ID
def getClusterList(clusters,nfp):
	l = ['']*nfp
	count = 1
	for i in clusters:
		for j in i:
			l[j] = count
		count += 1
	return l


# Calculate similarity percentage between two fingerprints
def similPercentage(fp1,fp2):
	from rdkit import DataStructs
	percentage = DataStructs.FingerprintSimilarity(fp1,fp2) * 100

	return (percentage)


# Cluster with Tanimoto and Butina(lower cutoff, more similar)
def ClusterFps(fps,cutoff):
	cutoff = 1-cutoff
	from rdkit.ML.Cluster import Butina

	# Generate distance matrix
	dists = []
	nfps = len(fps)
	
	for i in range(1,nfps):
		sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
		#sims = DataStructs.BulkDiceSimilarity(fps[i],fps[:i]) # Alternative similarity metric
		dists.extend([1-x for x in sims])
	
	# Obtain clusters with distance matrix
	cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True,reordering=True)
	cs = getClusterList(cs,len(fps))
	
	return cs


# Main structure
def getPrediction(type,prefix,infile,path_output,sep,treshold_simil,infile_pred):
	msg = 'Processing in silico prediction\n'

	# Check input file validity
	try:
		f=open(infile,'r')
		f.close()
	except:
		msg += 'Error: Reading '+ infile + ' file\n'
		return(False,msg)

	
	# Create dataframe with reference (known drug-target pairs) and user input files (query file)
	## Read antimalarial or input prediction set, if alternative prediction file was provided
	### Please not that alternative separator character must be the same as the input file, if different from default (default: '\t')
	if infile_pred:
		df_antim = pd.read_csv(infile_pred, header=0, sep=sep)
	else:
		df_antim = pd.read_csv('ref/antimalarials-annotations.txt',header=0,sep='\t')

	df_antim['Type'] = 'reference'
	df_antim.columns.values[0] = 'Compound Name'
	df_antim.columns.values[1] = 'SMILES'

	df_in = pd.read_csv(infile,header=0,sep=sep)

	# Include additional columns provided by user's input file (other than SMILES, and names if provided)
	if(len(df_in.columns.tolist())>=type):
		df_in = df_in.add_suffix('_user')
	
	if (type ==2):
		df_in.rename(columns={df_in.columns[0]:'Compound Name',df_in.columns[1]:'SMILES'}, inplace = True)
	else:
		df_in.rename(columns={df_in.columns[0]:'SMILES'}, inplace = True)
		df_in['Compound Name'] = 'unknown_cmp '+ (df_in.index.astype('int')+1).astype(str)
	df_in['Type'] = 'input'

	data_df = pd.concat([df_antim,df_in],sort=False,ignore_index=True)
	names = data_df['Compound Name'].tolist()

	# Transform SMILES for clustering
	can_smi = mol.getCanonical(data_df['SMILES'].tolist()) # Obtain canonical SMILES
	msset = mol.getMS(can_smi) # Obtain MOL format
	cmp_mw = mol.getMW(msset) # Obtain molecular weight
	fpset = mol.getFingerprint(msset) # Obtain fingerprints
	#fpset = mol.getMorganFingerprint(msset) # If desired to use morgan fingerprints instead, uncomment this line

	cmp_cl = ClusterFps(fpset, treshold_simil) # Obtain clusters
	#data_df['Canonical SMILES'] = data_df['SMILES'] # To return original SMILES, uncomment this section
	data_df['Canonical SMILES'] = can_smi # To only return provided SMILES, comment this section
	data_df['Cluster'] = cmp_cl
	data_df.to_excel(path_output+'/'+prefix+'_predUntrimmed.xlsx',sheet_name='Prediction',header=True,index=False) # To omit displaying all clusters, comment this line

	# Trim dataset to show only clusters with input query compounds
	list_clus = data_df[data_df['Type'] == 'input'].Cluster.unique()
	data_df = data_df.loc[data_df['Cluster'].isin(list_clus)]
	data_df = data_df.sort_values(by=['Cluster','Type'],ascending=True)

	# Add similarity percentages
	first_clust = data_df.groupby('Cluster').first()['Compound Name'].tolist()
	ind_first_clust = [names.index(i) for i in first_clust]
	
	count = 0
	perc_list = []
	for clust, group in data_df.groupby('Cluster'):
		for i, row in group.iterrows():
			perc = similPercentage(fpset[ind_first_clust[count]],fpset[i])
			perc_list.append(perc)
		count += 1
	data_df['Percentage'] = perc_list	


	# Generate report
	data_df.to_excel(path_output+'/'+prefix+'_pred.xlsx',sheet_name='Prediction',header=True,index=False)

	msg += 'Success: Complete in silico prediction\n'
	return (True,msg)




