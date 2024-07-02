# JSON Parser functions
import json
import urllib.parse as urllib
import requests
import re
#from urllib.error import HTTPError


# VARIABLES
#ommitSyn = ['MCULE-','SMR','MLS','AKOS','SR-','HMS','EU','OPERA','OPREA','MAYBRIDGE','ZINC','IDI']
ommitSyn = ['SMR','MLS','SR-','HMS','EU','OPERA','OPREA','IDI']

def getPubchemCID(cmp,smiles):
	query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smiles.replace('/','.'))+'/cids/JSON?MaxRecords=20'
	result = requests.get(query).json()
	if('Fault' in result and cmp!=''):
		query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.quote(cmp)+'/cids/JSON?MaxRecords=20'
		result = requests.get(query).json()
	if (bool(result) and 'IdentifierList' in result):
		ids = filter(lambda a: a != 0, result['IdentifierList']['CID'])
		return (';'.join(str(x).replace("'",'') for x in ids))
	else:
		return (' ')


def getPubmedID(cmp,db,api,keywords):
	import re
	full_list = []
	for item in cmp.split('\n'):
		try:
			val = [s for s in ommitSyn if item[:re.search(r"\d", item).start()].upper() in s]
		except:
			val = []
			
		if (len(val) == 0):
			query = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?&retmode=json&db='+db+'&term='+urllib.quote(item)+'+('+' or '.join(keywords)+')'+'&retmax=100'
			if (bool(api)): query += '&api='+api
			result = requests.get(query).json()
			if ('esearchresult' in result.keys()):
				if ('ERROR' not in result['esearchresult'] and int(result['esearchresult']['count'])>0):
					for id in result['esearchresult']['idlist']:
						if (id not in s for s in full_list):
							if (db == 'pubmed'):
								full_list.append(item+'\tPMID'+id+'\t https://pubmed.ncbi.nlm.nih.gov/'+id)
							elif (db == 'pmc'):
								full_list.append(item+'\tPMC'+id+'\t https://www.ncbi.nlm.nih.gov/pmc/articles/PMC'+id)

	if not full_list: return(' ')
	return ('\n'.join(full_list))

def getDOI(db,id):
	import re

	doi = ''
	query = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db='+db+'&id='+str(id)+'&rettype=medline'
	result = requests.get(query).text
	
	ind = re.search('doi:',result).start()
	doi = result[ind : re.search('. ',result[ind:]).start()]

	return 	(doi)


def getPubmedID_CID(id):
	full_list = []
	ids = []
	for item in id.split(';'):
		query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(item)+'/xrefs/PubMedID/JSON'
		result = requests.get(query).json()
		if (bool(result) and 'InformationList' in result):
			for cid in result['InformationList']['Information']:
				for i in cid['PubMedID']:
					if i not in ids:
						full_list.append('CID:'+item+'\t'+'PMID'+str(i)+'\t https://pubmed.ncbi.nlm.nih.gov/'+str(i))
						ids.append(i)

	return ('\n'.join(full_list))

def getPatent_CID(id):
	full_list = []
	for item in id.split(';'):
		query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(item)+'/xrefs/PatentID/JSON'
		result = requests.get(query).json()
		if (bool(result) and 'InformationList' in result):
			for cid in result['InformationList']['Information']:
				for pat in cid['PatentID']:
					pat = pat.strip()
					if str(pat) not in full_list : full_list.append(str(pat))
	if not full_list: return(' ')
	if (len(full_list) > 100) : full_list = full_list[:100]
	return ('\n'.join(full_list))

def getEBIID(id):
	import re
	full_list = []
	tmp = []
	for item in id.split('\n'):
		try:
			val = [s for s in ommitSyn if item[:re.search(r"\d", item).start()].upper() in s]
		except:
			val = []

		if (len(val) == 0):
			query = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&query='+'"'+str(item)+'"'
			result = requests.get(query).json()
			if (bool(result) and ('errCode' not in result) and len(result.keys())>1):
				if(result['hitCount']>0):
					for i in result['resultList']['result']:
						if (str(i['id']) not in s for s in full_list): full_list.append((item+'\tPMID'+str(i['id'])+'\t https://pubmed.ncbi.nlm.nih.gov/'+str(i['id']) if str(i['id']).isdigit() else item+'\t'+str(i['id'])))

	return ('\n'.join(full_list))


def getChemblMech(id):
	full_list = []
	ids = []
	for item in id.split('\n'):
		if (item != ''):
			query = 'https://www.ebi.ac.uk/chembl/api/data/mechanism.json?limit=100&molecule_chembl_id='+item
			result = requests.get(query).json()
			if (bool(result) and 'mechanisms' in result):
				for mec in result['mechanisms']:
					for mec_id in mec['mechanism_refs']:
						if mec_id['ref_id']+mec_id['ref_type'] not in ids:
							full_list.append('MOA: '+mec['mechanism_of_action']+' REF: '+mec_id['ref_url'])
							ids.append(mec_id['ref_id']+mec_id['ref_type'])

	return (str('\n'.join(full_list)))


def getChemblDoc(id):
	query = 'https://www.ebi.ac.uk/chembl/api/data/document/' + id +'.json'
	result = requests.get(query).json()
	if (bool(result)):
		line = result['doc_type'] + ' ID: '+result['document_chembl_id']
		if result['title'] : line += ' Title: '+result['title']
		if result['doi'] : line += ' DOI: '+result['doi'] + ' https://www.doi.org/'+result['doi']

		return(line)



def getChemblLit(id):
	full_list = []
	ids = []
	for item in id.split('\n'):
		if (item != ''):
			query = 'https://www.ebi.ac.uk/chembl/api/data/compound_record.json?limit=100&molecule_chembl_id='+item
			result = requests.get(query).json()

			if (bool(result) and 'compound_records' in result):
				for id in result['compound_records']:
					if (id['document_chembl_id'] not in ids):
						line = getChemblDoc(id['document_chembl_id'])
						full_list.append(line)
						ids.append(id['document_chembl_id'])

	return('\n'.join(full_list))


def getBindingDB(smi, simil):
	query = 'http://www.bindingdb.org/axis2/services/BDBService/getTargetByCompound?smiles=' + urllib.quote(str(smi)) \
			+ '&cutoff=' + str(simil) + '&response=application/json'
	result = requests.get(query).json()
	full_list = []
	if (bool(result) and 'bdb.getTargetByCompoundResponse' in result):  # If correct smiles
		if (result['bdb.getTargetByCompoundResponse']['bdb.hit'] > 0):  # If hits found
			full_list = []
			ids = []

			if type(result['bdb.getTargetByCompoundResponse']['bdb.affinities']) != list:
				result['bdb.getTargetByCompoundResponse']['bdb.affinities'] = [
					result['bdb.getTargetByCompoundResponse']['bdb.affinities']]

			for item in result['bdb.getTargetByCompoundResponse']['bdb.affinities']:
				if item['bdb.inhibitor'] + item['bdb.target'] not in ids:
					key = item['bdb.inhibitor'] + '\t' + item['bdb.target']
					ids.append(key)
					full_list.append(item['bdb.inhibitor'] + '\t' + item['bdb.species'] + '\t' + 'Target:' + item[
						'bdb.target'] + '\t' + 'Tanimoto similarity:' + item['bdb.tanimoto'])
	if full_list:
		return ('\n'.join(full_list))
	else:
		return ('')

def getBindingDBSimil(smi,simil):
	query = 'http://www.bindingdb.org/axis2/services/BDBService/getTargetByCompound?smiles=' + urllib.quote(str(smi)) \
			+ '&cutoff=' + str(simil) + '&response=application/json'

	result = requests.get(query).json()
	full_list = []
	ids = []
	if (bool(result) and 'bdb.getTargetByCompoundResponse' in result):  # If correct smiles
		print(result['bdb.getTargetByCompoundResponse']['bdb.hit'])
		if (result['bdb.getTargetByCompoundResponse']['bdb.hit'] == 1):
			item = result['bdb.getTargetByCompoundResponse']['bdb.affinities']
			if item['bdb.inhibitor'] + item['bdb.smiles'] not in ids:
				key = item['bdb.inhibitor'] + '\t' + item['bdb.smiles']
				ids.append(key)
				full_list.append([item['bdb.inhibitor'], item['bdb.smiles']])
		elif (result['bdb.getTargetByCompoundResponse']['bdb.hit'] > 1):  # If hits found
			for item in result['bdb.getTargetByCompoundResponse']['bdb.affinities']:
				if item['bdb.inhibitor'] + item['bdb.smiles'] not in ids:
					key = item['bdb.inhibitor'] + '\t' + item['bdb.smiles']
					ids.append(key)
					full_list.append([item['bdb.inhibitor'],item['bdb.smiles']])

	return(full_list)




def getChemblSmiID(smi,simil_perc):
	ids = []
	query = 'https://www.ebi.ac.uk/chembl/api/data/similarity/'+urllib.quote(smi)+'/'+str(simil_perc)+'.json'
	result = requests.get(query).json()

	if (bool(result) and 'molecules' in result):
		for mol in result['molecules']:
			if mol['molecule_chembl_id'] not in ids: ids.append(mol['molecule_chembl_id'])

	return (str('\n'.join(ids)))



def getChemblSimilPerc(smi,treshold_simil):
	#Structure [mol ID,canonical smiles, similarity percentage]
	pairs = []

	query= 'https://www.ebi.ac.uk/chembl/api/data/similarity/'+urllib.quote(smi)+'/'+str(treshold_simil)+'.json'
	result = requests.get(query).json()

	if (bool(result) and 'molecules' in result):
		for mol in result['molecules']:

			if (mol['molecule_structures']['canonical_smiles'] != smi and int(float(mol['similarity'])) < 100):
				pairs.append([mol['molecule_chembl_id'],mol['molecule_structures']['canonical_smiles'],mol['similarity']])
	return pairs

def getCIDSmiles(cid):
	query='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+urllib.quote(str(cid))+'/property/IsomericSmiles/JSON'
	result = requests.get(query).json()
	#id = result['PropertyTable']['Properties'][0]['IsomericSMILES'] # To obatin SMILES (potential duplication between CID SMILES, uncomment this line
	id = result['PropertyTable']['Properties'][0]['IsomericSMILES']

	return(id)


def getPubChemSubstructure(smi,treshold_simil):
	#[mol ID,canonical smiles]
	substruc = []
	# Limit to return 100 analogs per search, this number can be increased if desired
	query= 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/'+urllib.quote(smi)+'/cids/JSON?Threshold='+str(treshold_simil)+'&MaxRecords=100'
	result = requests.get(query).json()
	if (bool(result) and 'IdentifierList' in result):
		for p in result['IdentifierList']['CID']:
			substruc.append([p,getCIDSmiles(p)])

	return (substruc)

def getSynonym(cmp,smi):
	syn_list = []
	query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smi.replace('/','.'))+'/synonyms/json'
	result = requests.get(query).json()

	# In case SMILES issue but valid, 
	if('Fault' in result and cmp!=''):
		query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.quote(cmp)+'/synonyms/json'
		result = requests.get(query).json()
	if (bool(result) and 'InformationList' in result):
		for item in result['InformationList']['Information']:
			for syn in item['Synonym']:
				if syn not in syn_list and syn != cmp : syn_list.append(syn)

	return str('\n'.join(syn_list))


def getMW(smi):
	query = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'+urllib.quote(smi)+'/property/MolecularWeight/json'
	#print(query)
	result = requests.get(query).json()
	if (bool(result) and 'PropertyTable' in result):
		result = result['PropertyTable']['Properties'][0]
		return (str(result['MolecularWeight']))
	else:
		return('')





