#Karla Godinez

# DECLARATIONS
from config import *
import parseInputFunc as pi

# Import libraries
from rdkit import Chem
from rdkit.Chem import Descriptors


# Find duplicates
def checkDuplic(names,smiset):
	dup_ms = []
	try:
		can = [Chem.CanonSmiles(smi) for smi in smiset]
	except:
		return(-1)
	for i in range(len(can)):
		if can[i] in can[i+1:]:
			if names != '':
				index = [names[j] for j, x in enumerate(can) if x == can[i]]
			else:
				index = [can[j] for j, x in enumerate(can) if x == can[i]]
			if index not in dup_ms:
				dup_ms.append(index)
		
	return (dup_ms)



# Validate correct syntaxis and chemistry
def checkSMILES(names,smiset):
	inv_ms = []
	try:
		ms = [Chem.MolFromSmiles(smi,sanitize=False) for smi in smiset]
	except:
		return(-1)

	for i in range(len(ms)):
		if ms[i] is None:
			if names != '':
				inv_ms.append(names[i]) #syntaxis
			else:
				inv_ms.append(smiset[i])
		else:
			try:
				Chem.SanitizeMol(ms[i])
			except:
				if names != '':
					inv_ms.append(names[i]) #chemistry
				else:
					inv_ms.append(smiset[i])
	return (inv_ms)

# Check single SMILES
def getStatusSMILES(smi):
	stat = True
	try:
		m = Chem.MolFromSmiles(smi)
	except:
		stat = False
	if m is None:
		stat = False
	return(stat)


def validateCheck(infile,sep):
	msg = 'Performing Quality Check\n'
	type = 1 # 1 -> SMILES, 2 -> Names + SMILES
	try:	
		f = open(infile,'r')
		header = f.readline()
		reads = f.readlines()
		f.close()
		if len(header.split(sep)) == 1:
			type = 1
			#msg += 'Incorrect file separator\n'
			#return (msg)
		else:
			type = 2
	except:
		msg += 'Error: Reading '+ infile + '\n'
		return (0,msg)

	# Input file check
	header = pi.getHeader(header,sep)
	msg += 'Total of ' + str(len(header)) + ' columns\n'
	msg += 'Columns: '+ ',\t'.join(header)+'\n'

	# Input format
	replace_read = []
	
	for read in reads:
		if len(read.split(sep)) != len(header) or '' in read.split(sep): replace_read.append(read)
	
	if replace_read:
		msg += 'Error: Format error, please fill out empty cells: \n' + '\n'.join(replace_read) + '\n'
		return (0,msg)
	

	# Parse check

	if type == 2:
		cmp_name = pi.getCmpName(reads,sep)
	else:
		cmp_name = ''

	if type == 2:
		cmp_smi = pi.getCmpSmi(reads,sep)
	else:
		cmp_smi = [x.strip() for x in reads]
	#print(cmp_smi)

	msg += 'Total of ' + str(len(cmp_smi)) + ' molecules\n'
	
	# SMILES check
	tmp = checkSMILES(cmp_name,cmp_smi)
	if tmp != -1:
		msg += 'Total of ' + str(len(tmp)) + ' invalid SMILES\n'
	else:
		msg += 'Unable to map SMILES\n'

	msg += ',\t'.join(tmp)
	if (len(tmp)>0): msg += '\n'
	
	# Duplicate check
	tmp = checkDuplic(cmp_name,cmp_smi)
	if tmp != -1:
		tmp = ',\t'.join((', '.join(x) for x in tmp))
		msg += 'Total of '
		msg += '0' if (tmp.split(',')[0]=='') else  str(len(tmp.split(',')))
		msg += ' duplicates\n'
		msg += tmp + '\n'
	
		if (len(tmp)>0): msg += '\n'

	#else:
	#msg += 'Unable to access file\n'
	
	return (type,msg)


