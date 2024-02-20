#Karla Godinez

# Import libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors


# Convert SMILES to fingerprint
def getMS(smiles):
	if type(smiles) != list:
		return Chem.MolFromSmiles(smiles)

	# Check list size, get molecule accordingly
	ms = [Chem.MolFromSmiles(x) for x in smiles]
	
	'''
	if len(list_smiles) <= 500:
		ms = [Chem.MolFromSmiles(x) for x in list_smiles]
	
	else:
		ms = [Chem.MolFromSmiles(x) for x in list_smiles if MolWT(x) <= 300]
		ms = [Chem.MolFromSmiles(x) for x in list_smiles if MolWT(x) > 300]
	'''
	
	return ms


# Get Molecular weight
def getMW(ms):
	mw = [round(Descriptors.MolWt(x),3) for x in ms]

	return mw

def getMWfromSMI(smi):
	mw = Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))
	return mw

# Get Tanimoto similarity
def getTanSim(fp1,fp2):
	sim = float(DataStructs.TanimotoSimilarity(fp1, fp2))
	return sim


# Convert SMILES to fingerprint
# RDKIT Fingerprint
def getFingerprint(ms):
	if type(ms) != list:
		return Chem.RDKFingerprint(ms)

	fp = [Chem.RDKFingerprint(x) for x in ms]
	
	return fp

def getMorganFingerprint(ms,rad=4,nbit=2048):
	if type(ms) != list:
		return AllChem.GetMorganFingerprintAsBitVect(ms, useChirality=True, radius=rad, nBits=nbit)

	fp = [AllChem.GetMorganFingerprintAsBitVect(x, useChirality=True, radius=rad, nBits=nbit) for x in ms]

	return fp
	
def getFingerprintParts(n,idset,mwset,msset):
	fp = []
	ind_fp = []
	fp_dic = {501:[]}
	
	for i in range(1,n+1):
		fp_dic[(500/n)*i] = []
	
	buckets = list(fp_dic.keys())
	buckets.sort()

	# Determine bucket for each compound	
	for i in range(len(msset)):
		buck = -1
		for j in range(len(buckets)-1):
			if (mwset[i] <= buckets[j] and buck == -1):
				buck=j

		fp_dic[buckets[buck]].append(msset[i])
		ind_fp.append([i,buckets[buck]])

	# Get fingerprints per bucket
	fp_buck = {}
	for i in range(len(buckets)):
		print(buckets[i],len(fp_dic[buckets[i]]))
		fp_buck[buckets[i]] = [Chem.RDKFingerprint(x) for x in fp_dic[buckets[i]]] #possibly get len and do 2048 bits or else
	
	# Return fingerprints to msset position
	ind_buckets = [0]*len(buckets)
	print(len(fp))
	for element in ind_fp:
		i_ms = element[0]
		i_b = buckets.index(element[1])
		
		fp.append(fp_buck[element[1]][i_b])
		ind_buckets[i_b] += 1

	return (fp)
	


# Morgan Fingerprint
def getMorganFingerprint4(ms):

	from rdkit.Chem import rdMolDescriptors
	fp = [rdMolDescriptors.GetMorganFingerprint(x,4,nBits=2048,useFeatures=True) for x in ms]
	
	return fp
	
# ECFP Fingerprint with Features and bit
def getMorganBit(ms):
	from rdkit.Chem import AllChem
	fp = [AllChem.GetMorganFingerprintAsBitVect(x,4,nBits=2048,useFeatures=True) for x in ms]

	return fp

# ECFP Morgan Fingerprint with Features no bit
def getMorgan(ms):
	from rdkit.Chem import AllChem
	fp = [AllChem.GetMorganFingerprint(x,4,useFeatures=True) for x in ms] 

# Get canonical SMILES
def getCanonical(smiset):
	#ms = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiset]
	#print(smiset)
	ms = []
	for smi in smiset:
		try:
			mol = Chem.MolFromSmiles(smi)
			molsmi = Chem.MolToSmiles(mol,True)
			ms.append(molsmi)
		except:
			#print(smi)
			x=1

	return ms


# Get similarity percentage SMILES
def similPercentSmiles(smi1,smi2):
	perc = 0
	mol1 = Chem.MolFromSmiles(smi1)
	mol2 = Chem.MolFromSmiles(smi2)

	if (mol1 != None and mol2 != None):
		fp1 = Chem.RDKFingerprint(mol1)
		fp2 = Chem.RDKFingerprint(mol2)

		perc = DataStructs.FingerprintSimilarity(fp1,fp2)*100

	return(perc)


	
def smi2svg(smi):
	from rdkit.Chem.Draw import rdMolDraw2D
	from rdkit.Chem import AllChem

	mol = Chem.MolFromSmiles(smi)
	try:
		Chem.rdmolops.Kekulize(mol)
	except:
		pass
	drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
	AllChem.Compute2DCoords(mol)
	drawer.DrawMolecule(mol)
	drawer.FinishDrawing()
	svg = drawer.GetDrawingText().replace("svg:", "")
	return svg

	
# Generate image from SMILES
def smi2image(smi):
	from urllib import parse
	svg_string = smi2svg(smi)
	impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg_string, safe="")
	return impath


# Generate SDF from Mol
def mol2sdf(out_path,molse):
	w = Chem.SDWriter(out_path,molset)
	for m in molset:
		w.write(m)


