# Input file Parser functions
# Karla Godinez


def getHeader(infile,sep):
	#print ('head\t',infile.strip().split(sep))
	return infile.splitlines()[0].strip().split(sep)

def getCmpName(infile,sep):
	cmp_name = []
	for line in infile:
		tmp = line.strip().split(sep)
		cmp_name.append(tmp[0])
	return cmp_name


def getCmpSmi(infile,sep):
	cmp_smi = []
	for line in infile:
		tmp = line.strip().split(sep)
		cmp_smi.append(tmp[1])
	return cmp_smi

def getCmpExtra(infile,sep,head):
	cmp_extra = []
	for line in infile:
		tmp = line.strip().split(sep)
		if (len(tmp) < len(head)): tmp += ['']*(len(head)-len(tmp))
		if (len(tmp)>2):
			cmp_extra.append(tmp[2:])
		else:
			cmp_extra.append('')

	return cmp_extra
	

def annotationRef2Dic(infile):
	dic = {}
	f = open()
	
	f_out = open('antimalRefDic','w')
	for key,value in dic.items():
		f_out.write('')
	f_out.close()
	



