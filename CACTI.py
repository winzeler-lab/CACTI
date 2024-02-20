## Command line compound search, literature
### Karla Godinez


# DECLARATIONS
import os
import argparse
import time
import sys

from config import *
import getLiterature as lit
import getSimilarCompounds as net
import getPrediction as pred
import qualityCheck as qual


# SET DEFAULTS
sep = '\t'
similthreshold = 0.8 # similarity threshold
analysis = ['lit','simnet','simnetlit','pred','antimval']
keywords = []


# REPLACE VARS BY ARGPARSE
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Compound Search Pipeline - MalDA/UCSD - Winzeler lab')
parser.add_argument('-in', metavar="FILENAME", dest="infile", required=True, help='Input txt file with compound(s) of interest')
parser.add_argument('-predin', metavar="FILENAME", dest="predinfile", help="Prediction with custom input file")
parser.add_argument('-valin', metavar="FILENAME", dest="valinfile", help="Scaffold similarity with custom input file")
parser.add_argument('-outdir', metavar="STRING", dest='outdir', help='PATH to results file',default='./')
parser.add_argument('-jobId', metavar="STRING", dest='jobId', help='Predefined job ID')
parser.add_argument('-p', metavar="STRING", dest='prefix', help='Predefined file name')
parser.add_argument('-k', metavar="STRING", dest='keywords', help='Literature key word search')

#parameters -> type of analysis, needed files
parser.add_argument('-similthreshold', metavar="INT",type=float, dest="similthreshold", help='Threshold for compound similarity. Number between 0,1')
parser.add_argument('-sep', metavar="STRING", dest="sep", help='Indicate input file separator if other than "\\t"')
parser.add_argument('-apincbi', metavar="STRING", dest="apincbi", help='Provide NCBI api if available')


# analysis to perform / figures to obtain
parser.add_argument('-all', action='store_true', dest='all', help='Perform every analysis')
parser.add_argument('-lit', action='store_true', dest='lit', help='Obtain literature search for each compound')
parser.add_argument('-simnet', action='store_true', dest='simnet', help='Obtain similarity network with compound list')
parser.add_argument('-simnetf', action='store_true', dest='simnetf', help='Obtain similarity network figures')
parser.add_argument('-pred', action='store_true', dest='pred', help='Obtain target prediction') ##WRITE UP
parser.add_argument('-antimval', action='store_true', dest='antimval', help='Similarity percentage to antimalarials chemically validated')

# configuration for analysis
parser.add_argument('-sdb', action='store_true', dest='simnetdb', help='Similarity network with dataset and database search')
parser.add_argument('-slit', action='store_true', dest='simnetlit', help='Obtain literature of similar compounds from database lookup')

# set variables
results = parser.parse_args()

if (results.infile):
	infile = results.infile
if (results.predinfile):
	predinfile = results.predinfile
else:
	predinfile = ''
if (results.valinfile):
	valinfile = results.valinfile
else:
	valinfile = ''
if (results.similthreshold):
	similthreshold = results.similthreshold
if (results.sep):
	sep = results.sep
if (results.jobId):
	job_ID = results.jobId
else:
	job_ID = str(time.time())
if (results.prefix and results.prefix != ''):
	prefix = results.prefix
else:
	prefix = job_ID
#output_path = results.outdir+prefix
output_path = results.outdir
if (results.apincbi):
        api_ncbi = results.apincbi
else:
        api_ncbi = ''
if (results.keywords):
	keywords = results.keywords.split(' ')


# CHECK FLAGS FOR ANALYSIS
tmp = []
if (not results.all):
	if (results.lit):
		tmp.append('lit')
	if (results.simnet):
		tmp.append('simnet')
	if (results.pred):
		tmp.append('pred')
	if (results.antimval):
		tmp.append('antimval')

if (len(tmp)!=0):
	analysis = tmp

# CHECK FLAGS FOR ANALYSIS CONFIGURATION
sim_conf = []
if (results.simnetdb):
	sim_conf.append('simnetdb')
if (results.simnetlit):
	sim_conf.append('simnetlit')
if (results.simnetf):
	figs = True
else:
	figs = False



## To have UNIQUE output files, uncomment section below
'''
# OUTPUT FOLDER
if (os.path.exists(output_path)):
	print('Folder not valid\nPlease remove prior folder or rename and resubmit')
	log = open(output_path+'/'+job_ID+'.log', 'w')
	log.write('Invalid output folder, remove prior folder or rename')
	log.close()	
	sys.exit()
'''

#os.mkdir(output_path)
log = open(output_path+'/'+job_ID+'.log', 'w')
log.write('Output folder '+output_path+'\n')

# ERROR LOG
log.write('Processing of submission has started.\n')
#log.close()

type,msg = qual.validateCheck(infile,sep)
log.write(msg)

if type == 0: sys.exit()
	
# CALL FUNCTIONS
for i in analysis:
	if (i == 'lit'):
		print('Processing Literature Search')
		stat,msg=lit.getLit(type,prefix,infile,output_path,sep,api_ncbi,keywords)
		log.write(msg)
	if (i == 'simnet'):
		print('Processing Similarity Network')
		stat,msg=net.getNetwork(type,prefix,infile,output_path,sep,similthreshold,sim_conf,figs)
		log.write(msg)
	if (i == 'pred'):
		print('Processing In Silico Target Prediction')
		stat,msg=pred.getPrediction(type,prefix,infile,output_path,sep,similthreshold,predinfile)
		log.write(msg)
	if (i == 'antimval'):
		print('Processing Similarity Percentages to Validated Antimalarials')
		stat,msg=pred.getAntimSimil(type,prefix,infile,output_path,sep,valinfile)
		log.write(msg)


# CHECK FOR END RUN
print('Complete')
log.write('\nComplete')
log.close()



