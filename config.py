import sys

REF_DIR = 'ref/' #Reference folder
PYTHON_DIR = 'python/' #Functionality scripts folder
BASEFUNC_DIR = 'baseFunc/' #Functionality scripts folder
#Add python paths
sys.path.insert(0,'%s'%PYTHON_DIR)
sys.path.insert(0,'%s'%REF_DIR)
sys.path.insert(0,'%s'%BASEFUNC_DIR)
