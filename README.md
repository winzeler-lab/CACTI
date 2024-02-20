<img src="https://drive.google.com/uc?id=12A2gMVLFWFhG94D5kn5YSoMueQ8Mnj_8">
Open-source annotation and target prediction tool that explores some of the largest chemical and biological databases, mining these for identification of common name, synonyms, and structurally similar molecules.

## Getting Started

Scripts can be executed using Linux-based command line. Please note that reference folders must be present for CACTI correct execution.

### Prerequisites

Python			3.0
xlrd 			1.2
matplotlib		2.0
scipy			1.4
Rdkit			2018.03.4
networkx			2.3
numpy			1.15
pandas			0.23
requests			2.12
scikit-learn		0.20



### Conda environment
This project can be executed with preinstalled packages, or using a conda environment. To install the Anaconda environment, please type
```
conda env create --file metadata/CACTI-env.txt --name CACTI python=3
```
Activate the directory by typing
```
conda activate CACTI
```
## Execute
1. Change directory into that containing the base script (CACTI.py). Please note that if using anaconda environment, it must be activated first
2. Execute command line script

### Full analysis
To execute full analysis (synonyms, literature, similar analogs and clustering analysis), please type
```
python CACTI.py -in filename.txt
```
To save the output to a particular directory, one can type the path using:
```
python CACTI.py -in filename.txt -out /outputdir
```
To add an optional output name:
```
python CACTI.py -in filename.txt -out /outputdir -p outname
```

Description of available options:
```
python CACTI.py -h
```

### Partial analysis
#### Module 1: Synonyms and literature

```
python CACTI.py -in filename.txt -lit
```
To obtain results faster and if available, one could use their own NBCI API code when executing the literature module
```
python CACTI.py -in filename.txt -lit -apincbi apiValue
```
If desired to find literature evidence with a keyword, please type the following
```
python CACTI.py -in filename.txt -lit -k keyword1,keyword2
```
#### Module 2: Close analogs
To  find close analogs and perform similarity analysis between them, with 80% Tanimoto similarity
```
python CACTI.py -in filename.txt -simnet -sdb
```
To change the similarity threshold, add the corresponding flag, where X is a numerical value between 0-1
```
python CACTI.py -in filename.txt -simnet -sdb -similthreshold X
```

#### Module 3: Clustering analysis
The clustering analysis is performed using the drug-target pair collection, provided in the metadata folder.
```
python CACTI.py -in filename.txt -pred -similthreshold X
```
If there is an particular clustering reference file to use instead, one can provide the path by typing
```
python CACTI.py -in filename.txt -pred -predin alternative_textfile.txt -similthreshold X
```










