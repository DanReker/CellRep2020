######
# imports
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import csv, math
import pylab as pl
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import sys, re, os

######
# setting global variables
directory="chembl22"
mf=None
n_trees=500;
descr = Descriptors._descList
calc = [x[1] for x in descr]

######
# define function to calc descriptors
def descr_calc(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=2048)  # MorganFP
    ds_n = []
    for d in calc:
        v = d(mol)
        if math.isnan(v):  # postprocess descriptors for NaN or freak large values
            ds_n.append(np.float32(0.0))
        elif v > np.finfo(np.float32).max: 
            ds_n.append(np.finfo(np.float32).max)
        else:
            ds_n.append(np.float32(v)) 
    
    fp_list = []
    fp_list.extend(fp.ToBitString())
    fp_expl = [float(x) for x in fp_list]
    
    return fp_expl + list(ds_n)


######
# read activity data for a specific target
#tid = sys.argv[1]
tid = "11398"
f = open(directory + "/" + tid + "_data.tsv", 'r');

fingerprints = [];
activities = [];
chembl_ids = [];
reader = csv.reader(f,delimiter='\t');

next(reader, None) # skip the header
for entry in reader:
    mol = Chem.MolFromSmiles(entry[2])
    if not mol is None:
        fingerprints += [descr_calc(mol)];
        activities += [float(entry[3])];
        chembl_ids += [entry[0]];

fingerprints = np.array(fingerprints,dtype=np.float32) # convert from list to numpy array
activities = np.array(activities)



######
# read structures for GRAS / IIG compounds
GRASIIG_smiles_file = open("gras_iig.tsv",'r')
GRASIIG_fingerprints = []
GRASIIG_names = []

reader2 = csv.reader(GRASIIG_smiles_file,delimiter='\t');
for entry in reader2:
    mol = Chem.MolFromSmiles(entry[2])
    if not mol is None:
        GRASIIG_fingerprints += [descr_calc(mol)];
        GRASIIG_names += [entry[0]]

GRAS_fingerprints = np.array(GRASIIG_fingerprints,dtype=np.float32)

######
# read structures for background compounds
BG_smiles_file = open("background.tsv",'r')
BG_fingerprints = []

reader3 = csv.reader(BG_smiles_file,delimiter=' ');
for entry in reader3:
    mol = Chem.MolFromSmiles(entry[0])
    if not mol is None:
        BG_fingerprints += [descr_calc(mol)]

BG_fingerprints = np.array(BG_fingerprints,dtype=np.float32)

######
# train model
RF = RandomForestRegressor(n_estimators=n_trees,n_jobs=1,max_features=mf);
RF.fit(fingerprints,activities)

######
# predict GRAS/IIG and background compounds and normalize 
predictions = RF.predict(GRASIIG_fingerprints)
predictions_bg = RF.predict(BG_fingerprints)
z_scores = (predictions - np.mean(predictions_bg)) / np.std(predictions_bg) 

######
# print top predictions
for (name, score) in zip(np.array(GRASIIG_names)[np.argsort(-z_scores)], -np.sort(-z_scores)):
	if score > 1.5:
		print name, score


