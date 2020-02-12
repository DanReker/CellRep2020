# CellRep2020

Data and code for Reker et al Cell Rep 2020.

#### /Raw bio data
contains raw experimental data of the *in vitro*, *ex vivo*, *in vivo* validations.

#### /chembl22
contains raw chembl22 data for model training. in this repository we only include P-gp data as a single example to run prediction code

#### prediction.py
Python code to run predictions of GRAS/IIG bioactivities

#### gras_iig.tsv
SMILES structures of GRAS/IIG compounds used in predictions

#### background.tsv
SMILES structures of molecular weight-matched background library of random chemicals

#### test_data.tsv
known IC50/Ki activities of GRAS/IIG compounds from ChEMBL22, serves as validation data for model optimization/selection
