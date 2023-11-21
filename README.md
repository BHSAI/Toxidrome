# Toxidrome - Graph Convolutional Neural Network Models

##### Supporting information for paper:
Title: "Rapid Screening of Chemicals for their Potential to Cause Specific Toxidromes"
Authors: Ruifeng Liu, Mohamed Diwan M AbdulHameed, Zhen Xu, Benjamin Clancy, Valmik Desai, Anders Wallqvist
Journal: Frontiers in Drug Discovery, section In silico Methods and Artificial Intelligence for Drug Discovery
### Intro
This repository contains the data and models used to make Toxidrome's predictions and has sorted this data into 4 sections:
- An excel workbook containing the compounds used in training our graph convolutional neural network (GCNN) models as well as running our similarity ensemble approach (SEA)
- A folder containing the CMPNN models and script
- A folder containing the DMPNN mdoels and script
- A folder containing the SEA script

### Compound Data Excel Workbook
This workbook contains 8 pages, one for each toxidrome we make predictions for (with one exception: cholinergic and anticholinergic predictions are made using the same compounds).
The pages used in training our GCNN models (cholinergic, convolsant, opioid, and sympathomimetic) contain the each compound's SMILES, model, value, and unit, while the pages used in running the SEA (Anticoagulant, Irritant-corrosive, knockdown, and Solvents-anesthetics-sedatives) only contain the SMILES for each compound.

### CMPNN Models and Script
[what is different about the CMPNN]

This folder contains 2 main parts: a folder of models called "toxidromes" and a script called run_cmpnn.py.

This script runs the CMPNN using these models.
### DMPNN Models and Script
[what is different about the DMPNN]

This folder contains 2 main parts: a folder of models called "toxidromes" and a script called run_dmpnn.py.

This script runs the DMPNN using these models.
### SEA Script
This folder contians a script to run the similarity ensemble approach and an excel sheet which has the list of compounds used in running the SEA script (these are the same compounds listed in the SEA pages the workbook). 