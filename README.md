# Toxidrome - Graph Convolutional Neural Network Models

##### Supporting information for paper:
Title: "Rapid Screening of Chemicals for their Potential to Cause Specific Toxidromes"
Authors: Ruifeng Liu, Mohamed Diwan M AbdulHameed, Zhen Xu, Benjamin Clancy, Valmik Desai, Anders Wallqvist
Journal: Frontiers in Drug Discovery, section In silico Methods and Artificial Intelligence for Drug Discovery
### Intro
This repository contains the data and models used to make Toxidrome's predictions and has sorted this data into 4 sections:
- An excel workbook containing the compounds used in training our graph convolutional neural network (GCNN) models as well as running our similarity ensemble approach (SEA)
- A folder containing the CMPNN models and script
- A folder containing the DMPNN models and script
- A folder containing the SEA script

### Compound Data Excel Workbook
This workbook contains 8 pages, one for each toxidrome we make predictions for (with one exception: cholinergic and anticholinergic predictions are made using the same compounds).
The pages used in training our GCNN models (cholinergic, convulsant, opioid, and sympathomimetic) contain the each compound's SMILES, model, value, and unit, while the pages used in running the SEA (Anticoagulant, Irritant-corrosive, knockdown, and Solvents-anesthetics-sedatives) only contain the SMILES for each compound.

### CMPNN Models and Script
We used Communicative message passing neural network (CMPNN) approach developed by _Song et al. Proceedings of the Twenty-Ninth International Joint Conference on Artificial Intelligence (IJCAI-20): 2020 (https://github.com/SY575/CMPNN)._

This folder contains 2 main parts: a folder of models called "toxidromes" and a script called run_cmpnn.py.

This script runs the CMPNN using these models.
### DMPNN Models and Script
We used Directed message passing neural network (DMPNN) approach developed by _Yang et al. Analyzing learned molecular representations for property prediction. J Chem Inf Model 2019, 59(8):3370-3388 (https://chemprop.readthedocs.io/)_.

This folder contains 2 main parts: a folder of models called "toxidromes" and a script called run_dmpnn.py.

This script runs the DMPNN using these models.
### SEA Script
This folder contains a script to run the similarity ensemble approach and an excel sheet which has the list of compounds used in running the SEA script (these are the same compounds listed in the SEA pages the workbook). 

### Running the python script
To run the python script, you will first need Anaconda installed. From an Anaconda prompt, sut up a new environment using the following commands:

`conda create -n toxidrome python=3.7`
`conda activate toxidrome`

Next, navigate to the "command_line_tool" folder and enter the following command to install Toxidrome's dependencies:

`pip install -r requirements.txt`

Once everything is installed, you can then run the script by running `python toxidrome_clt.py` followed by any of the following tags (including at least 1 that adds compounds to the job):
- `-h` or `--help`: Shows a help message explaining these tags.
- `-i [INPUT]` or `--input [INPUT]`: The file location of a .CSV file whose first column is 'Name' and whose second is 'SMILES' and contains the list of SMILES to be submitted. This tag counts for adding compounds to the job.
- `-ln [NAMES]` or `--names [NAMES]`: A delimited list of compound names.
- `-ls [SMILES]` or `--smiles [SMILES]`: A delimited list of compound SMILES. This tag counts for adding compounds to the job.
- `-o [OUTPUT]` or `--output [OUTPUT]`: The output file path
- `-f [FORMAT]` or `--format [FORMAT]`: File format for the output. The options are 'json', 'csv', or 'database'. The json format produces a json file and the csv and database formats produce a csv file. Note: the 'csv' format does not output all information, only the main information in a readable format.
- `-e [ERROR]` or `--error [ERROR]`: The file path to output error messages.
- `-sf [STATE]` or `--state [STATE]`: The file path for outputting state info.
- `-m [MODE]` or `--mode [MODE]`: The mode of execution: 'performance' or 'robust'. default: performance
- `-b [BULK]` or `--bulk [BULK]`: The amount of compounds that go through Main Predictions each time in order to avoid a memory issue.
- `-nnType [NNTYPE]`: Type of neural network to run. Either 'cmpnn' or 'dmpnn'. default: dmpnn

Here is an example prompt:

`python toxidrome_clt.py -i /saved_extra_files/input.csv -o /saved_extra_files/output.csv -f database -nnType cmpnn`

This prompt takes the input csv from "/saved_extra_files" and outputs the output csv from the same directory.