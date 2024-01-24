import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
import math
from chembl_structure_pipeline import standardizer
from rdkit.Chem.MolStandardize import rdMolStandardize
import base64
from io import BytesIO
import os
import time
import json
import pandas as pd
from rdkit.Chem import Descriptors

from toxidrome_lib.constants import *

class Channels:
    errorChannel = sys.stderr
    stateChannel = None
    outputChannel = sys.stdout

    def writeStateMsg(self, msg):
        if self.stateChannel != None:
            print(msg, file= self.stateChannel)
            self.stateChannel.flush()

#- read toxidrome models that are used in FP similarity calculation in SEA
# Create a mapping of Toxidrome to molecules
# at the time of writing the toxidrome column contains three values:
# - Anticoagulant
# - Irritant, Corrosive
# - Solvents, anesthetics, sedatives
#- calculate the FP for the molecules in the previous model
#Assume the predefined smiles can go through the FP calculation without any problem
def loadToxidromeData_forSEA(dataPath):
    toxidromeMols_data = pd.read_excel(dataPath)
    toxidromeMols = toxidromeMols_data.groupby('Toxidrome')['Molecule'].agg(list).to_dict()
    toxidromeFPs = {}
    for toxidrome in toxidromeMols:
        toxidromeFPs[toxidrome] = []
        for smiles in toxidromeMols[toxidrome]:
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2, nBits=2048)
            toxidromeFPs[toxidrome].append(fp)
    return toxidromeMols, toxidromeFPs


largest_Fragment = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()
def standardize_one_smiles(smiles):
    if smiles is not None:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            #get the larget fragment
            largest_mol = largest_Fragment.choose(molecule)
            m_no_salt = standardizer.get_parent_mol(largest_mol)
            molecule_std = standardizer.standardize_mol(m_no_salt[0])
            # try to neutralize molecule
            #uncharged_mol = uncharger.uncharge(molecule_std)
            smiles_std = Chem.MolToSmiles(molecule_std)
            return smiles_std
    return None

def parse_input_csv(csv_path):
    with open(csv_path, 'r') as file:
        csvreader = csv.reader(file)
        names = []
        smiles = []
        for row in csvreader:
            if len(row) == 1:
                smiles.append(row[0])
            else:
                names.append(row[0])
                smiles.append(row[1])
        if len(names) > 0:
            names.pop(0)
        if len(smiles) > 0:
            smiles.pop(0)
        return names, smiles

#TODO what is rs or n
#Calculate what? reference maybe?
def calcP(rs, n):
    z = (rs - (0.000424 * n)) / (0.00449 * (n ** 0.665))
    xz = -(math.e ** ((-1 * z * math.pi) / (math.sqrt(6) - 0.577215665)))
    if z <= 28:
        return (1 - (math.e ** xz))
    else:
        return ((-1 * xz) - ((xz ** 2) / 2) - ((xz ** 3) / 6))

# build and print the data into channel
# the channel is assumed to be already opened by the caller
def write_json(datapassList, channel):
    for datapass in datapassList:
        for toxidrome in MAIN_TOXIDROMES:
            models = TOX_CATEGORY[toxidrome]
            for model in models:
                acc = datapass.getToxidromeAccuracy(toxidrome, model)
                if acc != None:
                    datapass.setToxidromeAD(toxidrome, model, acc < 1)
        for toxidrome in SEA_TOXIDROMES:
            val = datapass.toxidromes[toxidrome]
            datapass.setToxidromeAD(toxidrome, None, val)
    jsonStr = json.dumps(datapassList, default = lambda x: x.__dict__)
    print(jsonStr, file=channel)


# Create a csv file where each row represents one data entry
# which can be directly read into the Java Result entity
# Header: name, toxidrome, subcategory, value, accuracy
# subcategory refers to 'model' in the Java entity
# name refers to the 'compoundid' in the Java entity
def write_database(datapassList, channel):
    writer = csv.writer(channel)
    headers = ['name', 'smiles', 'toxidrome', 'subcategory','value','accuracy','appDomain']
    _ = writer.writerow(headers)
    # two classes of toxidromes
    # toxidromes in MAIN_TOXIDROMES have both value and accuracy
    # the data structure for them:
    # toxidrome => subcategory => <data>
    # where data can be value or accuracy depending on the mapping field
    for datapass in datapassList:
        for toxidrome in MAIN_TOXIDROMES:
            for subcategory in TOX_CATEGORY[toxidrome]:
                value, accuracy = datapass.getToxidromeData(toxidrome, subcategory)
                if (accuracy != None): accuracy = round(accuracy, ROUND_DIGIT)
                appDomain = accuracy != None and accuracy < AD_CUTOFF
                newRow = [datapass.name, datapass.smiles, toxidrome, subcategory, value, accuracy, appDomain]
                _ = writer.writerow(newRow)
        for toxidrome in SEA_TOXIDROMES:
            value = datapass.getToxidromeValue(toxidrome)
            appDomain = value != None and value < SEA_THRESHOLD
            newRow = [datapass.name, datapass.smiles, toxidrome, None, value, None, appDomain]
            _ = writer.writerow(newRow)

#fileToWrite can be stdout
#up to the caller to handle file open and close
#sd controls if toxidrome model data or accurracy shall be printed
def write_csv(datapassList, channel, sd):
    writer = csv.writer(channel)

    models = TOX_CATEGORY
    main_toxidromes = MAIN_TOXIDROMES
    other_toxidromes = SEA_TOXIDROMES

    # document why two titles
    title_row1 = [NAME_TITLE,SMILE_TITLE]
    title_row2 = [NAME_TITLE,SMILE_TITLE]
    for toxidrome in main_toxidromes:
        title_row1 += ([toxidrome] * len(models[toxidrome]))
        title_row2 += models[toxidrome]
    title_row1 += other_toxidromes
    title_row2 += other_toxidromes

    # the '_ = writer...' ensures the return of writerow does not 
    # get written into stdout. The writer can be writing to stdout 
    _ = writer.writerow(title_row1) 
    _ = writer.writerow(title_row2)

    for datapass in datapassList:
        newRow = [datapass.name, datapass.original]
        for toxidrome in main_toxidromes:
            for model in models[toxidrome]:
                value, accuracy = datapass.getToxidromeData(toxidrome, model)
                if sd:
                    newRow.append(accuracy)
                else:
                    newRow.append(value)
        for toxidrome in other_toxidromes:
            newRow.append(datapass.getToxidromeValue(toxidrome))
        _ = writer.writerow(newRow)


def isOrganicAndHasWeight(smiles, weight):
    mol = Chem.MolFromSmiles(smiles)
    if(mol == None):
        return False
    mw = Descriptors.MolWt(mol)
    # Check for inorganic atoms
    inorganic = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ["H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]:
            inorganic = True
            break
            
    # Add the row to the filtered DataFrame if it passes all criteria
    return not inorganic and mw >= weight
