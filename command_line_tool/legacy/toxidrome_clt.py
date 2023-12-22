import argparse
import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
import pandas as pd
import math
from cmpnn_toxidrome.run_cmpnn import run_cmpnn
from dmpnn_toxidrome.run_dmpnn import run_dmpnn
from chembl_structure_pipeline import standardizer
import base64
from io import BytesIO

import time
import json

NAME_TITLE = 'ID'
SMILE_TITLE = 'SMILES'
NO_RESULT_TEXT = '(None)'

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", nargs=1, help="File location of input csv")
    parser.add_argument("-ln", "--names", nargs=1, help="delimited string of compound names")
    parser.add_argument("-ls", "--smiles", nargs=1, help="delimited string of compound smiles")
    parser.add_argument("-o", "--output", nargs=1, help="File location of where to save output csv")
    parser.add_argument("-p", "--print", action='store_true', help="print the output in json")
    parser.add_argument("-nnType", nargs=1, help="Type of Neural network to run. Either \"cmpnn\" or \"dmpnn\"")
    
    args = parser.parse_args()

    if args.input == None and args.smiles == None:
        parser.error("Either an input file or a list of smiles must be provided.")
    if args.input != None and args.smiles != None:
        parser.error("Cannot provide both an input file and a list of smiles")
    if args.input != None and args.names != None:
        parser.error("Cannot provide custom compound names when providing an input file")
    names = []
    smiles = []
    if args.output == None and not args.print:
        parser.error("No valid output, use at least one of \"-o\" (file location of output csv) and \"-p\" (print output)")
    
    if args.names != None:
        names = args.names[0].split(",")
    if args.smiles != None:
        smiles = args.smiles[0].split(",")
    else:
        names = None
    
    if names != None and len(names) != len(smiles):
        parser.error("Unequal number of names and smiles provided")
    

    if args.nnType == None:
        args.nnType = ["dmpnn"]
    else:
        if args.nnType[0] != "cmpnn" and args.nnType[0] != "dmpnn":
            parser.error(f"Invalid model neural network \"{args.nnType[0]}\". Available neural network types for --nnType: cmpnn, dmpnn")

    if smiles == []:
        names, smiles = parse_input_csv(args.input[0])

    return (names, smiles, args.output, args.print, args.nnType[0])

class molecule:
    def __init__(self, name, smiles):
        self.name = name
        self.smiles = smiles
        self.fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2, nBits=2048)
        self.toxidromes = {}
        self.toxidromes["Opioid"] = {}
        self.toxidromes["Sympathomimetic"] = {}
        self.toxidromes["Convulsant"] = {}
        self.toxidromes["Cholinergic"] = {}
        self.toxidromes["Anticholinergic"] = {}
        self.accuracy = {}
        self.accuracy["Opioid"] = {}
        self.accuracy["Sympathomimetic"] = {}
        self.accuracy["Convulsant"] = {}
        self.accuracy["Cholinergic"] = {}
        self.accuracy["Anticholinergic"] = {}
        #self.img = gen_mol_img(smiles)

    def setToxidrome(self, toxidrome, l):
        if toxidrome == "Opioid":
            self.setOpioid(l)
        elif toxidrome == "Sympathomimetic":
            self.setSympathomimetic(l)
        elif toxidrome == "Convulsant":
            self.setConvulsant(l)
        elif toxidrome == "Cholinergic":
            self.setCholinergic(l)

    def setToxidromeAccuracy(self, toxidrome, l):
        if toxidrome == "Opioid":
            self.setOpioidAccuracy(l)
        elif toxidrome == "Sympathomimetic":
            self.setSympathomimeticAccuracy(l)
        elif toxidrome == "Convulsant":
            self.setConvulsantAccuracy(l)
        elif toxidrome == "Cholinergic":
            self.setCholinergicAccuracy(l)

    def setOpioid(self, l):
        self.toxidromes["Opioid"]["pEC50_Mu"] = l[3]
        self.toxidromes["Opioid"]["pEC50_Delta"] = l[4]
        self.toxidromes["Opioid"]["pEC50_Kappa"] = l[5]

    def setSympathomimetic(self, l):
        self.toxidromes["Sympathomimetic"]["NET_pIC50_rat"] = l[0]
        self.toxidromes["Sympathomimetic"]["DAT_pIC50_rat"] = l[1]
        self.toxidromes["Sympathomimetic"]["DAT_pIC50_human"] = l[2]
        self.toxidromes["Sympathomimetic"]["B1_pEC50_human"] = l[3]
        self.toxidromes["Sympathomimetic"]["B3_pEC50_human"] = l[4]
        self.toxidromes["Sympathomimetic"]["A1a_pEC50_human"] = l[5]
        self.toxidromes["Sympathomimetic"]["B2_pEC50_human"] = l[6]

    def setConvulsant(self, l):
        self.toxidromes["Convulsant"]["pEC50_GluN_e2z1_rat"] = l[0]
        self.toxidromes["Convulsant"]["pEC50_GluN_e3z1_rat"] = l[1]
        self.toxidromes["Convulsant"]["pEC50_GlyT1_human"] = l[2]
        self.toxidromes["Convulsant"]["pEC50_GlyT2_human"] = l[3]
        self.toxidromes["Convulsant"]["pIC50_GABAa1b2g2_human"] = l[4]
        self.toxidromes["Convulsant"]["pEC50_GluA_human"] = l[5]

    def setCholinergic(self, l):
        self.toxidromes["Cholinergic"]["pEC50_rM1"] = l[0]
        self.toxidromes["Cholinergic"]["pEC50_hM1"] = l[1]
        self.toxidromes["Cholinergic"]["pEC50_hM2"] = l[2]
        self.toxidromes["Cholinergic"]["pEC50_hM4"] = l[3]
        self.toxidromes["Cholinergic"]["pEC50_rM4"] = l[4]
        self.toxidromes["Cholinergic"]["pEC50_hM5"] = l[5]
        self.toxidromes["Anticholinergic"]["pIC50_hM2"] = l[9]
        self.toxidromes["Anticholinergic"]["pIC50_rM2"] = l[10]
        self.toxidromes["Anticholinergic"]["pIC50_hM3"] = l[11]
        self.toxidromes["Anticholinergic"]["pIC50_hM4"] = l[12]
        self.toxidromes["Anticholinergic"]["pIC50_hM5"] = l[13]
        self.toxidromes["Anticholinergic"]["pIC50_hM1"] = l[16]
        self.toxidromes["Anticholinergic"]["pIC50_rM1"] = l[17]
        self.toxidromes["Cholinergic"]["pIC50_AChE_HouseMouse"] = l[21]
        self.toxidromes["Cholinergic"]["pIC50_AChE_human"] = l[22]
        self.toxidromes["Cholinergic"]["pIC50_AChE_Cow"] = l[23]
        self.toxidromes["Cholinergic"]["pIC50_AChE_eel"] = l[24]
        self.toxidromes["Cholinergic"]["pIC50_AChE_BrownRat"] = l[25]
        self.toxidromes["Cholinergic"]["pIC50_BChE_human"] = l[26]

    def setOpioidAccuracy(self, l):
        self.accuracy["Opioid"]["pEC50_Mu"] = l[3]
        self.accuracy["Opioid"]["pEC50_Delta"] = l[4]
        self.accuracy["Opioid"]["pEC50_Kappa"] = l[5]

    def setSympathomimeticAccuracy(self, l):
        self.accuracy["Sympathomimetic"]["NET_pIC50_rat"] = l[0]
        self.accuracy["Sympathomimetic"]["DAT_pIC50_rat"] = l[1]
        self.accuracy["Sympathomimetic"]["DAT_pIC50_human"] = l[2]
        self.accuracy["Sympathomimetic"]["B1_pEC50_human"] = l[3]
        self.accuracy["Sympathomimetic"]["B3_pEC50_human"] = l[4]
        self.accuracy["Sympathomimetic"]["A1a_pEC50_human"] = l[5]
        self.accuracy["Sympathomimetic"]["B2_pEC50_human"] = l[6]

    def setConvulsantAccuracy(self, l):
        self.accuracy["Convulsant"]["pEC50_GluN_e2z1_rat"] = l[0]
        self.accuracy["Convulsant"]["pEC50_GluN_e3z1_rat"] = l[1]
        self.accuracy["Convulsant"]["pEC50_GlyT1_human"] = l[2]
        self.accuracy["Convulsant"]["pEC50_GlyT2_human"] = l[3]
        self.accuracy["Convulsant"]["pIC50_GABAa1b2g2_human"] = l[4]
        self.accuracy["Convulsant"]["pEC50_GluA_human"] = l[5]

    def setCholinergicAccuracy(self, l):
        self.accuracy["Cholinergic"]["pEC50_rM1"] = l[0]
        self.accuracy["Cholinergic"]["pEC50_hM1"] = l[1]
        self.accuracy["Cholinergic"]["pEC50_hM2"] = l[2]
        self.accuracy["Cholinergic"]["pEC50_hM4"] = l[3]
        self.accuracy["Cholinergic"]["pEC50_rM4"] = l[4]
        self.accuracy["Cholinergic"]["pEC50_hM5"] = l[5]
        self.accuracy["Anticholinergic"]["pIC50_hM2"] = l[9]
        self.accuracy["Anticholinergic"]["pIC50_rM2"] = l[10]
        self.accuracy["Anticholinergic"]["pIC50_hM3"] = l[11]
        self.accuracy["Anticholinergic"]["pIC50_hM4"] = l[12]
        self.accuracy["Anticholinergic"]["pIC50_hM5"] = l[13]
        self.accuracy["Anticholinergic"]["pIC50_hM1"] = l[16]
        self.accuracy["Anticholinergic"]["pIC50_rM1"] = l[17]
        self.accuracy["Cholinergic"]["pIC50_AChE_HouseMouse"] = l[21]
        self.accuracy["Cholinergic"]["pIC50_AChE_human"] = l[22]
        self.accuracy["Cholinergic"]["pIC50_AChE_Cow"] = l[23]
        self.accuracy["Cholinergic"]["pIC50_AChE_eel"] = l[24]
        self.accuracy["Cholinergic"]["pIC50_AChE_BrownRat"] = l[25]
        self.accuracy["Cholinergic"]["pIC50_BChE_human"] = l[26]
    
    def setAnticoagulant(self, p):
        self.toxidromes["Anticoagulant"] = p
    def setIrritantCorrosive(self, p):
        self.toxidromes["Irritant, corrosive"] = p
    def setKnockdown(self, p):
        self.toxidromes["Knockdown"] = p
    def setSolventsAnestheticsSedatives(self, p):
        self.toxidromes["Solvents, anesthetics, sedatives"] = p

def gen_mol_img(smiles):
    img_pil = Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(600,600))
    buffered = BytesIO()
    img_pil.save(buffered, format="JPEG")
    img_bytes = buffered.getvalue()
    img_b64 = base64.b64encode(img_bytes)

    n=500
    strParts = [img_b64[i:i+n] for i in range(0, len(img_b64), n)]

    return img_b64.decode("utf-8")


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
        if len(names) > 0 and names[0] == NAME_TITLE:
            names.pop(0)
            smiles.pop(0)
        elif smiles[0] == SMILE_TITLE:
            smiles.pop(0)
            names = None
        else:
            raise Exception("Invalid CSV format")
        return names, smiles
    

   
def getToxidromeMols():
    toxidromeMols = pd.read_excel('app/toxidrome_mols.xlsx') # Format: Toxidrome, SMILES
    return toxidromeMols.groupby('Toxidrome')['Molecule'].agg(list).to_dict()

def getToxidromeFPs(toxidromeMols):
    toxidromeFPs = {}
    for toxidrome in toxidromeMols:
        toxidromeFPs[toxidrome] = []
        for smiles in toxidromeMols[toxidrome]:
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2, nBits=2048)
            toxidromeFPs[toxidrome].append(fp)
    return toxidromeFPs

def standardize(smiles_list):
    count = len(smiles_list)
    smiles_std_list = [None] * count
    for index in range(0,count):
        smiles = smiles_list[index]
        if smiles is not None:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is not None:
                m_no_salt = standardizer.get_parent_mol(molecule)
                molecule_std = standardizer.standardize_mol(m_no_salt[0])
                if molecule_std is not None:
                    smiles_std = Chem.MolToSmiles(molecule_std)
                    smiles_std_list[index] = smiles_std
    return smiles_std_list

def calcP(rs, n):
    z = (rs - (0.000424 * n)) / (0.00449 * (n ** 0.665))
    xz = -(math.e ** ((-1 * z * math.pi) / (math.sqrt(6) - 0.577215665)))
    if z <= 28:
        return (1 - (math.e ** xz))
    else:
        return ((-1 * xz) - ((xz ** 2) / 2) - ((xz ** 3) / 6))

def run_sea(molecules):
    toxidromeMols = getToxidromeMols()
    toxidromeFPs = getToxidromeFPs(toxidromeMols)

    for mol in molecules:
        a_rs = 0
        for toxFP in toxidromeFPs["Anticoagulant"]:
            ts = DataStructs.FingerprintSimilarity(mol.fp, toxFP, metric=DataStructs.TanimotoSimilarity)
            if ts >= 0.57:
                a_rs += ts
        mol.setAnticoagulant(round(calcP(a_rs, len(toxidromeMols["Anticoagulant"])), 5))

        i_rs = 0
        for toxFP in toxidromeFPs["Irritant, Corrosive"]:
            ts = DataStructs.FingerprintSimilarity(mol.fp, toxFP, metric=DataStructs.TanimotoSimilarity)
            if ts >= 0.57:
                i_rs += ts
        mol.setIrritantCorrosive(round(calcP(i_rs, len(toxidromeMols["Irritant, Corrosive"])), 5))

        k_rs = 0
        for toxFP in toxidromeFPs["Knockdown"]:
            ts = DataStructs.FingerprintSimilarity(mol.fp, toxFP, metric=DataStructs.TanimotoSimilarity)
            if ts >= 0.57:
                k_rs += ts
        mol.setKnockdown(round(calcP(k_rs, len(toxidromeMols["Knockdown"])), 5))

        s_rs = 0
        for toxFP in toxidromeFPs["Solvents, anesthetics, sedatives"]:
            ts = DataStructs.FingerprintSimilarity(mol.fp, toxFP, metric=DataStructs.TanimotoSimilarity)
            if ts >= 0.57:
                s_rs += ts
        mol.setSolventsAnestheticsSedatives(round(calcP(s_rs, len(toxidromeMols["Solvents, anesthetics, sedatives"])), 5))




def run_main_models(mols, nnType):

    smiles = []
    for m in mols:
        smiles.append(m.smiles)

    all_preds = {}
    all_acc = {}

    if nnType == "dmpnn":
        all_preds, all_acc = run_dmpnn(smiles)
    elif nnType == "cmpnn":
        all_preds, all_acc = run_cmpnn(smiles)
        
    for toxidrome in all_preds:
        for m in range(len(mols)):
            mols[m].setToxidrome(toxidrome, all_preds[toxidrome][m])
            mols[m].setToxidromeAccuracy(toxidrome, all_acc[toxidrome][m])

def stan_dev_to_error(mols):
    for mol in mols:
        for toxidrome in mol.accuracy.keys():
            for model in mol.accuracy[toxidrome].keys():
                mol.accuracy[toxidrome][model] = (1.186 * mol.accuracy[toxidrome][model]) + 0.252
        
def write_output(mols, outputfile, sd):
    with open(outputfile, 'w', newline='') as file:
        writer = csv.writer(file)

        models = {
            "Opioid" :          ["pEC50_Mu", "pEC50_Delta", "pEC50_Kappa"],
            "Sympathomimetic" : ["NET_pIC50_rat", "DAT_pIC50_rat", "DAT_pIC50_human", "B1_pEC50_human", "B3_pEC50_human", "A1a_pEC50_human", "B2_pEC50_human"],
            "Convulsant" :      ["pEC50_GluN_e2z1_rat", "pEC50_GluN_e3z1_rat", "pEC50_GlyT1_human", "pEC50_GlyT2_human", "pIC50_GABAa1b2g2_human", "pEC50_GluA_human"],
            "Cholinergic" :     ["pEC50_rM1", "pEC50_hM1", "pEC50_hM2", "pEC50_hM4", "pEC50_rM4", "pEC50_hM5", "pIC50_AChE_HouseMouse", "pIC50_AChE_human", "pIC50_AChE_Cow", "pIC50_AChE_eel", "pIC50_AChE_BrownRat", "pIC50_BChE_human"],
            "Anticholinergic" : ["pIC50_hM2", "pIC50_rM2", "pIC50_hM3", "pIC50_hM4", "pIC50_hM5", "pIC50_hM1", "pIC50_rM1"]
        }
        main_toxidromes = ["Opioid", "Sympathomimetic", "Convulsant", "Cholinergic", "Anticholinergic"]
        other_toxidromes = ["Anticoagulant", "Irritant, corrosive", "Knockdown", "Solvents, anesthetics, sedatives"]

        title_row1 = [NAME_TITLE,SMILE_TITLE]
        title_row2 = [NAME_TITLE,SMILE_TITLE]
        for toxidrome in main_toxidromes:
            title_row1 += ([toxidrome] * len(models[toxidrome]))
            title_row2 += models[toxidrome]
        title_row1 += other_toxidromes
        title_row2 += other_toxidromes

        writer.writerow(title_row1)
        writer.writerow(title_row2)

        for mol in mols:
            newRow = [mol.name, mol.smiles]
            for toxidrome in main_toxidromes:
                for model in models[toxidrome]:
                    if sd:
                        newRow.append(mol.accuracy[toxidrome][model])
                    else:
                        newRow.append(mol.toxidromes[toxidrome][model])
            for toxidrome in other_toxidromes:
                newRow.append(mol.toxidromes[toxidrome])
            writer.writerow(newRow)
        
def build_json(mols):
    job_dict = {}
    for i in range(len(mols)):
        job_dict[i] = {}
        job_dict[i]["name"] = mols[i].name
        job_dict[i]["smiles"] = mols[i].smiles
        #job_dict[i]["img"] = mols[i].img
        job_dict[i]["toxidromes"] = {}
        for toxidrome in mols[i].toxidromes:
            if isinstance(mols[i].toxidromes[toxidrome], float):
                job_dict[i]["toxidromes"][toxidrome] = {
                    None: {
                        "val": mols[i].toxidromes[toxidrome],
                        "acc": None
                    }
                }
            else:
                job_dict[i]["toxidromes"][toxidrome] = {}
                for model in mols[i].toxidromes[toxidrome]:
                    job_dict[i]["toxidromes"][toxidrome][model] = {}
                    job_dict[i]["toxidromes"][toxidrome][model]["val"] = mols[i].toxidromes[toxidrome][model]
                    job_dict[i]["toxidromes"][toxidrome][model]["acc"] = mols[i].accuracy[toxidrome][model]

    sys.stdout.write(json.dumps(job_dict))
    sys.stdout.flush() 
        
def main(names, smiles, output_file, doPrint, nnType):#, is_json = False):

    smiles = standardize(smiles)

    mols = []
    for i in range(len(smiles)):
        n = ""
        if names == None:
            n = f"Compound {i + 1}"
        else:
            n = names[i]
        mols.append(molecule(n, smiles[i]))
        
    run_sea(mols)

    run_main_models(mols, nnType)

    stan_dev_to_error(mols)

    if doPrint:
        build_json(mols)
    if output_file != None:
        write_output(mols, output_file[0], False)
 
if __name__ == "__main__":
    names, smiles, output_file, doPrint, nnType = parse_args()
    main(names, smiles, output_file, doPrint, nnType)