import math
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import pandas as pd

SMILES = "C1C(CC2=CC=CC=C2C1C3=C(OC4=CC=CC=C4C3=O)O)C5=CC=C(C=C5)C6=CC=C(C=C6)Br" #Place SMILES you want to run SEA on here.

ROUND_DIGIT = 2
TOX_Anticoagulant = "Anticoagulant"
TOX_SAS = "Solvents, anesthetics, sedatives"
TOX_IrritantCorrosive = "Irritant, corrosive"
TOX_Knockdown = "Knockdown"
SEA_TOXIDROMES = [TOX_Anticoagulant, TOX_IrritantCorrosive, TOX_Knockdown, TOX_SAS]

def run_sea(mol_fingerprint, toxidromeMols, toxidromeFPs):
	if mol_fingerprint == None:
		return # no need to do anything if FP is None
	for toxidrome_sea in SEA_TOXIDROMES:
		rs = 0
		for toxFP in toxidromeFPs[toxidrome_sea]:
			ts = DataStructs.FingerprintSimilarity(mol_fingerprint, toxFP, metric=DataStructs.TanimotoSimilarity)
			if ts >= 0.57:
				rs += ts
		roundedValue = calcP(rs, len(toxidromeMols[toxidrome_sea]))
		print(f"{toxidrome_sea}: {roundedValue}")
		
def calcP(rs, n):
    z = (rs - (0.000424 * n)) / (0.00449 * (n ** 0.665))
    xz = -(math.e ** ((-1 * z * math.pi) / (math.sqrt(6) - 0.577215665)))
    if z <= 28:
        return (1 - (math.e ** xz))
    else:
        return ((-1 * xz) - ((xz ** 2) / 2) - ((xz ** 3) / 6))
	
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

if __name__ == "__main__":
      toxidromeMols, toxidromeFPs = loadToxidromeData_forSEA("sea_mols.xlsx")
      fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(SMILES), 2, nBits=2048)
      run_sea(fp, toxidromeMols, toxidromeFPs)
