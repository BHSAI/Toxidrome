'''
This file contains constants and configuration used by the python scripts at the same directory level
'''

# A line that is printed before the error log
# ERROR_START_LINE = '---------------------error log-----------------------'
# MAX_ERROR_LOG_LINE = 100
ROUND_DIGIT = 2
ORGANIC_COMPOUND_MIN_WEIGHT = 50

# Constants used in the CSV output file
NAME_TITLE = 'Name'
SMILE_TITLE = 'SMILES'
NO_RESULT_TEXT = '(None)'

OUTPUT_CHOICE = ['csv', 'json', 'database']

#At the time of writing
#The order of the following list matches the order presented in the web
#TOX_SAS is 'Acute Exposure'
TOX_Opioid = "Opioid"
TOX_Sympathomimetic = "Sympathomimetic"
TOX_Convulsant= "Convulsant"
TOX_Cholinergic = "Cholinergic"
TOX_Anticholinergic = "Anticholinergic"
TOX_Anticoagulant = "Anticoagulant"
TOX_SAS = "Solvents, anesthetics, sedatives"
TOX_IrritantCorrosive = "Irritant, corrosive"
TOX_Knockdown = "Knockdown"

# This is used to assign prediction results into the supporting model.
# It is a 2D mapping from toxidromes to subcategories to 
# [toxidrome_model, index] where toxidrome_model could differ from the 
# toxidrome. e.g. data for Anticholinergic will be found under 
# toxidrome_model 'Cholinergic'
# The toxidrome_model list can be found under dmpnn_toxidrome/toxidromes

# The make_prediction method returns the following structure
# 	toxidrome_model -> 
# 	array:
# 		index: corresponding to an input smiles
# 		value: predictions (an array of numbers)
# The number below corresponds to the element in the 'predictions' array 	

# The order of the subcategories is preserved though the list structure.
TOX_CATEGORY_DATA_POSITION ={
	# TOX_Opioid:          3-5
	# TOX_Sympathomimetic: 0-6
	# TOX_Convulsant:      0-5 
	TOX_Opioid :          [{ "pEC50_Mu":               [TOX_Opioid, 3]},
						   { "pEC50_Delta":            [TOX_Opioid, 4]},
						   { "pEC50_Kappa":            [TOX_Opioid, 5]}],
	TOX_Sympathomimetic : [{ "NET_pIC50_rat":          [TOX_Sympathomimetic, 0]},
						   { "DAT_pIC50_rat":          [TOX_Sympathomimetic, 1]},
						   { "DAT_pIC50_human":        [TOX_Sympathomimetic, 2]},
						   { "B1_pEC50_human":         [TOX_Sympathomimetic, 3]},
						   { "B3_pEC50_human":         [TOX_Sympathomimetic, 4]},
						   { "A1a_pEC50_human":        [TOX_Sympathomimetic, 5]},
						   { "B2_pEC50_human":         [TOX_Sympathomimetic, 6]}],
	TOX_Convulsant :      [{ "pEC50_GluN_e2z1_rat":    [TOX_Convulsant, 0]},
						   { "pEC50_GluN_e3z1_rat":    [TOX_Convulsant, 1]},
						   { "pEC50_GlyT1_human":      [TOX_Convulsant, 2]},
						   { "pEC50_GlyT2_human":      [TOX_Convulsant, 3]},
						   { "pIC50_GABAa1b2g2_human": [TOX_Convulsant, 4]},
						   { "pEC50_GluA_human":       [TOX_Convulsant, 5]}],
	# data for TOX_Cholinergic and TOX_Anticholinergic are in one array
	# TOX_Cholinergic:     0-5,21-26
	# TOX_Anticholinergic: 9-13,16,17
	TOX_Cholinergic :     [{ "pEC50_rM1":              [TOX_Cholinergic, 0]},
						   { "pEC50_hM1":              [TOX_Cholinergic, 1]},
						   { "pEC50_hM2":              [TOX_Cholinergic, 2]},
						   { "pEC50_hM4":              [TOX_Cholinergic, 3]},
						   { "pEC50_rM4":              [TOX_Cholinergic, 4]},
						   { "pEC50_hM5":              [TOX_Cholinergic, 5]},
						   { "pIC50_AChE_HouseMouse":  [TOX_Cholinergic, 21]},
						   { "pIC50_AChE_human":       [TOX_Cholinergic, 22]},
						   { "pIC50_AChE_Cow":         [TOX_Cholinergic, 23]},
						   { "pIC50_AChE_eel":         [TOX_Cholinergic, 24]},
						   { "pIC50_AChE_BrownRat":    [TOX_Cholinergic, 25]},
						   { "pIC50_BChE_human":       [TOX_Cholinergic, 26]}],
	TOX_Anticholinergic : [{ "pIC50_hM2":              [TOX_Cholinergic, 9]},
						   { "pIC50_rM2":              [TOX_Cholinergic, 10]},
						   { "pIC50_hM3":              [TOX_Cholinergic, 11]},
						   { "pIC50_hM4":              [TOX_Cholinergic, 12]},
						   { "pIC50_hM5":              [TOX_Cholinergic, 13]},
						   { "pIC50_hM1":              [TOX_Cholinergic, 16]},
						   { "pIC50_rM1":              [TOX_Cholinergic, 17]}]
}

TOX_CATEGORY = {key : [list(row.keys())[0] \
	for row in TOX_CATEGORY_DATA_POSITION[key]] \
	for key in TOX_CATEGORY_DATA_POSITION.keys() }

# The straight forward assignment
# TOX_CATEGORY = {
#     TOX_Opioid :          ["pEC50_Mu", "pEC50_Delta", "pEC50_Kappa"],
#     TOX_Sympathomimetic : ["NET_pIC50_rat", "DAT_pIC50_rat", "DAT_pIC50_human", "B1_pEC50_human", "B3_pEC50_human", "A1a_pEC50_human", "B2_pEC50_human"],
#     TOX_Convulsant :      ["pEC50_GluN_e2z1_rat", "pEC50_GluN_e3z1_rat", "pEC50_GlyT1_human", "pEC50_GlyT2_human", "pIC50_GABAa1b2g2_human", "pEC50_GluA_human"],
#     TOX_Cholinergic :     ["pEC50_rM1", "pEC50_hM1", "pEC50_hM2", "pEC50_hM4", "pEC50_rM4", "pEC50_hM5", "pIC50_AChE_HouseMouse", "pIC50_AChE_human", "pIC50_AChE_Cow", "pIC50_AChE_eel", "pIC50_AChE_BrownRat", "pIC50_BChE_human"],
#     TOX_Anticholinergic : ["pIC50_hM2", "pIC50_rM2", "pIC50_hM3", "pIC50_hM4", "pIC50_hM5", "pIC50_hM1", "pIC50_rM1"]
# }

def __revser_map__ (map):
	ret = {}
	for toxidrome in TOX_CATEGORY_DATA_POSITION.keys():
		for row in TOX_CATEGORY_DATA_POSITION[toxidrome]:
			# only 1 item per row
			subcategory = list(row.keys())[0]
			model_index = list(row.values())[0]
			model = model_index[0]
			index = model_index[1]
			#print(f'{model} {index} {toxidrome} {subcategory}')
			if not model in ret.keys():
				ret[model] = {}
			ret[model][index] = [toxidrome, subcategory]
	return ret

# a mapping from toxidrome_model, index => [toxidrome, subcategory]
DATA_POSITION_TOX = __revser_map__(TOX_CATEGORY_DATA_POSITION)


MAIN_TOXIDROMES = [TOX_Opioid, TOX_Sympathomimetic, TOX_Convulsant, TOX_Cholinergic, TOX_Anticholinergic]
SEA_TOXIDROMES = [TOX_Anticoagulant, TOX_IrritantCorrosive, TOX_Knockdown, TOX_SAS]

SEA_THRESHOLD = 0.05
AD_CUTOFF = 1


