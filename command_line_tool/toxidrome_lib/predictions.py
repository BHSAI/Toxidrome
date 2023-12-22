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
import os
import time
import json

from toxidrome_lib.constants import *
from toxidrome_lib.util import *


# This class contains
# 1. one compound input
# 2. intermediate data like standardized smiles and FP
# 3. outputs
class PredictionPass:
	def __init__(self, name, original):
		self.name = name
		self.original = original
		self.smiles = None
		self.fp = None
		# create a mapping from a toxidrome to data which can be
		# 1. another mapping based on subcategories (models)
		# 2. value (numeric)
		self.toxidromes = {toxidrome: {} for toxidrome in TOX_CATEGORY.keys()}
		self.accuracies = {toxidrome: {} for toxidrome in TOX_CATEGORY.keys()}
		self.appDomain = {toxidrome: {} for toxidrome in TOX_CATEGORY.keys()}

	def setToxidromeValue(self, toxidrome, value):
		self.toxidromes[toxidrome] = value

	def setToxidromeAD(self, toxidrome, category, applicabilityDomain):
		if category == None:
			self.appDomain[toxidrome] = applicabilityDomain
		else:
			self.appDomain[toxidrome][category] = applicabilityDomain

	def setToxidromeData(self, toxidrome, category, value, accuracy):
		self.toxidromes[toxidrome][category] = value
		self.accuracies[toxidrome][category] = accuracy

	def getToxidromeValue(self, toxidrome):
		if toxidrome in self.toxidromes.keys():
			return self.toxidromes[toxidrome]
		return None
	
	def getToxidromeAccuracy(self, toxidrome, category):
		if toxidrome in self.accuracies.keys():
			if category in self.accuracies[toxidrome].keys():
				return self.accuracies[toxidrome][category]
		return None
	
	def getToxidromeData(self, toxidrome, category):
		if toxidrome in self.toxidromes.keys():
			categoryToData = self.toxidromes[toxidrome]
			if category in categoryToData.keys():
				return categoryToData[category], self.accuracies[toxidrome][category]
		return None, None

def run_sea_one_compound(datapass, toxidromeMols, toxidromeFPs):
	if datapass.fp == None:
		return # no need to do anything if FP is None
	for toxidrome_sea in SEA_TOXIDROMES:
		rs = 0
		for toxFP in toxidromeFPs[toxidrome_sea]:
			ts = DataStructs.FingerprintSimilarity(datapass.fp, toxFP, metric=DataStructs.TanimotoSimilarity)
			if ts >= 0.57:
				rs += ts
		roundedValue = round(calcP(rs, len(toxidromeMols[toxidrome_sea])), ROUND_DIGIT)
		datapass.setToxidromeValue(toxidrome_sea, roundedValue)


def run_main_models_multiple_compounds(datapassList, nnType):
	#the input requires an array
	smilesList = [datapass.smiles for datapass in datapassList] 

	all_preds = {}
	all_acc = {}

	if nnType == "dmpnn":
		all_preds, all_acc = run_dmpnn(smilesList)
	elif nnType == "cmpnn":
		all_preds, all_acc = run_cmpnn(smilesList)
		

	'''
	all_preds and all_acc have the following structure:
	toxidrome_model -> 
	array:
		index: corresponding to an input smiles
		value: predictions (an array of numbers)	

	'''
	for toxidrome_model in all_preds.keys():
		pred_model = all_preds[toxidrome_model]
		accu_model = all_acc[toxidrome_model]
		#DATA_POSITION_TOX contains the info mapping from toxidrome_model to 
		#index to [toxidrome, subcategory]
		indexToPositionData = DATA_POSITION_TOX[toxidrome_model]
		for input_index in range(len(pred_model)):
			predictions = pred_model[input_index]
			accuracies = accu_model[input_index]
			datapass = datapassList[input_index]
			for pred_index in range(len(predictions)): 
				if not pred_index in indexToPositionData.keys():
					continue # the element at the pred_index is not used
				positionData = indexToPositionData[pred_index]
				toxidrome = positionData[0]
				category = positionData[1]
				# assuming either number or None
				if predictions[pred_index] == None:
					pred = None
				else:
					pred = round(predictions[pred_index], ROUND_DIGIT)
				if accuracies[pred_index] == None:
					accu = None
				else:
					accu = round(accuracies[pred_index], ROUND_DIGIT)
				datapass.setToxidromeData(toxidrome, category, pred, accu)

def run_main_models_one_compound(datapass, nnType):
	#the input requires an array
	smilesList = [datapass.smiles]

	all_preds = {}
	all_acc = {}

	if nnType == "dmpnn":
		all_preds, all_acc = run_dmpnn(smilesList)
	elif nnType == "cmpnn":
		all_preds, all_acc = run_cmpnn(smilesList)

	'''
	all_preds and all_acc have the following structure:
	toxidrome_model -> 
	array:
		index: corresponding to an input smiles
		value: predictions (an array of numbers)	

	'''
	for toxidrome_model in all_preds.keys():
		predictions = all_preds[toxidrome_model][0]
		accuracies = all_acc[toxidrome_model][0]
		if not toxidrome_model in DATA_POSITION_TOX.keys():
			# might want to raise exception or log it
			# does not expect this
			# can be name mismatching at data file level
			# with whats inside the constants.py
			continue
		indexToPositionData = DATA_POSITION_TOX[toxidrome_model]
		# not every element in the array is used
		for index in range(len(predictions)):
			if not index in indexToPositionData.keys():
				continue # the element at the index is not used
			positionData = indexToPositionData[index] # [toxidrome, category]
			toxidrome = positionData[0]
			category = positionData[1]
			if predictions[pred_index] == None:
				pred = None
			else:
				pred = round(predictions[pred_index], ROUND_DIGIT)
			if accuracies[pred_index] == None:
				accu = None
			else:
				accu = round(accuracies[pred_index], ROUND_DIGIT)
			datapass.setToxidromeData(toxidrome, category, pred, accu)

def stan_dev_to_error_one_compound(datapass):
	for toxidrome in datapass.accuracies.keys():
		for subcategory in datapass.accuracies[toxidrome].keys():
			value = datapass.accuracies[toxidrome][subcategory]
			if value != None:
				errorValue = (1.186 * datapass.accuracies[toxidrome][subcategory]) + 0.252
				datapass.accuracies[toxidrome][subcategory] = errorValue