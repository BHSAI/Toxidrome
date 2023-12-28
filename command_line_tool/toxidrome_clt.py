import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__))) 

import argparse
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs

import math
from cmpnn_toxidrome.run_cmpnn import run_cmpnn
from dmpnn_toxidrome.run_dmpnn import run_dmpnn
from chembl_structure_pipeline import standardizer
import base64
from io import BytesIO
import time
import json
import traceback
import contextlib

from toxidrome_lib.constants import *
from toxidrome_lib.predictions import *
from toxidrome_lib.util import *


'''
This file serves as the entry point and master loop which does
1. parse inputs
2. prepare data
    a. standardization
    b. FP calculation
    c. read data file
3. call prediction procedures
4. create outputs
Each pass handles one compound.
'''
class InputBundle:
    names = None         # -ln      list
    smiles = None        # -ls      list
    outputPath = None    # -o       string?
    fileFormat = 'csv'   # -f       string
    errorPath = None     # -e       string?
    statePath = None     # -sf      string?
    mode = 'performance' # -m       string
    bulk = None
    nnType = 'dmpnn'     # -nnType  string

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs=1, help="File location of input csv")
    parser.add_argument("-ln", "--names", nargs=1, help="delimited string of compound names")
    parser.add_argument("-ls", "--smiles", nargs=1, help="delimited string of compound smiles")
    parser.add_argument("-o", "--output", nargs=1, help="Output file path")
    parser.add_argument("-f", "--format", nargs=1, help="File format for the output. either csv, json or database. default: csv")
    parser.add_argument("-e", "--error", nargs=1, help="The file path for error messages")
    parser.add_argument("-sf", "--state", nargs=1, help="The file path for outputting state info.")
    parser.add_argument("-m", "--mode", nargs=1, help="The mode of execution. 'performance' or 'robust'. default: 'performance'")
    parser.add_argument("-b", "--bulk", nargs=1, help="The amount of compounds that go through Main Predictions each time in order to avoid any memory issue.")
    parser.add_argument("-nnType", nargs=1, help="Type of Neural network to run. Either \"cmpnn\" or \"dmpnn\". default: dmpnn")
    args = parser.parse_args()

    inputs = InputBundle()
    # if no output path is provided, the output goes to the stdout

    # ensure smiles are provided either through cmd or an input file
    if args.input == None and args.smiles == None:
        parser.error("Either an input file or a list of smiles must be provided.")
    if args.input != None and args.smiles != None:
        parser.error("Cannot provide both an input file and a list of smiles")
    if args.input != None and args.names != None:
        parser.error("Cannot provide custom compound names when providing an input file")
    # Process input smiles and names
    if args.names != None:
        inputs.names = args.names[0].split(",")
    if args.smiles != None:
        inputs.smiles = args.smiles[0].split(",")
    if inputs.names != None and len(inputs.names) != len(inputs.smiles):
        parser.error("Unequal number of names and smiles provided")
    if args.smiles == None: # not provided through cmd
        inputs.names, inputs.smiles = parse_input_csv(args.input[0])
    if inputs.names == None or len(inputs.names) == 0:
        inputs.names = ['Compound '+str(i) for i in range(0, len(inputs.smiles))]
    if inputs.smiles == None or len(inputs.smiles) == 0:
        parser.error("Empty input file.")

    # Process paths for outputs
    if args.output != None and len(args.output) > 0:
        inputs.outputPath = args.output[0]
    if args.error != None and len(args.error) > 0:
        inputs.errorPath = args.error[0]
    if args.state != None and len(args.state) > 0:
        inputs.statePath = args.state[0]

    # Process configurations
    if args.format != None and len(args.format) > 0:
        inputs.fileFormat = args.format[0]
        if inputs.fileFormat not in ['csv', 'json','database']:
            parser.error(f"Invalid file format \"{inputs.fileFormat}\". Expecting csv, json or database")
    if args.mode != None and len(args.mode) > 0:
        if args.mode[0] != "performance" and args.mode[0] != "robust":
            parser.error(f"Invalid mode \"{args.mode[0]}\". Available modes -m performance or robust")
        else:
            inputs.mode = args.mode[0]
    if args.nnType != None and len(args.nnType) > 0:
        if args.nnType[0] != "cmpnn" and args.nnType[0] != "dmpnn":
            parser.error(f"Invalid model neural network \"{args.nnType[0]}\". Available neural network types for --nnType: cmpnn, dmpnn")
        else:
            inputs.nnType = args.nnType[0]
    if args.bulk != None and len(args.bulk) > 0:
        if args.bulk[0].isdigit():
            inputs.bulk = int(args.bulk[0])
            if inputs.bulk < 100:
                parser.error("the 'bulk' needs to be an integer and minimum is 100")
        else:
            parser.error("the 'bulk' needs to be an integer and minimum is 100")

    # return (names, smiles, output_file, error_file, output_format, args.nnType[0])
    # print(inputs.__dict__)
    return inputs

# handle one compound at a time
def multiPass(inputs, channels):
    #prepare data
    #- read toxidrome models that are used in FP similarity calculation in SEA
    dir_path = os.path.dirname(os.path.realpath(__file__))
    toxidromeMols, toxidromeFPs = loadToxidromeData_forSEA(dir_path + '/toxidrome_mols.xlsx')

    datapassList = []
    total = len(inputs.smiles)
    for i in range(total):
        channels.writeStateMsg("Working on #"+str(i+1)+" out of "+str(total))
        smiles_original = inputs.smiles[i]
        name = inputs.names[i]
        datapass = PredictionPass(name, smiles_original)
        datapassList.append(datapass)
        try:
            datapass.smiles = standardize_one_smiles(smiles_original)
            if datapass.smiles != None:
                datapass.fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(datapass.smiles), 2, nBits=2048)
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
        if datapass.fp == None:
            msg = f"Failed to generate fingerprints for smiles: {smiles_original}"
            print(msg, file=channels.errorChannel)
            continue

        try:
            run_sea_one_compound(datapass,toxidromeMols, toxidromeFPs)
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            msg = f"Failed to run SEA predictions on smiles: {smiles_original}"
            print(msg, file=channels.errorChannel)

        try:
            run_main_models_multiple_compounds(datapassList, inputs.nnType)
            stan_dev_to_error_one_compound(datapass)
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            msg = f"Failed to run {inputs.nnType} predictions on smiles: {smiles_original}"
            print(msg, file=channels.errorChannel)

    return datapassList

#One stage at a time
def multiStages(inputs, channels):
    #prepare data
    #- read toxidrome models that are used in FP similarity calculation in SEA
    dir_path = os.path.dirname(os.path.realpath(__file__))
    toxidromeMols, toxidromeFPs = loadToxidromeData_forSEA(dir_path + '/toxidrome_mols.xlsx')

    datapassList = []
    total = len(inputs.smiles)

    #Stage 1: process compounds
    channels.writeStateMsg("Processing compounds (standardization, FP calculation) ... ")
    comp_success = 0
    for i in range(total):
        smiles_original = inputs.smiles[i]
        name = inputs.names[i]
        datapass = PredictionPass(name, smiles_original)
        datapassList.append(datapass)
        try:
            datapass.smiles = standardize_one_smiles(smiles_original)
            if datapass.smiles == None:
                datapass.smiles = datapass.original
            datapass.fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(datapass.smiles), 2, nBits=2048)
            comp_success += 1
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            msg = f"RDKit failed to process smiles: {smiles_original}"
            print(msg, file=channels.errorChannel)
    channels.writeStateMsg(str(comp_success) + " out of "+ str(total) +" processed successfully. ")

    #Stage 2: SEA predictions
    sea_pred_success = 0
    channels.writeStateMsg("Running SEA predictions ... ")
    for i in range(len(datapassList)):
        datapass = datapassList[i]
        try:
            run_sea_one_compound(datapass,toxidromeMols, toxidromeFPs)
            sea_pred_success += 1
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            msg = f"Failed to run SEA predictions on smiles: {datapass.original}"
            print(msg, file=channels.errorChannel)
    channels.writeStateMsg("SEA predictions is completed.")        

    #Stage 3: Main predictions
    if inputs.bulk == None:
        smilesList_filtered = []
        for datapass in datapassList:
            try:
                isOrganic = isOrganicAndHasWeight(datapass.smiles, ORGANIC_COMPOUND_MIN_WEIGHT)
            except Exception as e:
                print(e, file=sys.stderr)
                print(traceback.format_exc(), file=sys.stderr)
                isOrganic = False

            if isOrganic :
                smilesList_filtered.append(datapass)
            else:
                msg = "Filtered inorganic compound or weight less than "+str(ORGANIC_COMPOUND_MIN_WEIGHT)+": "+datapass.smiles
                print(msg, file=channels.errorChannel)

        channels.writeStateMsg("Running neural network predictions for all compounds ... ")
        try:
            run_main_models_multiple_compounds(smilesList_filtered, inputs.nnType)
            for i in range(len(smilesList_filtered)):
                stan_dev_to_error_one_compound( smilesList_filtered[i] )
        except Exception as e:
            print(e, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            msg = f"Failed to run neural network predictions for all compounds ..."
            print(msg, file=channels.errorChannel)
            print(e, file=channels.errorChannel)
    else:
        channels.writeStateMsg("Running neural network predictions bulk by bulk ")
        for bulkIndex in range(0, len(datapassList), inputs.bulk):
            startIndex = bulkIndex
            endIndex = min(startIndex + inputs.bulk, len(datapassList))
            #print(startIndex)
            #print(endIndex)

            subDatapassList = datapassList[startIndex:endIndex]

            smilesList_filtered = []
            for datapass in subDatapassList:
                try:
                    isOrganic = isOrganicAndHasWeight(datapass.smiles, ORGANIC_COMPOUND_MIN_WEIGHT)
                except Exception as e:
                    print(e, file=sys.stderr)
                    print(traceback.format_exc(), file=sys.stderr)
                    isOrganic = False

                if isOrganic :
                    smilesList_filtered.append(datapass)
                else:
                    msg = "Filtered inorganic compound or weight less than "+str(ORGANIC_COMPOUND_MIN_WEIGHT)+": "+datapass.smiles
                    print(msg, file=channels.errorChannel)

            channels.writeStateMsg("Running neural network predictions for compounds #" + str(startIndex)+" to #"+str(endIndex))
            
            try:
                run_main_models_multiple_compounds(smilesList_filtered, inputs.nnType)
                for i in range(len(smilesList_filtered)):
                    stan_dev_to_error_one_compound(smilesList_filtered[i])
            except Exception as e:
                print(e, file=sys.stderr)
                print(traceback.format_exc(), file=sys.stderr)
                msg = f"Failed to run neural network predictions"
                print(msg, file=channels.errorChannel)
                print(e, file=channels.errorChannel)

    return datapassList

# assume inputs contain all the needed data including smiles and names 
def main(inputs):
    # the sub module like dmpnn outputs a lot of stuffs
    # dont wanna suppress it in its entirety. yet some clean updates for
    # both states and errors are needed.
    channels = Channels()
    if inputs.outputPath != None:
        channels.outputChannel = open(inputs.outputPath, 'w')
    if inputs.errorPath != None:
        channels.errorChannel = open(inputs.errorPath, 'w')
    if inputs.statePath != None:
        channels.stateChannel = open(inputs.statePath, 'w')

    # Execution based on the mode
    if inputs.mode == 'performance':
        datapassList = multiStages(inputs, channels)
    elif inputs.mode == 'robust':
        datapassList = multiPass(inputs, channels)

    # write the data in the correct format into the output channel
    if inputs.fileFormat == 'csv':
        write_csv(datapassList, channels.outputChannel, False)
    elif inputs.fileFormat == 'json':
        write_json(datapassList, channels.outputChannel)
    else: # format is database
        write_database(datapassList, channels.outputChannel)

    channels.writeStateMsg("Execution finished")

    # close the channel if it was previous opened
    if inputs.errorPath != None:
        channels.errorChannel.close()
    if inputs.statePath != None:
        channels.stateChannel.close()
    if inputs.outputPath != None:
        channels.outputChannel.close()

if __name__ == "__main__":
    inputs = parse_args()
    main(inputs)
