from os import listdir
import os
import chemprop
import math

from chemprop.train import make_predictions

def run_dmpnn(smiles):
    smiles = [[s] for s in smiles]
    dir_path = os.path.dirname(os.path.realpath(__file__))
    #Assuming the folder toxidromes exists along side the current file
    toxidromes = listdir(dir_path+"/toxidromes")

    all_preds = {}
    all_acc = {}

    for toxidrome in toxidromes:
        arguments = [
            '--test_path', '/dev/null',
            '--preds_path', '/dev/null',
            '--checkpoint_dir', dir_path+'/toxidromes/'+toxidrome,
            '--uncertainty_method', 'ensemble'
        ]
        args = chemprop.args.PredictArgs().parse_args(arguments)
        #preds, smiles, acc = make_predictions(args, smiles)
        preds, acc = make_predictions(args, smiles, return_uncertainty = True)

        for i in range(len(acc)):
            for j in range(len(acc[i])):
                val_acc = acc[i][j]
                if isinstance(val_acc, int) or isinstance(val_acc, float): #ignore complex
                    acc[i][j] = math.sqrt(val_acc)
                else:
                    acc[i][j] = None
                val_pred = preds[i][j]
                if not (isinstance(val_pred, int) or isinstance(val_pred, float)): #ignore complex
                    preds[i][j] = None # in case if 'invalid smiles'

        all_preds[toxidrome] = preds
        all_acc[toxidrome] = acc
    return all_preds, all_acc