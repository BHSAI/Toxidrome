from os import listdir
import os

from cmpnn_toxidrome.chemprop_cmpnn.train import make_predictions

from cmpnn_toxidrome.chemprop_cmpnn.parsing import parse_train_args
from cmpnn_toxidrome.chemprop_cmpnn.parsing import modify_train_args

def run_cmpnn(smiles):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    toxidromes = listdir(dir_path+"/toxidromes")

    all_preds = {}
    all_acc = {}

    for toxidrome in toxidromes:
        args = parse_train_args()
        args.checkpoint_dir = dir_path+'/toxidromes/'+toxidrome
        modify_train_args(args)
        preds, smiles, acc = make_predictions(args, smiles)
        
        all_preds[toxidrome] = preds
        all_acc[toxidrome] = acc
    return all_preds, all_acc