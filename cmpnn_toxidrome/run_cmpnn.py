from os import listdir

from cmpnn_toxidrome.chemprop_cmpnn.train import make_predictions

from cmpnn_toxidrome.chemprop_cmpnn.parsing import parse_train_args
from cmpnn_toxidrome.chemprop_cmpnn.parsing import modify_train_args

def run_cmpnn(smiles):
    toxidromes = listdir(f"C:/Users/bclancy/Desktop/projects/toxidrome_clt/cmpnn_toxidrome/toxidromes")

    all_preds = {}
    all_acc = {}

    for toxidrome in toxidromes:
        args = parse_train_args()
        args.checkpoint_dir = f'C:/Users/bclancy/Desktop/projects/toxidrome_clt/cmpnn_toxidrome/toxidromes/{toxidrome}'
        modify_train_args(args)
        preds, smiles, acc = make_predictions(args, smiles)
        
        all_preds[toxidrome] = preds
        all_acc[toxidrome] = acc
    return all_preds, all_acc