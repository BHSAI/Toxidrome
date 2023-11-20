from os import listdir
import chemprop
import math

from chemprop.train import make_predictions

def run_dmpnn(smiles):
    smiles = [[s] for s in smiles]
    toxidromes = listdir(f"C:/Users/bclancy/Desktop/projects/toxidrome_clt/dmpnn_toxidrome/toxidromes")

    all_preds = {}
    all_acc = {}

    for toxidrome in toxidromes:
        arguments = [
            '--test_path', '/dev/null',
            '--preds_path', '/dev/null',
            '--checkpoint_dir', f'C:/Users/bclancy/Desktop/projects/toxidrome_clt/dmpnn_toxidrome/toxidromes/{toxidrome}',
            '--uncertainty_method', 'ensemble'
        ]
        args = chemprop.args.PredictArgs().parse_args(arguments)
        #preds, smiles, acc = make_predictions(args, smiles)
        preds, acc = make_predictions(args, smiles, return_uncertainty = True)

        for i in range(len(acc)):
            for j in range(len(acc[i])):
                acc[i][j] = math.sqrt(acc[i][j])

        all_preds[toxidrome] = preds
        all_acc[toxidrome] = acc
    return all_preds, all_acc