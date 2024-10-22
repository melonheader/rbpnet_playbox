#!/usr/bin/env python
import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from rbpnet_playbox.common import process_fasta
from rbpnet import models


"""
A script to run RBPnet.predict for selected RBPs over sequences in the provided fasta file.
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Run RBPnet.predict for selected RBPs over sequences in the provided fasta file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input fasta file")
    parser.add_argument("-m", "--model-dir", type=str, required=True, help="Path to the directories with RBPnet models")
    parser.add_argument("-r", "--rbps", type=str, required=True, help="Names of the RBPs from the collection specified in -m")
    parser.add_argument("-o", "--output-dir", type=str, default=".", help="Path to the directory to write predictions (one file per RBP)")
    parser.add_argument("-n", "--n--cores", type=int, default=1, help="Numeber of threads to parallelise over.")
    args = parser.parse_args()
    return args

def predict_binding(
        path_to_input: str,
        models_dir: str,
        rbps: dict,
        output_dir: str = '.',
        num_workers: int = 1
    ):
    """Predict binding in parallel or sequentially based on user input.
    Args:
        path_to_wmasks (str): Path to the fasta file with window masks.
        rbps (dict): Dictionary of RNA binding proteins (RBPs).
        output_dir (str, optional): Directory to save the output. Defaults to '.'.
        models_dir (str, optional): Directory where models are located. Defaults to 'gasoline/resources/models_rbpnet'.
        num_workers (int, optional): Number of workers to use for parallel execution. Defaults to 1.
    """
    seqs = process_fasta(path_to_input)
    def process_rbp(rbp):
        """Helper function to process RBP."""
        model_file = f'{rbp}_HepG2.model.h5'
        model_path = os.path.join(models_dir, model_file)
        model = models.load_model(model_path)
        output_path = os.path.join(output_dir, f'{rbp}_preds.csv')
        with open(output_path, 'w') as preds:
            for seqname, seq in seqs.items():             
                pred = model.predict_from_sequence(seq)[f'{rbp}_HepG2_profile_target']
                pred_num = [t.numpy() for t in pred]
                row = [seqname] + list(pred_num)
                preds.write(','.join(map(str, row)) + '\n')
    if num_workers > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_rbp, rbp): rbp for rbp in rbps}
            for future in as_completed(futures):
                rbp = futures[future]
                try:
                    future.result() 
                    print(f'RBP {rbp} processed successfully.')
                except Exception as e:
                    print(f'Error processing RBP {rbp}: {e}')
    else:
        for rbp in rbps:
            process_rbp(rbp)
            print(f'RBP {rbp} processed successfully.')

if __name__ == "__main__":
    args = parse_args()
    predict_binding(
        path_to_input = args.input,
        models_dir = args.model_dir,
        rbps = args.rbps,
        output_dir = args.output_dir,
        num_workers = args.n_cores
    )



