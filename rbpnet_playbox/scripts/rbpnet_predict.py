#!/usr/bin/env python
import argparse
from math import exp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from rbpnet_playbox.common import process_fasta
from typing import Dict
from rbpnet import models


"""
A script to run RBPnet.predict for selected RBPs over sequences in the provided fasta file.
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Run RBPnet.predict for selected RBPs over sequences in the provided fasta file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input fasta file")
    parser.add_argument("-m", "--model-dir", type=str, required=True, help="Path to the directories with RBPnet models")
    parser.add_argument("-r", "--rbps", nargs='+', required=True, help="Names of the RBPs from the collection specified in -m")
    parser.add_argument("-o", "--output-dir", type=str, default=".", help="Path to the directory to write predictions (one file per RBP)")
    parser.add_argument("-n", "--n-cores", type=int, default=1, help="Numeber of threads to parallelise over.")
    args = parser.parse_args()
    return args
def process_rbp(
            rbp: str, seqs: Dict[str, str], models_dir, output_dir=''
        ) -> None:
        """
        Helper function to process an individual RBP and predict binding.
        Args:
            rbp (str): Name of the RNA binding protein.
            seqs (dict): Dictionary of sequences with sequence names as keys and sequences as values.
        """
        # Construct the path to the model file for the given RBP
        model_file = f'{rbp}_HepG2.model.h5'
        model_path = os.path.join(models_dir, model_file)
        # Load the RBP model
        model = models.load_model(model_path)
        # Define the output file path for the predictions
        output_path = os.path.join(output_dir, f'{rbp}_preds.csv')
        # Write predictions to the output file
        with open(output_path, 'w') as preds:
            header_written = False
            for seqname, seq in seqs.items():
                # Write the header once
                if not header_written:
                    header = ['id'] + [str(i + 1) for i in range(len(seq))]
                    preds.write(','.join(header) + '\n')
                    header_written = True
                # Make predictions for the sequence
                pred = model.predict_from_sequence(seq)[f'{rbp}_HepG2_profile_target']
                # Convert predictions to a list of numeric values
                pred_num = [t.numpy() for t in pred]
                # Create a row with the sequence name and predictions
                row = [seqname] + list(pred_num)
                # Write the row to the CSV file
                preds.write(','.join(map(str, row)) + '\n')
def predict_binding(
        path_to_input: str,
        models_dir: str,
        rbps: list,
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
    if num_workers > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_rbp, rbp, seqs, models_dir, output_dir): rbp for rbp in rbps}
            for future in as_completed(futures):
                rbp = futures[future]
                try:
                    future.result() 
                    print(f'RBP {rbp} processed successfully.')
                except Exception as e:
                    print(f'Error processing RBP {rbp}: {e}')
    else:
        for rbp in rbps:
            try:
                process_rbp(rbp, seqs, models_dir)
                print(f'RBP {rbp} processed successfully.')
            except Exception as e:
                print(f'Error processing RBP {rbp}: {e}')

if __name__ == "__main__":
    args = parse_args()
    predict_binding(
        path_to_input = args.input,
        models_dir = args.model_dir,
        rbps = args.rbps,
        output_dir = args.output_dir,
        num_workers = args.n_cores
    )