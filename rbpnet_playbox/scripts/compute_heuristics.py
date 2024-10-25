#!/usr/bin/env python
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import kurtosis

"""
A script to compute heuristics from RBPnet predictions.
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Compute per-sequence heuristics from RBPnet predictions.")
    parser.add_argument("-i", "--input", required=True, help="Path to the directory with RBPnet predictions")
    parser.add_argument("-o", "--output-dir", type=str, default=None, help="Path to the directory to write heuristics (defaults to input directory)")
    args = parser.parse_args()
    return args
def compute_metrics(preds, preds_id):
        """Computes and returns relevant metrics for RBPnet scores."""
        # subtract an average to remove noise
        preds_centered = preds - np.mean(preds)
        preds_centered[preds_centered < 0] = 0
        # heuristics to compute can be adjusted at one's own volition
        metrics = {
            "id": preds_id,
            "seq_length": len(preds_centered),
            "mean": np.mean(preds_centered),
            "variance": np.var(preds_centered),
            "kurtosis": kurtosis(preds_centered)
        }
        return metrics
def compute_heuristics(path_to_preds_dir: str, output_dir: str = None):
    if not output_dir:
        output_dir = path_to_preds_dir
    results = []
    for model_out in os.listdir(path_to_preds_dir):
        if '_preds.csv' in model_out:
            with open(os.path.join(path_to_preds_dir, model_out), 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith("id"):
                        continue
                    preds_values = np.array([float(x) for x in line.split(',')[1:]])
                    id = line.split(',', 1)[0]
                    result = compute_metrics(preds_values, id)
                    results.append(result)
            rbp = model_out.replace('_preds.csv', '')
            path_to_write = os.path.join(output_dir, f'{rbp}_heuristics.csv')
            df = pd.DataFrame(results)
            results = []  
            df.to_csv(path_to_write, index=False)
            print(f'Summarised {model_out} predictions to {rbp}_heuristics.csv')
            
if __name__ == "__main__":
    args = parse_args()
    compute_heuristics(
        path_to_preds_dir = args.input,
        output_dir = args.output_dir
    )