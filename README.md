# rbpnet_playbox

## Foreword
This repository contains a small toolkit designed to comfortably run RNPnet.predict_from_sequence for a desired set of RBPs for every sequence in the input fasta file. The predicitons can be parallelised. 

## Usage

### Setup
One would need to setup the environment and install pysster. Due to old dependencies, it is recommended to setup a separate environemt for pysster to avoid conflict with newer python packages:
```bash
git clone https://github.com/melonheader/rbpnet_playbox.git
cd rbpnet_playbox
conda env create -f rbpnet_playbox/auxiliary/rbpnet.yaml
conda activate rbpnet
```
### Main tool
The package contains two scripts to:
1) Run predictions over sequences in a fasta file
2) Compute heuristics for the output predictions
```bash
python -m rbpnet_playbox.scripts.rbpnet_predict --help   
usage: rbpnet_predict.py [-h] -i INPUT -m MODEL_DIR -r RBPS [-o OUTPUT_DIR] [-n N__CORES]

Run RBPnet.predict for selected RBPs over sequences in the provided fasta file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input fasta file
  -m MODEL_DIR, --model-dir MODEL_DIR
                        Path to the directories with RBPnet models
  -r RBPS, --rbps RBPS  Names of the RBPs from the collection specified in -m
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Path to the directory to write predictions (one file per RBP)
  -n N__CORES, --n--cores N__CORES
                        Numeber of threads to parallelise over.

python -m rbpnet_playbox.scripts.compute_heuristics --help
usage: compute_heuristics.py [-h] -i INPUT [-o OUTPUT_DIR]

Compute per-sequence heuristics from RBPnet predictions.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the directory with RBPnet predictions
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Path to the directory to write heuristics (defaults to input directory)
```