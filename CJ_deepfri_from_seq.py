import subprocess
from CJ_argparse import args
import os

pathExist = os.path.exists(f"./results/{args.output}")
if pathExist == 0:
    os.mkdir(f"./results/{args.output}")

seq2pdbchain = f"python ./seq2pdb/seq2pdbchain.py -s {args.sequence} > './results/{args.output}/{args.output}.pdb'"
subprocess.run(seq2pdbchain, shell=True)

pdb2deepfri = f"python ./Deepfri/predict.py -pdb './results/{args.output}/{args.output}.pdb' -ont mf bp cc ec -v -o './results/{args.output}/{args.output}' {args.additional_args}"
subprocess.run(pdb2deepfri, shell=True)
