import argparse

parser = argparse.ArgumentParser(description = 'Inputs a protein sequence file for processing by seq2pdb')
parser.add_argument('-s', '--sequence', required=True, help='Input protein amino acid sequence')
parser.add_argument('-o', '--output', default='protein', help='Name of pdb output file')
parser.add_argument('-x', '--additional_args', default='', help='Any additional args for Deepfri must be typed in full here')

args = parser.parse_args()