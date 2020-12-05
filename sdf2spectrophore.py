#!/usr/bin/python
# Authors: Hans De Winter
# Authors: Fabio Mendes

# General Python
from __future__ import print_function
import numpy as np

# RDkit stuff
from rdkit import Chem
from rdkit.Chem import AllChem

# Spectrophores
from spectrophores.spectrophore import spectrophoreCalculator

# Download this from http://pypi.python.org/pypi/futures
from concurrent import futures

# Download this from http://pypi.python.org/pypi/progressbar
import progressbar

# Command-line arguments
import argparse


# Function to calculate spectrophores
def calculateSpectrophore(molecule, calculator, label):
	spec = calculator.calculate(molecule)
	if not np.all(spec): return None
	specString = np.array2string(spec, max_line_width=1000, suppress_small=True, formatter={'float':lambda x: "%.5f" % x})
	specString = label + " " + specString[1:-1]
	return specString
 

# Function to check the value of the resolution command-line argument
def __check_resolution(value):
	v = float(value)
	if v <= 0: raise argparse.ArgumentTypeError("should be larger than 0")
	return v


# Function to parse the command-line arguments
def __processCommandline():
    parser = argparse.ArgumentParser(description = "Calculate spectrophores",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--norm', 
                        dest = 'norm', 
                        help = 'normalization setting', 
                        default = 'none', 
                        choices = ['none','mean','all','std'])
    parser.add_argument('-s', '--stereo', 
                        dest = 'stereo', 
                        help = 'stereo setting', 
                        default = 'none', 
                        choices = ['none','unique','mirror','all'])
    parser.add_argument('-a', '--accuracy', 
                        dest = 'accuracy', 
                        help = 'accuracy setting', 
                        type = int, 
                        default = 20, 
                        choices = [1, 2, 5, 10, 15, 20, 30, 36, 45, 60])
    parser.add_argument('-r', '--resolution', 
                        dest = 'resolution', 
                        help = 'resolution setting (>0)', 
                        default = 3, 
                        type = __check_resolution)
    parser.add_argument('-p', '--np',
                        dest = 'max_workers',
                        help = 'number of processors to use; -1 is all processors',
                        default = -1,
                        type = int)
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--in', 
                               dest = 'infile', 
                               help = 'input sdf file', 
                               required = True)
    requiredNamed.add_argument('-o', '--out', 
                               dest = 'outfile', 
                               help = 'output spectrophore file', 
                               required = True)
    return parser.parse_args()


# Main
if __name__ == "__main__":
    args = __processCommandline()

    calculator = spectrophoreCalculator(normalization = args.norm, 
                                        stereo = args.stereo, 
                                        accuracy = args.accuracy, 
                                        resolution = args.resolution)
    print ("Normalization: ", args.norm)
    print ("Stereo:        ", args.stereo)
    print ("Accuracy:      ", args.accuracy)
    print ("Resolution:    ", args.resolution)
    print ("Input file:    ", args.infile)
    print ("Output file:   ", args.outfile)
    print ("Processors:    ", args.max_workers)
    supplier = Chem.SDMolSupplier(args.infile, removeHs=False)
    of = open(args.outfile, 'w')
    mw = None
    if args.max_workers > 0: mw = args.max_workers
    with futures.ProcessPoolExecutor(max_workers = mw) as executor:
	    jobs = []
	    for mol in supplier:
	    	if mol:
		    	label = mol.GetProp("_Name")
		    	job = executor.submit(calculateSpectrophore, mol, calculator, label)
		    	jobs.append(job)
	
	    widgets = ["Generating spectrophores; ", progressbar.Percentage(), " ", progressbar.ETA(), " ", progressbar.Bar()]
	    pbar = progressbar.ProgressBar(widgets = widgets, maxval = len(jobs))
	    for job in pbar(futures.as_completed(jobs)):
		    spec = job.result()
		    if spec is not None: of.write("%s\n" % (spec))

    of.close()

