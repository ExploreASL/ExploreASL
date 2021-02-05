#xASL_ut_generate_DRO Script to generate a DRO test dataset
#
# INPUT:        testdir - xASL testing directory
#
# OUTPUT:       Store default ASL DRO input parameters and the default test dataset
#
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DESCRIPTION:  This is a script to generate a DRO test dataset for QC and testing in general.
# 				Install the ASL DRO python package before running this script.
# 				https://pypi.org/project/asldro/
# 				pip install asldro
#
# EXAMPLE:      python xASL_ut_generate_DRO.py --testdir path\to\testdir
# 				(insert the path to the test dir without typewriter apostrophes)
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright 2015-2021 ExploreASL

import os
import asldro
import zipfile
import argparse

# Initialize argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--testdir', help='Path to the xASL testing directory')

# Parse arguments
args = parser.parse_args()

# Print testing directory
print('Testing directory:\t' + str(args.testdir))

# Create output directory
if not os.path.exists(os.path.join(str(args.testdir),"DRO")):
		os.makedirs(os.path.join(str(args.testdir),"DRO"))

# Create default ASL DRO parameter file
inputParmsFile = os.path.join(str(args.testdir),"DRO","input_params.json")
os.system("asldro output params "+inputParmsFile)

# Run default ASL DRO workflow
testDRO = os.path.join(str(args.testdir),"DRO","test_dataset.zip")
os.system("asldro generate "+testDRO)

# Unzip folder
testPatient = os.path.join(str(args.testdir),"DRO","test_patient")
with zipfile.ZipFile(testDRO, 'r') as zip_ref:
    zip_ref.extractall(testPatient)

# Delete zip file
os.remove(testDRO) 


