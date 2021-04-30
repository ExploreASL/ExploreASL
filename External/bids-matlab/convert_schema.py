"""
This script runs through the bids-schema YAML files of the BIDS-specification
and converts them to JSON.

Created by Remi Gau
"""

from ruamel.yaml import YAML
import json
import glob
import os

input_dir = "bids-specification/src/schema"
output_dir = "schema"

print("\n\nCONVERTING SCHEMA\n\n")

# create output directories
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

yaml = YAML(typ="safe")

# list all yaml files in an iterator
file_ls = glob.glob(os.path.join(input_dir, "**", "*.yaml"), recursive=True)

for in_file in file_ls:

    print(in_file)

    # create output directory and filename
    path, fname = os.path.split(in_file)

    path = path.replace(input_dir, output_dir)
    if not os.path.exists(path):
        os.makedirs(path)

    out_file = os.path.join(path, fname.replace("yaml", "json"))

    # convert to json
    with open(in_file) as fpi:
        data = yaml.load(fpi)

    with open(out_file, "w") as fpo:
        json.dump(data, fpo, indent=2)
