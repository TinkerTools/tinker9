#!/bin/bash

# Run a python script for all YAML files in this directory to generate their cuda kernals

for f in *.yaml; do
    python3 ../ck3.py -c $f | bash
done
