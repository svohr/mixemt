#! /bin/bash

PYTHONPATH=`pwd`:${PYTHONPATH}

pylint mixemt phylotree.py preprocess.py em.py assemble.py stats.py

python -m unittest discover -s test -p "*_test.py"
