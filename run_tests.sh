#! /bin/bash

PYTHONPATH=`$pwd`:${PYTHONPATH}

pylint mixemt phylotree.py preprocess.py em.py assemble.py

python test/test_phylotree.py

