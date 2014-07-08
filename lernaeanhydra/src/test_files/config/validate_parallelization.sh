#!/bin/bash


#to validate a parallelization, one should:
# 1) run the test case in serial and save the data
# 2) run the test case in parallel
# 3) compare the files obtained


#give the input file
input='inputs.txt'


#configure the source code for serial case
./config.py -i inputs.txt -c

#run the serial case

#configure the source code for the parallel case
./config.py -i inputs.txt -c

#run the parallel case
