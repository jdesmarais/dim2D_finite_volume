#!/bin/bash

doxygen DoxygenConfigFortran #1> /dev/null 2>&1

if [ -d "./doc" ]
then
rm -rf ./doc
fi 
mkdir ./doc
mv ./html ./doc