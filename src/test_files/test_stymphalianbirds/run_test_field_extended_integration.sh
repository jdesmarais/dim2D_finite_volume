#!/bin/bash

TEMP=temp.f

#clean directory from previous compilation
#and data files
make cleanall
rm *.dat


#---------------------------------------
#compile the test case for the interface
#with one piece domain
#---------------------------------------
#domain_decomposition -> one_piece
INPUT=test_field_extended_integration.f
change_parameter.sh -i $INPUT -o $TEMP -p domain_decomposition -v one_piece
mv $TEMP $INPUT

#nx -> 40
#ny -> 40
INPUT=$augeanstables/src/parameters/parameters_input.f
change_parameter.sh -i $INPUT -o $TEMP -p nx -v 40
mv $TEMP $INPUT

change_parameter.sh -i $INPUT -o $TEMP -p ny -v 40
mv $TEMP $INPUT

#compile and run one-piece test case
make test_field_extended_integration
./test_field_extended_integration


#---------------------------------------
#compile the test case for the interface
#with four pieces domain
#---------------------------------------
#domain_decomposition -> four_pieces
INPUT=test_field_extended_integration.f
change_parameter.sh -i $INPUT -o $TEMP -p domain_decomposition -v four_pieces
mv $TEMP $INPUT

#nx -> 20
#ny -> 20
INPUT=$augeanstables/src/parameters/parameters_input.f
change_parameter.sh -i $INPUT -o $TEMP -p nx -v 20
mv $TEMP $INPUT

change_parameter.sh -i $INPUT -o $TEMP -p ny -v 20
mv $TEMP $INPUT

#compile and run four-pieces test case
make test_field_extended_integration
./test_field_extended_integration


#---------------------------------------
#run the visualization using python
#---------------------------------------
folder_path=$(pwd)
pushd $augeanstables/scripts_py/scripts_erymanthianboar
python test_field_extended_integration.py -f $folder_path

##clean original directory
#pushd
#echo ''
#echo 'Cleaning directory...'
#make cleanall 1>/dev/null 2>&1
#rm *.dat
#rm *.nc