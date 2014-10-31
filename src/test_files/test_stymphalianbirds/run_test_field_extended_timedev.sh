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
INPUT=test_field_extended_timedev.f
change_parameter.sh -i $INPUT -o $TEMP -p domain_decomposition -v one_piece
mv $TEMP $INPUT

#nx -> 40
#ny -> 40
INPUT=$augeanstables/src/parameters/parameters_input.f
change_parameter.sh -i $INPUT -o $TEMP -p ntx -v 40
mv $TEMP $INPUT

change_parameter.sh -i $INPUT -o $TEMP -p nty -v 40
mv $TEMP $INPUT

#compile and run one-piece test case
make test_field_extended_timedev
mv test_field_extended_timedev test_field_extended_timedev_onepiece
./test_field_extended_timedev_onepiece


#---------------------------------------
#compile the test case for the interface
#with four pieces domain
#---------------------------------------
#domain_decomposition -> four_pieces
INPUT=test_field_extended_timedev.f
change_parameter.sh -i $INPUT -o $TEMP -p domain_decomposition -v four_pieces
mv $TEMP $INPUT

#nx -> 24
#ny -> 24
INPUT=$augeanstables/src/parameters/parameters_input.f
change_parameter.sh -i $INPUT -o $TEMP -p ntx -v 24
mv $TEMP $INPUT

change_parameter.sh -i $INPUT -o $TEMP -p nty -v 24
mv $TEMP $INPUT

#compile and run four-pieces test case
make test_field_extended_timedev
mv test_field_extended_timedev test_field_extended_timedev_fourpieces
./test_field_extended_timedev_fourpieces


#---------------------------------------
#run the visualization using python
#---------------------------------------
folder_path=$(pwd)
pushd $augeanstables/scripts_py/scripts_erymanthianboar
python test_field_extended_timedev.py -f $folder_path

##clean original directory
#pushd
#echo ''
#echo 'Cleaning directory...'
#make cleanall 1>/dev/null 2>&1
#rm *.dat
#rm *.nc