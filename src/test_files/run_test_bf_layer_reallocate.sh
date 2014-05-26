#!/bin/bash

#default values
TEST_RELATIVE_SIZE='1'
TEST_RELATIVE_DISTANCE='1'
TEST_FINAL_ALIGNMENT='1'
TEST_OVER_ALLOCATED='0'

#manage the options
while getopts s:d:a:o:h option
  do
  case $option in

      #size
      s)
	  TEST_RELATIVE_SIZE=$OPTARG
	  ;;

      #distance
      d)
	  TEST_RELATIVE_DISTANCE=$OPTARG
	  ;;

      #final alignmnent
      a)
	  TEST_FINAL_ALIGNMENT=$OPTARG
	  ;;

      #over allocated
      o)
	  TEST_OVER_ALLOCATED=$OPTARG
	  ;;

      #help
      h)
	  echo "run the test_bf_layer.sh"
	  echo "-s : size test case"
	  echo "-d : distance test case"
	  echo "-f : final alignment test case"
	  echo "-o : over allocated test case"
	  echo "-h : help"
	  ;;

      #default
      \?)
	  echo "run the test_bf_layer.sh"
	  echo "-s : size test case"
	  echo "-d : distance test case"
	  echo "-f : final alignment test case"
	  echo "-o : over allocated test case"
	  echo "-h : help"
	  ;;
  esac
done


#clean directory from previous compilation
#and data files
make cleanall 1>/dev/null 2>&1
rm *.dat 1>/dev/null 2>&1


#change the test case tested
PATH_INPUT='test_bf_layer_reallocate_prog.f'
PATH_OUTPUT='temp.f'

sed -e s/' test_relative_size '[' ']*'='[' ']*[-0-9]*/' test_relative_size = '$TEST_RELATIVE_SIZE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' test_relative_distance '[' ']*'='[' ']*[-0-9]*/' test_relative_distance = '$TEST_RELATIVE_DISTANCE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' test_final_alignment '[' ']*'='[' ']*[-0-9]*/' test_final_alignment = '$TEST_FINAL_ALIGNMENT/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' over_allocated '[' ']*'='[' ']*[-0-9]*/' over_allocated = '$TEST_OVER_ALLOCATED/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT


#clean directory from previous compilation
#and data files
make cleanall 1> /dev/null 2>&1
rm *.dat 1> /dev/null 2>&1

#compile the test case for buffer layer
make test_bf_layer_reallocate_prog

#run the test case
./test_bf_layer_reallocate_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_layer_reallocate.py -f $folder_path
