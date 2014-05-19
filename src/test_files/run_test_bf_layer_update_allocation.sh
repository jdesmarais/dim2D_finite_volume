#!/bin/bash

#default values
TEST_CASE_ID='11'
TEST_MAINLAYER_ID='1'
TEST_RANDOM_SEED='86456'

#manage the options
while getopts t:m:r:h option
  do
  case $option in

      #test case id
      t)
	  TEST_CASE_ID=$OPTARG
	  ;;

      #mainlayer_id
      m)
	  TEST_MAINLAYER_ID=$OPTARG
	  ;;

      #random seed
      r)
	  TEST_RANDOM_SEED=$OPTARG
	  ;;

      #help
      h)
	  echo "run the test_interface_merge.sh"
	  echo "-t: set test case ID"
	  echo "-m: set mainlayer id"
	  echo "-r: set random seed"
	  ;;

      #default
      \?)
	  echo "run the test_interface_merge.sh"
	  echo "-t: set test case ID"
	  echo "-m: set mainlayer id"
	  echo "-r: set random seed"
  esac
done


#clean directory from previous compilation
#and data files
make cleanall 1>/dev/null 2>&1
rm *.dat 1>/dev/null 2>&1

#change the test case tested

PATH_INPUT='test_bf_layer_update_allocation_prog.f'
PATH_OUTPUT='temp.f'

sed -e s/' test_case_id '[' ']*'='[' ']*[-0-9]*/' test_case_id = '$TEST_CASE_ID/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' mainlayer_id '[' ']*'='[' ']*[-0-9]*/' mainlayer_id = '$TEST_MAINLAYER_ID/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' random_seed '[' ']*'='[' ']*[-0-9]*/' random_seed = '$TEST_RANDOM_SEED/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

#compile the test case for the interface
make test_bf_layer_update_allocation_prog

#run the test case
./test_bf_layer_update_allocation_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_layer_update_allocation.py -f $folder_path
