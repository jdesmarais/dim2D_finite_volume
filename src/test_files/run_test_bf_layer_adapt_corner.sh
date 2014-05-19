#!/bin/bash

#default values
TEST_CORNER='1'
TEST_CORNER_ORDER='1'
TEST_CASE_TESTED='21'
TEST_BF_LAYER_DISTANCE='0'
TEST_OVER_ALLOCATED='0'
TEST_CURRENT_PATH='0'

#manage the options
while getopts c:o:i:b:a:p:h option
  do
  case $option in

      #corner ID
      c)
	  TEST_CORNER=$OPTARG
	  ;;

      #corner order
      o)
	  TEST_CORNER_ORDER=$OPTARG
	  ;;

      #get the ID of the test case
      i)
	  TEST_CASE_TESTED=$OPTARG
	  ;;

      #get the buffer layer distance
      b)
	  TEST_BF_LAYER_DISTANCE=$OPTARG
	  ;;

      #over_allocated parameter
      a)
	  TEST_OVER_ALLOCATED=$OPTARG
	  ;;

      #current path
      p)
	  TEST_CURRENT_PATH=$OPTARG
	  ;;

      #display help
      h)
	  echo "run the test_bf_later_corner_check.sh"
	  echo "-c: ID of the corner tested"
	  echo "-o: corner priority for test case i=1"
	  echo "-i: ID of the test case tested"
	  echo "-b: distance b\w the buffer layer and the corner"
	  echo "-a: over allocated parameter"
	  echo "-p: use of current path"
	  ;;

      #default
      \?)
          echo "invalid option"
	  echo "-c: ID of the corner tested"
	  echo "-o: corner priority for test case i=1"
	  echo "-i: ID of the test case tested"
	  echo "-b: distance b\w the buffer layer and the corner"
	  echo "-a: over allocated parameter"
	  echo "-p: use of current path"
      
  esac
done


#clean directory from previous compilation
#and data files
make cleanall 1> /dev/null 2>&1
rm *.dat 1> /dev/null 2>&1

#change the test case tested

PATH_INPUT='test_bf_layer_adapt_corner_prog.f'
PATH_OUTPUT='temp.f'

sed -e s/'corner_tested '*'='[' ']*[-0-9]*/'corner_tested = '$TEST_CORNER/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/'corner_order '*'='[' ']*[-0-9]*/'corner_order = '$TEST_CORNER_ORDER/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/'test_case_id '*'='[' ']*[-0-9]*/'test_case_id = '$TEST_CASE_TESTED/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/'bf_corner_distance '*'='[' ']*[-0-9]*/'bf_corner_distance = '$TEST_BF_LAYER_DISTANCE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/'over_allocated '*'='[' ']*[-0-9]*/'over_allocated = '$TEST_OVER_ALLOCATED/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/'current_path_use '*'='[' ']*[-0-9]*/'current_path_use = '$TEST_CURRENT_PATH/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT


#compile the test case for the interface
make test_bf_layer_adapt_corner_prog

#run the test case
./test_bf_layer_adapt_corner_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_layer_adapt_corner.py -f $folder_path
