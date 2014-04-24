#!/bin/bash

#default values
TEST_NEIGHBOR_CASE='1'
TEST_SIZE_CASE='1'
TEST_DISTANCE_CASE='1'
TEST_RANDOM_SEED='86456'
TEST_OVER_ALIGNMENT_CASE_X='1'
TEST_OVER_ALIGNMENT_CASE_Y='1'
TEST_INVERSE_CASE='1'
TEST_INVERSE_SIZE_CASE='1'

#manage the options
while getopts n:s:d:r:x:y:i:j:h option
  do
  case $option in

      #neighbors
      n)
	  TEST_NEIGHBOR_CASE=$OPTARG
	  ;;

      #size
      s)
	  TEST_SIZE_CASE=$OPTARG
	  ;;

      #distance
      d)
	  TEST_DISTANCE_CASE=$OPTARG
	  ;;

      #random_seed
      r)
	  TEST_RANDOM_SEED=$OPTARG
	  ;;

      #over alignment case x
      x)
	  TEST_OVER_ALIGNMENT_CASE_X=$OPTARG
	  ;;

      #over alignment case y
      y)
	  TEST_OVER_ALIGNMENT_CASE_Y=$OPTARG
	  ;;

      #inverse case
      i)
	  TEST_INVERSE_CASE=$OPTARG
	  ;;

      #inverse size case
      j)
	  TEST_INVERSE_SIZE_CASE=$OPTARG
	  ;;

      #help
      h)
	  echo "run the test_bf_layer.sh"
	  echo "-n : neighbors test case"
	  echo "-s : size test case"
	  echo "-d : distance test case"
	  echo "-r : random seed"
	  echo "-x : over alignment for the merge along x-direction"
	  echo "-y : over alignment for the merge along y-direction"
	  echo "-i : inverse test for the two buffer layers "
	  echo "-j : decide which buffer is larger than the other "
	  echo "-h : help"
	  ;;

      #default
      \?)
	  echo "run the test_bf_layer.sh"
	  echo "-n : neighbors test case"
	  echo "-s : size test case"
	  echo "-d : distance test case"
	  echo "-r : random seed"
	  echo "-x : over alignment for the merge along x-direction"
	  echo "-y : over alignment for the merge along y-direction"
	  echo "-i : inverse test for the two buffer layers "
	  echo "-j : decide which buffer is larger than the other "
	  echo "-h : help"
	  ;;
  esac
done


#clean directory from previous compilation
#and data files
make cleanall 1>/dev/null 2>&1
rm *.dat 1>/dev/null 2>&1


#change the test case tested
PATH_INPUT='test_bf_layer_prog.f'
PATH_OUTPUT='temp.f'

sed -e s/' neighbor_case '[' ']*'='[' ']*[-0-9]*/' neighbor_case = '$TEST_NEIGHBOR_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' size_case '[' ']*'='[' ']*[-0-9]*/' size_case = '$TEST_SIZE_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' distance_case '[' ']*'='[' ']*[-0-9]*/' distance_case = '$TEST_DISTANCE_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' random_seed '[' ']*'='[' ']*[-0-9]*/' random_seed = '$TEST_RANDOM_SEED/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' over_alignment_case_x '[' ']*'='[' ']*[-0-9]*/' over_alignment_case_x = '$TEST_OVER_ALIGNMENT_CASE_X/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' over_alignment_case_y '[' ']*'='[' ']*[-0-9]*/' over_alignment_case_y = '$TEST_OVER_ALIGNMENT_CASE_Y/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' inverse_case '[' ']*'='[' ']*[-0-9]*/' inverse_case = '$TEST_INVERSE_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' inverse_size_case '[' ']*'='[' ']*[-0-9]*/' inverse_size_case = '$TEST_INVERSE_SIZE_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT


#clean directory from previous compilation
#and data files
make cleanall 1> /dev/null 2>&1
rm *.dat 1> /dev/null 2>&1

#compile the test case for buffer layer
make test_bf_layer_prog

#run the test case
./test_bf_layer_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_layer.py -f $folder_path
