#!/bin/bash

#default values
TEST_NO_ALIGNMENT='1'
TEST_CASE_ID='3'
TEST_NO_ALIGNMENT='1'
TEST_CASE_ID='3'
TEST_MAINLAYER_ID='1'
TEST_NB_SUBLAYERS='2'
TEST_SIZE_CASE='1'
TEST_DISTANCE_CASE='1'
TEST_RANDOM_SEED='86456'
TEST_MERGE_ID='1'
TEST_MERGE_INVERSE='0'
TEST_NEIGHBOR_CASE='1'

#manage the options
while getopts a:t:m:l:s:d:r:i:v:n:h option
  do
  case $option in

      #no_alignement
      a)
	  TEST_NO_ALIGNMENT=$OPTARG
	  ;;

      #test case id
      t)
	  TEST_CASE_ID=$OPTARG
	  ;;

      #mainlayer_id
      m)
	  TEST_MAINLAYER_ID=$OPTARG
	  ;;

      #nb_sublayers
      l)
	  TEST_NB_SUBLAYERS=$OPTARG
	  ;;

      #size_case study
      s)
	  TEST_SIZE_CASE=$OPTARG
	  ;;

      #distance_case study
      d)
	  TEST_DISTANCE_CASE=$OPTARG
	  ;;

      #random seed
      r)
	  TEST_RANDOM_SEED=$OPTARG
	  ;;

      #merge id
      i)
	  TEST_MERGE_ID=$OPTARG
	  ;;

      #merge inverse
      v)
	  TEST_MERGE_INVERSE=$OPTARG
	  ;;

      #neighbors
      n)
	  TEST_NEIGHBOR_CASE=$OPTARG
	  ;;

      #help
      h)
	  echo "run the test_interface_merge.sh"
	  echo "-a: set no alignement"
	  echo "-t: set test case ID"
	  echo "-m: set mainlayer id"
	  echo "-l: set number of sublayers"
	  echo "-s: set size_case study"
	  echo "-d: set distance_case study"
	  echo "-r: set random seed"
	  echo "-i: set merge_id"
	  echo "-v: set merge_inverse"
	  echo "-n: set neighbor_case"
	  ;;

      #default
      \?)
	  echo "run the test_interface_merge.sh"
	  echo "-a: set no alignement"
	  echo "-t: set test case ID"
	  echo "-m: set mainlayer id"
	  echo "-l: set number of sublayers"
	  echo "-s: set size_case study"
	  echo "-d: set distance_case study"
	  echo "-r: set random seed"
	  echo "-i: set merge_id"
	  echo "-v: set merge_inverse"
	  echo "-n: set neighbor_case"      
  esac
done


#clean directory from previous compilation
#and data files
make cleanall 1>/dev/null 2>&1
rm *.dat 1>/dev/null 2>&1

#change the test case tested

PATH_INPUT='test_bf_sublayers_merge_prog.f'
PATH_OUTPUT='temp.f'

sed -e s/' no_alignment '[' ']*'='[' ']*[-0-9]*/' no_alignment = '$TEST_NO_ALIGNMENT/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' test_case_id '[' ']*'='[' ']*[-0-9]*/' test_case_id = '$TEST_CASE_ID/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' mainlayer_id '[' ']*'='[' ']*[-0-9]*/' mainlayer_id = '$TEST_MAINLAYER_ID/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' nb_sublayers '[' ']*'='[' ']*[-0-9]*/' nb_sublayers = '$TEST_NB_SUBLAYERS/g \
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

sed -e s/' merge_id '[' ']*'='[' ']*[-0-9]*/' merge_id = '$TEST_MERGE_ID/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' merge_inverse '[' ']*'='[' ']*[-0-9]*/' merge_inverse = '$TEST_MERGE_INVERSE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

sed -e s/' neighbor_case '[' ']*'='[' ']*[-0-9]*/' neighbor_case = '$TEST_NEIGHBOR_CASE/g \
    < $PATH_INPUT > $PATH_OUTPUT
mv $PATH_OUTPUT $PATH_INPUT

#compile the test case for the interface
make test_bf_sublayers_merge_prog

#run the test case
./test_bf_sublayers_merge_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_sublayers_merge.py -f $folder_path
