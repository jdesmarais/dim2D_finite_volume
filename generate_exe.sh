#!/bin/bash

while getopts i: option
  do
  case $option in

      #input file
      i)
	  INPUT_PATH=$OPTARG

	  #check whether the input file exists
	  if [ -f $INPUT_PATH ];
	  then

	      #generate the executable
	      cd $AUGEANSTABLES_CONFIG
	      ./config.py -i $INPUT_PATH -c
	      cd ..

	      #copy the serial executable
	      #to the main folder
	      if [ -f sim_dim2d ];
	      then
		  mv sim_dim2d $augeanstables
	      fi

	      #copy the serial executable
	      #to the main folder
	      if [ -f sim_dim2d_par ];
	      then
		  mv sim_dim2d_par $augeanstables
	      fi
		  
	  else
	      echo "file does not exist"
	  fi
	  ;;

      #default
      \?)
	  echo "generate the executable corresponding to the input file"
	  echo "./generate_exe.sh -i <input_file>"
	  ;;
  esac
done