#!/bin/bash


PATH_INPUT='/home/jdesmarais/Code/nemeanlion/input/inputs_dim1d_interface_pml_cst.txt'


#manage the options
#--------------------------------------------
while getopts i:p: option
  do
  case $option in

      #get the path of the input file
      #submitted by the user
      i)
	  if test -e $OPTARG
	      then
	      PATH_INPUT=$OPTARG
	  else
	      echo "input file $OPTARG does not exist"
	      echo "default input file used"
	  fi
	  ;;

      #get the name of the parameter
      #extracted from the input file
      p)
	  PARAM_NAME=$OPTARG
	  ;;

      #default
      \?)
          echo "invalid option"
	  echo "-i: input file"
      
  esac
done



#extract the line where the parameter is mentionned
#--------------------------------------------------
line=$(grep $PARAM_NAME[' ']*'='[' ']*[A-Za-z0-9\\.]* $PATH_INPUT)


#extract the value out of the line
#---------------------------------
#value=$(echo $line | sed -e  s/[A-Za-z0-9_" "]*"="[" "]*/''/g)
value=$(echo $line | sed -e  s/[^"=""!"]*"="[" "]*/''/g)
value=$(echo $value | sed -e  s/[' ']*'!'[^'!']*/''/g)
echo $value
