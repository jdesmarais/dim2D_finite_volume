#!/bin/bash

PATH_OUTPUT='inputs.txt' #default path for the output
PATH_INPUT='inputs.txt'  #default path for the input
PARAM_NAME=''            #default name of the parameter changed
PARAM_VALUE=''           #default value of the parameter changed
QUOTE=false

#manage the options
#--------------------------------------------
while getopts i:o:p:v:qh option
  do
  case $option in

      #get the path of the input file
      #submitted by the user
      i)
	  if test -e $OPTARG
	  then
	      PATH_INPUT=$OPTARG
	  fi
	  ;;
      
      #get the path of the output file
      #submitted by the user
      o)
	  PATH_OUTPUT=$OPTARG
	  ;;

      #get the name of the parameter changed
      #by the user
      p)
	  PARAM_NAME=$OPTARG
	  ;;

      #get the value of the parameter changed
      #by the user
      v)
	  PARAM_VALUE=$OPTARG
	  ;;

      q)
	  QUOTE=true
	  ;;

      #display help
      h)
          echo "invalid option"
	  
	  echo "-i: input file path"
	  echo "-o: output file path"
	  echo "-p: name of the parameter changed"
	  echo "-v: value of the parameter changed"
	  echo "-q: add quotes around the parameter changed"
	  ;;

      #default
      \?)
          echo "invalid option"
	  
	  echo "-i: input file path"
	  echo "-o: output file path"
	  echo "-p: name of the parameter changed"
	  echo "-v: value of the parameter changed"
	  echo "-q: add quotes around the parameter changed"
      
  esac
done


#change the parameter by the value given by the user in the input file
#---------------------------------------------------------------------
if [ "$PATH_INPUT" == "$PATH_OUTPUT" ]
then
    TEMP_FILE=true
    PATH_OUTPUT='temp.txt'
else
    TEMP_FILE=false
fi

if $QUOTE
then
    sed -e s/$PARAM_NAME[' ']*'='[' ']*[\']*[-0-9a-zA-Z_\\.e-]*[\']*/$PARAM_NAME' = '\'$PARAM_VALUE\'/g \
    < $PATH_INPUT > $PATH_OUTPUT
else
    sed -e s/$PARAM_NAME[' ']*'='[' ']*[\']*[-0-9a-zA-Z_\\.e-]*[\']*/$PARAM_NAME' = '$PARAM_VALUE/g \
    < $PATH_INPUT > $PATH_OUTPUT
fi

if $TEMP_FILE
then
    mv 'temp.txt' $PATH_INPUT
fi
