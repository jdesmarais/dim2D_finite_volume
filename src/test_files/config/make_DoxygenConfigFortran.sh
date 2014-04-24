#!/bin/bash

PATH_INPUT='DoxygenConfigFortran'
PATH_TEMP='temp'

#modify the path for the source code depending
#on the relative path on the cluster
sed -e 's|^INPUT [ ]*=[ ]*[a-zA-Z0-9/\"]*|'"INPUT                  = \"$erymanthianboar/src\""'|g' < $PATH_INPUT > $PATH_TEMP
mv $PATH_TEMP $PATH_INPUT

#modify the path for the documentation pictures
#depending on the relative path on the cluster
sed -e 's|^IMAGE_PATH [ ]*=[ ]*[a-zA-Z0-9/\"]*|'"IMAGE_PATH             = \"$erymanthianboar/src/test_files/config/doc_pictures\""'|g' < $PATH_INPUT > $PATH_TEMP
mv $PATH_TEMP $PATH_INPUT