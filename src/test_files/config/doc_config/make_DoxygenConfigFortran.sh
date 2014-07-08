#!/bin/bash

PATH_INPUT='DoxygenConfigFortran'
PATH_TEMP='temp'

#modify the path for the source code depending
#on the relative path on the cluster
sed -e 's|^INPUT [ ]*=[ ]*[a-zA-Z0-9/\_"]*|'"INPUT                  = \"$augeanstables/src\""'|g' < $PATH_INPUT > $PATH_TEMP
mv $PATH_TEMP $PATH_INPUT

#modify the path for the documentation pictures
#depending on the relative path on the cluster
sed -e 's|^IMAGE_PATH [ ]*=[ ]*[a-zA-Z0-9/\_"]*|'"IMAGE_PATH             = \"$AUGEANSTABLES_CONFIG/doc_config/doc_pictures\""'|g' < $PATH_INPUT > $PATH_TEMP
mv $PATH_TEMP $PATH_INPUT

#exclude the surrogate_class.f as it makes the doxygen engine crash
sed -e 's|^EXCLUDE [ ]*=[ ]*[a-zA-Z0-9/\_"]*|'"EXCLUDE             = \"$augeanstables/src/field/surrogate_class.f\""'|g' < $PATH_INPUT > $PATH_TEMP
mv $PATH_TEMP $PATH_INPUT
