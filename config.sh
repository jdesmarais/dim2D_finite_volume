#!/bin/bash

#enable the script for generating executable
echo 'enable script to generate executable'
chmod 777 generate_exe.sh


#enable the script for generating documentation
echo 'enable script to generate documentation'
chmod 777 generate_doc.sh


#enable the scripts for configuring the source code
echo 'enable the scripts for configuring the source code'
cd $AUGEANSTABLES_CONFIG
chmod 777 *.sh
chmod 777 config.py


#enable the scripts for configuring the source code documentation
echo 'enable the scripts for configuring the source code documentation'
echo 'configure the generation of documentation'
cd $AUGEANSTABLES_CONFIG/doc_config
chmod 777 *.sh
./make_DoxygenConfigFortran.sh


#enable the scripts for generating default executables
echo 'enable the scripts for generating default executables'
cd $AUGEANSTABLES_CONFIG/default_gen_config
chmod 777 *.sh


#enable the scripts for testing the source code
echo 'enable the scripts for testing the source code'
cd $augeanstables/src/test_files
chmod 777 *.sh

