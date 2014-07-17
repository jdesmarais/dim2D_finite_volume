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


#enable the scripts for testing the erymanthianboar project
echo 'enable the scripts for testing the erymanthianboar project'
cd $augeanstables/src/test_files/test_erymanthianboar
chmod 777 *.sh


#modify the files for testing the erymanthianboar project
echo 'modify erymanthianboar test files for ifport module'
if [ "$AUGEANSTABLES_COMPILER" != "*ifort*" ]
    then
    sed -i 's/ use ifport/ !use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_layer_module.f
    sed -i 's/ use ifport/ !use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_layer_prog.f
    sed -i 's/ use ifport/ !use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_mainlayer_prog.f
else
    sed -i 's/!use ifport/ use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_layer_module.f
    sed -i 's/!use ifport/ use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_layer_prog.f
    sed -i 's/!use ifport/ use ifport/g' $augeanstables/src/test_files/test_erymanthianboar/test_bf_mainlayer_prog.f
fi

