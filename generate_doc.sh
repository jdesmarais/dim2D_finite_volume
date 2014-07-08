#!/bin/bash

original_dir=$(pwd)

echo $original_dir

if [ -d "./doc" ]
then
rm -rf ./doc
fi
cd $AUGEANSTABLES_CONFIG/doc_config
./generate_doc.sh
mv doc $original_dir
