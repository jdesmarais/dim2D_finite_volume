#!/bin/bash

#test for bf_layer
./run_test_bf_layer.sh

#test for bf_sublayer
/home/jdesmarais/local/runtest/runtest.sh test_bf_sublayer_prog
sleep 2

#test for bf_mainlayer
/home/jdesmarais/local/runtest/runtest.sh test_bf_mainlayer_prog

#test for interface_abstract
./run_test_interface_abstract.sh