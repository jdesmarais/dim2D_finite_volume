#!/bin/bash

inputs="../../parameters/parameters_input.f"
makefile="../makefile"

#makefile and 
./config.py -i inputs_tests.txt > /dev/null
./change_parameter.sh -i $inputs -o $inputs -p ne -v 4

#test field
echo ''
echo '**********************************'
echo 'test_field'
echo '**********************************'
cd ..
make clean > /dev/null
make test_field > /dev/null
./test_field
cd config

#test cg_operators
echo ''
echo '**********************************'
echo 'test_cg_operators'
echo '**********************************'
cd ..
make clean > /dev/null
make test_cg_operators > /dev/null
./test_cg_operators
cd config

#test dim2d_prim
echo ''
echo '**********************************'
echo 'test_dim2d_prim'
echo '**********************************'
cd ..
make clean > /dev/null
make test_dim2d_prim > /dev/null
./test_dim2d_prim
cd config

#test dim2d_fluxes
echo ''
echo '**********************************'
echo 'test_dim2d_fluxes'
echo '**********************************'
cd ..
make clean > /dev/null
make test_dim2d_fluxes > /dev/null
./test_dim2d_fluxes
cd config

#test dim2d_eq
echo ''
echo '**********************************'
echo 'test_dim2d_eq'
echo '**********************************'
cd ..
make clean > /dev/null
make test_dim2d_eq > /dev/null
./test_dim2d_eq
cd config

#test fv_operators
echo ''
echo '**********************************'
echo 'test_fv_operators'
echo '**********************************'
./change_parameter.sh -i $makefile -o $makefile -p pm_choice -v simpletest_choice
./change_parameter.sh -i $inputs -o $inputs -p ntx -v 10
./change_parameter.sh -i $inputs -o $inputs -p nty -v 6
./change_parameter.sh -i $inputs -o $inputs -p ne -v 1
./change_parameter.sh -i $inputs -o $inputs -p bc_choice -v periodic_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p bcx_type_choice -v bc_nodes_choice
./change_parameter.sh -i $inputs -o $inputs -p bcy_type_choice -v bc_nodes_choice
cd ..
make clean > /dev/null
make test_fv_operators > /dev/null 2>&1
./test_fv_operators
cd config

#test bc_periodic
echo ''
echo '**********************************'
echo 'test_bc_periodic'
echo '**********************************'
./config.py -i inputs_tests.txt > /dev/null
./change_parameter.sh -i $makefile -o $makefile -p bc_choice -v periodic_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p ntx -v 10
./change_parameter.sh -i $inputs -o $inputs -p nty -v 12
./change_parameter.sh -i $inputs -o $inputs -p ne -v 4
./change_parameter.sh -i $inputs -o $inputs -p bc_choice -v periodic_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p bcx_type_choice -v bc_nodes_choice
./change_parameter.sh -i $inputs -o $inputs -p bcy_type_choice -v bc_nodes_choice
cd ..
make clean > /dev/null
make test_bc_periodic > /dev/null 2>&1
./test_bc_periodic
cd config

#test bc_reflection
echo ''
echo '**********************************'
echo 'test_bc_reflection'
echo '**********************************'
./change_parameter.sh -i $makefile -o $makefile -p bc_choice -v reflection_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p ntx -v 10
./change_parameter.sh -i $inputs -o $inputs -p nty -v 12
./change_parameter.sh -i $inputs -o $inputs -p ne -v 4
./change_parameter.sh -i $inputs -o $inputs -p bc_choice -v reflection_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p bcx_type_choice -v bc_nodes_choice
./change_parameter.sh -i $inputs -o $inputs -p bcy_type_choice -v bc_nodes_choice
cd ..
make clean > /dev/null
make test_bc_reflection > /dev/null 2>&1
./test_bc_reflection
cd config

#test wall_xy_module
echo ''
echo '**********************************'
echo 'test_wall_xy_module'
echo '**********************************'
./change_parameter.sh -i $makefile -o $makefile -p bc_choice -v wall_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p ntx -v 10
./change_parameter.sh -i $inputs -o $inputs -p nty -v 12
./change_parameter.sh -i $inputs -o $inputs -p ne -v 4
./change_parameter.sh -i $inputs -o $inputs -p bc_choice -v wall_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p bcx_type_choice -v bc_fluxes_choice
./change_parameter.sh -i $inputs -o $inputs -p bcy_type_choice -v bc_fluxes_choice
cd ..
make clean > /dev/null
make test_wall_xy_module > /dev/null 2>&1
./test_wall_xy_module
cd config

#test bc_wall
echo ''
echo '**********************************'
echo 'test_bc_wall'
echo '**********************************'
cd ..
make clean > /dev/null
make test_bc_wall > /dev/null 2>&1
./test_bc_wall
cd config

#test rk3tvd
echo ''
echo '**********************************'
echo 'test_rk3tvd'
echo '**********************************'
./change_parameter.sh -i $makefile -o $makefile -p pm_choice -v simpletest_choice
./change_parameter.sh -i $makefile -o $makefile -p bc_choice -v periodic_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p ntx -v 10
./change_parameter.sh -i $inputs -o $inputs -p nty -v 6
./change_parameter.sh -i $inputs -o $inputs -p ne -v 1
./change_parameter.sh -i $inputs -o $inputs -p bc_choice -v periodic_xy_choice
./change_parameter.sh -i $inputs -o $inputs -p bcx_type_choice -v bc_nodes_choice
./change_parameter.sh -i $inputs -o $inputs -p bcy_type_choice -v bc_nodes_choice
cd ..
make clean > /dev/null
make test_rk3tvd > /dev/null 2>&1
./test_rk3tvd
cd config

#test nf90_operators
echo ''
echo '**********************************'
echo 'test_nf90_operators'
echo '**********************************'
./config.py -i inputs_tests.txt > /dev/null
./change_parameter.sh -i $inputs -o $inputs -p ne -v 4
cd ..
make clean > /dev/null
make test_nf90_operators > /dev/null 2>&1
./test_nf90_operators
cd config

#test dim2d_ic

#test rk3tvd_dim2d
