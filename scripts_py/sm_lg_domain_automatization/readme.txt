Author: J.L. Desmarais
email:  desmaraisjulien@gmail.com
----------------------------------------------------------------------
this folder contains python scripts useful to generate the comparison
between the small and large domain simulations
----------------------------------------------------------------------


automatization_csts.py:	  constants needed to generate the scripts
			  running the simualtions


library_sm_lg_inputs.py:  contain intermediate fcts useful to generate
			  the two input files needed to configure the
			  $augeanstables code to run the simulations
			  on the small and on the large domain


create_sm_lg_inputs.py:	  contain the function generating the two input
			  files needed to configure the $augeanstables
			  code to run simulations on the small and on the
			  large domain
			  can be also used an an executable using:

			  ex: python create_sm_lg_inputs.py -i model_input.txt
			      -s input_sm_domain.txt -l input_lg_domain.txt
			      -T 0.995 -v 0.1

			  ex: more info with python create_sm_lg_inputs.py -h


library_sm_lg_results.py: contain intermediate fcts useful to create the
			  main folder where the simulation on the small
			  and large domains will be run, generate the
			  executable for the small and large domains
			  simulations, generate the PBS scripts to run
			  the simulations and run the PBS scripts

			  - the main function can be used to create and
			    run simulations


library_sm_lg_error.py:	  contain intermediate fcts useful to generate the
			  error files by comparing the simulations on the
			  small and large domains
			  
			  - the main function can be used to create the error
			    graphs from the simulation output folders

library_sm_lg_graph.py:	  contain functions to contain a graph representing
			  the maximum of the error as function of time


