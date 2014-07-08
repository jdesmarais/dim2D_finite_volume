###########################
#AUGEANSTABLES VERSION 1.0
###########################

author: Desmarais Julien
email : desmaraisjulien@gmail.com


0. Get started
===========================
The commands for generating proper executable are controlled by:
    - the module file augeanstables in ~/modules/augeanstables
    - the path for the librairies (netcdf, mpi, ...) are defined
      in ./src/test_files/config/

The commands and librairies links used when generating the executables can be checked using:
> cd ./src/test_files/test_erymanthianboar
> make test

1. File description
===========================
src             : fortran source code
scripts_py      : python scripts needed to postprocess the results
config.sh       : configure the scripts generating executables and documentation
generate_doc.sh : generate the html documentation for the source code


2. Generate and open documentation
==================================
To generate the documentation, use:
> ./generate_doc.sh

To open the documentation, use:
> firefox ./doc/html/files.html
	  
