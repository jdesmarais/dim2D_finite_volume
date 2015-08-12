============================================================
Author: J.L. Desmarais (desmaraisjulien@gmail.com)
Date  : 12/08/2015
============================================================
This directory contains the files used to run the simulations
with vapor bubbles next to a wall using DIM in 2D

# briefly
#------------------------------------------------------------
ex: run_steady_state_simulations.py
    -> run all the simulations to generate the results
       presented in the paper for test case [A]

ex: run_non_steady_state_simulations.py
    -> run all the simulations to generate the results
       presented in the paper for test cases [B,C,D]


# test cases simulated
#------------------------------------------------------------
[A] Study of the steady state of a vapor bubble deposited on
    a substrate. The initial shape of the vapor bubble is a
    half-disc in 2D. Because the contact angle imposed on the
    substrate is not 90 degrees, the bubble converges to the
    spherical cap shape imposed by the contact angle.

[B] Study of the nucleation of a vapor bubble on a substrate
    for varying contact angles and heat fluxes from the wall.
    There is no inflow in the computational domain.

[C] Study of the detachment of a vapor bubble from a substrate.
    The vapor bubble deposited on the substrate is initially at
    equilibrium. It has the shape of a spherical cap imposed by
    the contact angle of the substrate. The detachment is 
    studied for varying contact angles and incoming flow velocity.
    Especially, the magnitude of the maximum velocity in the
    inflow and its profile are varied (parabolic and linear)

[D] Study of the nucleation of a vapor bubble on a substrate in the
    presence of a non-zero incoming flow. The nucleation is studied
    for a contact angle of 22.5 degrees favorising its detachment. The
    magnitude and the profile of the incoming flow is varied as well
    as the heat flux provided by the wall


# files automatizing the simulation runs
#------------------------------------------------------------
automatization_wall_cst.py : constants used to create the inputs when
			     compiling the fortran executables running
			     the simulations

library_wall_st_inputs.py : useful functions when creating the input file
			    submitted to config.py (see footnote [1])
			    These inputs are for simulations studying the
			    steady state of a half-disc shaped 2D vapor
			    bubble deposited on a substrate
			    (see test case [A]).

create_wall_st_inputs.py : python script used to generate the input for
			   simulating test case [A]

create_wall_nonst_inputs.py : python script used to generate the inputs
			      to simulate test cases [B,C,D]

library_wall_st_results.py : useful functions to generate the simulations
			     for test case [A]

library_wall_nonst_results.py : useful functions to generate the simulations
			        for test case [B,C,D]

run_steady_state_simulations.py : run all the simulations to generate the
				  results presented in the paper for test
				  case [A]. Switches inside the python
				  script are helpful to restrain the
				  simulations

run_non_steady_state_simulations.py : run all the simulations to
				  generate the results presented in the
				  paper for test cases [B,C,D]. Switches
				  inside the python script are helpful
				  to restrain the simulations
			    
			    
# footnotes
#------------------------------------------------------------
[1] config.py is a python script used in $augeanstables/src/config
    to compile the fortran executable to run the simulation. It requires
    a text input file to configure the source code and generate the
    corresponding fortran executable

