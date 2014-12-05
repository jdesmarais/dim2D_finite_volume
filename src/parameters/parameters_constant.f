      !> @file
      !> module containing the constants used in the
      !> program: they will be propagated at compilation
      !> time
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> definition of the constant used in the whole program
      !
      !> @date
      !> 08_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_constant

        !>program information and conventions
        character*(*) :: institut
        character*(*) :: prog_version
        character*(*) :: commit
        character*(*) :: ref
        character*(*) :: convention

        parameter (institut     = 'Eindhoven university of technology')
        parameter (prog_version = 'augeanstables V0.6')
        parameter (commit = '4696c130ba64ead184d8602706d740caf696cb66')
        parameter (ref          = 'desmaraisjulien@gmail.com')
        parameter (convention   = 'cf-1.6')
        

        !>main variable types
        integer, parameter :: scalar=0
        integer, parameter :: vector_x=1
        integer, parameter :: vector_y=2

        !>phase identification
        integer, parameter :: liquid=0
        integer, parameter :: vapor=1


        !>initial conditions choice for NS
        !integer, parameter :: steady_state=0
        integer, parameter :: peak=1
        integer, parameter :: vortex=2
        integer, parameter :: sym_x=3
        integer, parameter :: sym_y=4
        integer, parameter :: negative_spot=5

        !>initial conditions choice for DIM
        integer, parameter :: steady_state=0
        integer, parameter :: drop_retraction=1
        integer, parameter :: bubble_ascending=2
        integer, parameter :: homogeneous_liquid=3
        integer, parameter :: drop_collision=4
        integer, parameter :: phase_separation=5

        !>boundary conditions choice
        integer, parameter :: periodic_xy_choice=0
        integer, parameter :: reflection_xy_choice=1
        integer, parameter :: wall_xy_choice=2
        integer, parameter :: wall_x_reflection_y_choice=3
        integer, parameter :: hedstrom_xy_choice=4
        integer, parameter :: hedstrom_xy_corners_choice=5
        integer, parameter :: hedstrom_x_reflection_y_choice=6
        integer, parameter :: poinsot_xy_choice=7
        integer, parameter :: yoolodato_xy_choice=8

        !>boundary conditions type choice
        integer, parameter :: bc_nodes_choice=0
        integer, parameter :: bc_fluxes_choice=1
        integer, parameter :: bc_timedev_choice=2

        !>equations tuning choice
        integer, parameter :: no_gravity_choice=0
        integer, parameter :: earth_gravity_choice=1
        integer, parameter :: no_wave_forcing=0
        integer, parameter :: oscillatory_forcing=1
        integer, parameter :: intermittent_oscillatory_forcing=2

        !>i/o management choice
        integer, parameter :: netcdf_choice=0

        !>mpi constant
        integer, parameter :: N=1
        integer, parameter :: S=2
        integer, parameter :: E=3
        integer, parameter :: W=4
        integer, parameter :: interior=5

        integer, parameter :: x_direction=1
        integer, parameter :: y_direction=2
        integer, parameter :: xy_direction=3
        integer, parameter :: n1_direction=4
        integer, parameter :: n2_direction=5
        integer, parameter :: min_border=1
        integer, parameter :: max_border=2

        integer, parameter :: only_compute_proc=0
        integer, parameter :: compute_and_exchange_proc=1
        integer, parameter :: only_exchange_proc=2

        !>open b.c. constant
        logical, parameter :: left=.true.
        logical, parameter :: right=.false.

        logical, parameter :: inflow_type=.false.
        logical, parameter :: outflow_type=.true.

        integer, parameter :: always_inflow=0
        integer, parameter :: always_outflow=1
        integer, parameter :: ask_flow=2

      end module parameters_constant

