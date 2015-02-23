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
        parameter (commit = '1c4608fa5db3b574ab2f802e68f24a8bdc4a5e84')
        parameter (ref          = 'desmaraisjulien@gmail.com')
        parameter (convention   = 'cf-1.6')
        

        !>type of sd_operators
        integer, parameter :: sd_interior_type=0
        integer, parameter :: sd_L0_type=1
        integer, parameter :: sd_L1_type=2
        integer, parameter :: sd_R1_type=3
        integer, parameter :: sd_R0_type=4

        integer, parameter :: sd_interior_n_type=5
        integer, parameter :: sd_L0_n_type=6
        integer, parameter :: sd_L1_n_type=7
        integer, parameter :: sd_R1_n_type=8
        integer, parameter :: sd_R0_n_type=9


        !>main variable types
        integer, parameter :: scalar=0
        integer, parameter :: vector_x=1
        integer, parameter :: vector_y=2

        !>phase identification
        integer, parameter :: liquid=0
        integer, parameter :: vapor=1


        !>initial conditions choice for NS
        character(15), dimension(7), parameter :: ns2d_ic_code =[
     $       'steady_state   ',
     $       'peak           ',
     $       'vortex         ',
     $       'sym_x          ',
     $       'sym_y          ',
     $       'negative_spot  ',
     $       'not implemented']

        !integer, parameter :: steady_state=0
        integer, parameter :: peak=1
        integer, parameter :: vortex=2
        integer, parameter :: sym_x=3
        integer, parameter :: sym_y=4
        integer, parameter :: negative_spot=5

        !>initial conditions choice for DIM
        character(18), dimension(7), parameter :: dim2d_ic_code =[
     $       'steady_state      ',
     $       'drop_retraction   ',
     $       'bubble_ascending  ',
     $       'homogeneous_liquid',
     $       'drop_collision    ',
     $       'phase_separation  ',
     $       'bubble_transported']

        integer, parameter :: steady_state=0
        integer, parameter :: drop_retraction=1
        integer, parameter :: bubble_ascending=2
        integer, parameter :: homogeneous_liquid=3
        integer, parameter :: drop_collision=4
        integer, parameter :: phase_separation=5
        integer, parameter :: bubble_transported=6

        !>boundary conditions choice
        character(23), dimension(9), parameter :: bc_code =[
     $       'periodic_xy            ',
     $       'reflection_xy          ',
     $       'wall_xy                ',
     $       'wall_x_reflection_y    ',
     $       'hedstrom_xy            ',
     $       'hedstrom_xy_corners    ',
     $       'hedstrom_x_reflection_y',
     $       'poinsot_xy             ',
     $       'yoolodato_xy           ']

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
        integer, parameter :: moving_oscillatory_forcing=3

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


        !> open b.c. 
        !-------------------------------------------------------
        !control how the contribution of the outgoing waves 
        !is computed
        !-------------------------------------------------------
        ! obc_outgoing_cons: the contribution is computed using
        !                    the gradient of the conservative
        !                    variables
        !
        ! obc_outgoing_prim: the contribution is computed using
        !                    the gradient of the primitive
        !                    variables
        !-------------------------------------------------------
        integer, parameter :: obc_outgoing_cons = 0
        integer, parameter :: obc_outgoing_prim = 1


        !> open b.c. oneside fluxes
        !-------------------------------------------------------
        !control how the fluxes are computed at the edges
        !-------------------------------------------------------
        ! obc_edge_flux_capillarity:     the capillarity terms
        !                                are included in the
        !                                computation of the edge
        !                                fluxes
        !
        ! obc_edge_flux_no_capillarity:  the capillarity terms
        !                                are not_included in the
        !                                computation of the edge
        !                                fluxes
        !-------------------------------------------------------
        integer, parameter :: obc_edge_flux_capillarity = 0
        integer, parameter :: obc_edge_flux_no_capillarity = 1


        !> open b.c. eigenquantities
        !-------------------------------------------------------
        !control how the eigenquantities are computed for open b.c.
        !-------------------------------------------------------
        ! obc_eigenqties_bc:     the eigenquantities are computed
        !                        using the conservative variables
        !                        at the edge
        !
        ! obc_eigenqties_lin:    the eigenquantities are computed
        !                        using the values imposed in the
        !                        far field determined by the
        !                        initial conditions
        !
        ! obc_eigenqties_roe:    the eigenquantities are computed
        !                        using the Roe's average of the
        !                        values imposed in the far field
        !                        and the values at the edge
        !-------------------------------------------------------
        integer, parameter :: obc_eigenqties_bc  = 0
        integer, parameter :: obc_eigenqties_lin = 1
        integer, parameter :: obc_eigenqties_roe = 2


        !> open b.c. anti_corners
        !-------------------------------------------------------
        !control how the anti-corners are computed for open b.c.
        !-------------------------------------------------------
        ! obc_edge_xy_corner:    like there were corners
        !
        ! obc_edge_xy_flux:      the off-diagonal points are computed
        !                        as if there were N/S/E/W edge points
        !                        and the diagonal points are computed
        !                        like corners
        !
        ! obc_edge_xy_diag_flux: the quantities are computed in the
        !                        rotated frame such that the outward
        !                        direction to the boundary matches one
        !                        of the direction of the rotated frame
        !                        the terms in the outward direction
        !                        are approximated by the gradients and
        !                        in the transverse direction, there
        !                        are computed using the rotated
        !                        operators s_n1_L0...
        !-------------------------------------------------------
        integer, parameter :: obc_edge_xy_corner    = 0
        integer, parameter :: obc_edge_xy_flux      = 1
        integer, parameter :: obc_edge_xy_diag_flux = 2


        !> open b.c. edge_type
        !-----------------------------------------------------------
        !control how the type of open b.c. can be fixed at the edges
        !-----------------------------------------------------------
        ! always_inflow:    force the open b.c. to be always inflow
        !
        ! always_outflow:   force the open b.c. to be always outflow
        !
        ! ask_flow:         ask the open b.c to check the velocity
        !                   vector at the edge to know whether the
        !                   open b.c is inflow or outflow
        !-----------------------------------------------------------
        character(14), dimension(3), parameter :: obc_type_code =[
     $       'always_inflow ',
     $       'always_outflow',
     $       'ask_flow      ']
        integer, parameter :: always_inflow  = 0
        integer, parameter :: always_outflow = 1
        integer, parameter :: ask_flow       = 2


      end module parameters_constant

