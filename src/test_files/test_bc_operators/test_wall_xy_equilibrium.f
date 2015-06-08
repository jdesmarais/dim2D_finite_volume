      program test_wall_xy_equilibrium

        use check_data_module, only :
     $     is_real_validated,
     $     is_real_vector_validated

        use dim2d_parameters, only :
     $       cv_r,
     $       Re,Pr,We

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid

        use parameters_constant, only :
     $       right,
     $       sd_interior_type

        use parameters_kind, only :
     $       rkind

        use ridders_method_fcts_module, only :
     $       root_fct2,
     $       root_fct4

        use sd_operators_class, only : 
     $       sd_operators

        use wall_xy_equilibrium_module, only :
     $       dmddx,
     $       temperature,
     $       md_average,
     $       temperature_average,
     $       dwallInternalEnergy_dmd,
     $       
     $       wall_x_equilibrium_root_fct,
     $       wall_x_root_fct,
     $       get_wall_x_root_brackets,
     $       get_wall_x_md_ghost_cell,
     $       get_wall_heat_flux,
     $       get_wall_micro_contact_angle,
     $       
     $       compute_wall_E_ghost_cell,
     $       compute_wall_W_ghost_cell,
     $       compute_wall_N_ghost_cell,
     $       compute_wall_S_ghost_cell,
     $       
     $       flux_x_inviscid_momentum_x,
     $       flux_x_capillarity_momentum_x,
     $       flux_x_viscid_momentum_y,
     $       flux_x_capillarity_momentum_y,
     $       
     $       compute_wall_flux_x,
     $       compute_wall_flux_y


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        call test_inputs()

        test_loc = test_dmddx(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_dmddx: '',L1)', test_loc
        print '()'

        test_loc = test_temperature(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_temperature: '',L1)', test_loc
        print '()'

        test_loc = test_md_average(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_md_average: '',L1)', test_loc
        print '()'

        test_loc = test_temperature_average(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_temperature_average: '',L1)', test_loc
        print '()'

        test_loc = test_dwallInternalEnergy_dmd(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_dwallInternalEnergy_dmd: '',L1)', test_loc
        print '()'

        test_loc = test_wall_x_equilibrium_root_fct(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wall_x_equilibrium_root_fct: '',L1)', test_loc
        print '()'

        test_loc = test_wall_x_root_fct(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wall_x_root_fct: '',L1)', test_loc
        print '()'

        test_loc = test_get_wall_x_root_brackets(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_wall_x_root_brackets: '',L1)', test_loc
        print '()'

        test_loc = test_wall_x_equilibrium_root_fct(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wall_x_equilibrium_root_fct: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_E_ghost_cell(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_E_ghost_cell: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_W_ghost_cell(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_W_ghost_cell: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_N_ghost_cell(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_N_ghost_cell: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_S_ghost_cell(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_S_ghost_cell: '',L1)', test_loc
        print '()'

        test_loc = test_flux_x_inviscid_momentum_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_flux_x_inviscid_momentum_x: '',L1)', test_loc
        print '()'

        test_loc = test_flux_x_capillarity_momentum_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_flux_x_capillarity_momentum_x: '',L1)', test_loc
        print '()'

        test_loc = test_flux_x_viscid_momentum_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_flux_x_viscid_momentum_y: '',L1)', test_loc
        print '()'

        test_loc = test_flux_x_capillarity_momentum_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_flux_x_capillarity_momentum_y: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_flux_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_flux_x: '',L1)', test_loc
        print '()'

        test_loc = test_compute_wall_flux_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_wall_flux_y: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated

        contains


        function test_compute_wall_flux_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          type(sd_operators)            :: s

          real(rkind), dimension(4)     :: flux_x
          real(rkind), dimension(4)     :: flux_y

          
          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          flux_x = compute_wall_flux_x(nodes,s,4,3,0.0d0,
     $         [0.0d0,dx,2*dx,3*dx,4*dx],
     $         [0.0d0,dy,2*dy,3*dy,4*dy],
     $         right)


          call transpose_nodes(dx,dy,nodes)

          flux_y = compute_wall_flux_y(nodes,s,3,4,0.0d0,
     $         [0.0d0,dx,2*dx,3*dx,4*dx],
     $         [0.0d0,dy,2*dy,3*dy,4*dy],
     $         right)


          test_validated = is_real_vector_validated(
     $         flux_y,
     $         [flux_x(1), flux_x(3), flux_x(2), flux_x(4)],
     $         detailled)

        end function test_compute_wall_flux_y


        function test_compute_wall_flux_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          type(sd_operators)            :: s
          
          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          test_validated = is_real_vector_validated(
     $         compute_wall_flux_x(nodes,s,4,3,0.0d0,
     $         [0.0d0,dx,2*dx,3*dx,4*dx],
     $         [0.0d0,dy,2*dy,3*dy,4*dy],
     $         right),
     $         [0.0d0, 222.52820459905700d0, -82.09240176344290d0, 0.001d0],
     $         detailled)

        end function test_compute_wall_flux_x


        function test_flux_x_capillarity_momentum_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          type(sd_operators)            :: s
          
          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          test_validated = is_real_validated(
     $         flux_x_capillarity_momentum_y(nodes,s,4,3,dx,dy),
     $         786.25d0,
     $         detailled)

        end function test_flux_x_capillarity_momentum_y


        function test_flux_x_viscid_momentum_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          type(sd_operators)            :: s
          

          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          test_validated = is_real_validated(
     $         flux_x_viscid_momentum_y(nodes,s,4,3,dx,dy),
     $         17.3370088172145d0,
     $         detailled)

        end function test_flux_x_viscid_momentum_y


        function test_flux_x_capillarity_momentum_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          type(sd_operators)            :: s
          

          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          test_validated = is_real_validated(
     $         flux_x_capillarity_momentum_x(nodes,s,4,3,dx,dy),
     $         -2798.8046875d0,
     $         detailled)

        end function test_flux_x_capillarity_momentum_x


        function test_flux_x_inviscid_momentum_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          
          real(rkind), dimension(5) :: x_map
          real(rkind), dimension(5) :: y_map          

          call get_nodes_for_flux_test(
     $         nodes,dx,dy)

          x_map = [0.0d0, dx, 2.0d0*dx, 3.0d0*dx, 4.0d0*dx]
          y_map = [0.0d0, dy, 2.0d0*dy, 3.0d0*dy, 4.0d0*dy]

          test_validated = is_real_validated(
     $         flux_x_inviscid_momentum_x(
     $            0.0d0,x_map,y_map,nodes,dx,dy,4,3,right),
     $         -57.35226415094340d0,
     $         detailled)

        end function test_flux_x_inviscid_momentum_x


        subroutine get_nodes_for_flux_test(nodes,dx,dy)

          implicit none

          real(rkind), dimension(5,5,4), intent(out) :: nodes
          real(rkind)                  , intent(out) :: dx
          real(rkind)                  , intent(out) :: dy

          integer :: i,j,k

          do k=1,4
             do j=1,size(nodes,2)
                do i=1,size(nodes,1)
                   nodes(i,j,k) = -99.0e30
                end do
             end do
          end do

          nodes(2:5,2,1) = [-7.5d0, 8.1d0,4.6d0, 1.2d0]
          nodes(2:5,3,1) = [ 1.9d0,-2.3d0,5.1d0,-3.6d0]
          nodes(2:5,4,1) = [ 8.1d0, 1.0d0,3.2d0, 5.1d0]
          
          nodes(3:4,2,2) = [-2.3d0,3.2d0]
          nodes(3:4,3,2) = [ 0.5d0,6.8d0]
          nodes(3:4,4,2) = [-0.1d0,0.56d0]

          nodes(3:4,3,3) = [-0.2d0,9.5d0]
          nodes(3,3,4)   = -1.17991847826087d0

          dx = 0.1d0
          dy = 0.2d0

        end subroutine get_nodes_for_flux_test


        function test_compute_wall_S_ghost_cell(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(5)     :: x_map
          real(rkind), dimension(5)     :: y_map
          real(rkind), dimension(5,5,4) :: nodes

          real(rkind) :: velocity_x
          real(rkind) :: velocity_y

          real(rkind) :: T_guess
          real(rkind) :: md_guess

          real(rkind) :: md_ghost_cell
          real(rkind) :: qx_ghost_cell
          real(rkind) :: qy_ghost_cell
          real(rkind) :: md_ghost_cell2

          logical :: test_loc

          !compute ghost cell E
          call get_test_parameters_wall_ghost_cell(
     $         T_guess,
     $         md_guess,
     $         velocity_x,
     $         velocity_y,
     $         x_map,
     $         y_map,
     $         nodes)          
          
          call compute_wall_E_ghost_cell(
     $         4,3,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          md_ghost_cell  = nodes(4,3,1)
          qx_ghost_cell  = nodes(4,3,2)
          qy_ghost_cell  = nodes(4,3,3)
          md_ghost_cell2 = nodes(5,3,1)

          !test S ghost cell
          call transpose_nodes(
     $         velocity_x,
     $         velocity_y,
     $         nodes)
          
          call reflect_y(
     $         velocity_y,
     $         nodes)

          call compute_wall_S_ghost_cell(
     $         3,2,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          !test
          test_validated = .true.

          test_loc = is_real_validated(
     $         nodes(3,2,1),
     $         md_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,2,2),
     $         qy_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qx_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,2,3),
     $         -qx_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qy_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,1,1),
     $         md_ghost_cell2,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell2 failed'')'
          end if

        end function test_compute_wall_S_ghost_cell


        function test_compute_wall_N_ghost_cell(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(5)     :: x_map
          real(rkind), dimension(5)     :: y_map
          real(rkind), dimension(5,5,4) :: nodes

          real(rkind) :: velocity_x
          real(rkind) :: velocity_y

          real(rkind) :: T_guess
          real(rkind) :: md_guess

          real(rkind) :: md_ghost_cell
          real(rkind) :: qx_ghost_cell
          real(rkind) :: qy_ghost_cell
          real(rkind) :: md_ghost_cell2

          logical :: test_loc

          !compute ghost cell E
          call get_test_parameters_wall_ghost_cell(
     $         T_guess,
     $         md_guess,
     $         velocity_x,
     $         velocity_y,
     $         x_map,
     $         y_map,
     $         nodes)
          
          call compute_wall_E_ghost_cell(
     $         4,3,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          md_ghost_cell  = nodes(4,3,1)
          qx_ghost_cell  = nodes(4,3,2)
          qy_ghost_cell  = nodes(4,3,3)
          md_ghost_cell2 = nodes(5,3,1)

          !reflection along x to test N ghost cell
          call transpose_nodes(
     $         velocity_x,
     $         velocity_y,
     $         nodes)

          call compute_wall_N_ghost_cell(
     $         3,4,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          !test
          test_validated = .true.

          test_loc = is_real_validated(
     $         nodes(3,4,1),
     $         md_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,4,2),
     $         qy_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qx_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,4,3),
     $         qx_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qy_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(3,5,1),
     $         md_ghost_cell2,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell2 failed'')'
          end if


        end function test_compute_wall_N_ghost_cell


        function test_compute_wall_W_ghost_cell(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(5)     :: x_map
          real(rkind), dimension(5)     :: y_map
          real(rkind), dimension(5,5,4) :: nodes

          real(rkind) :: velocity_x
          real(rkind) :: velocity_y

          real(rkind) :: T_guess
          real(rkind) :: md_guess

          real(rkind) :: md_ghost_cell
          real(rkind) :: md_ghost_cell2
          real(rkind) :: qx_ghost_cell
          real(rkind) :: qy_ghost_cell

          logical :: test_loc

          !compute ghost cell E
          call get_test_parameters_wall_ghost_cell(
     $         T_guess,
     $         md_guess,
     $         velocity_x,
     $         velocity_y,
     $         x_map,
     $         y_map,
     $         nodes)          
          
          call compute_wall_E_ghost_cell(
     $         4,3,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          md_ghost_cell  = nodes(4,3,1)
          qx_ghost_cell  = nodes(4,3,2)
          qy_ghost_cell  = nodes(4,3,3)
          md_ghost_cell2 = nodes(5,3,1)


          !reflection along x to test W ghost cell
          call reflect_x(
     $         velocity_x,
     $         nodes)

          call compute_wall_W_ghost_cell(
     $         2,3,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          !test
          test_validated = .true.

          test_loc = is_real_validated(
     $         nodes(2,3,1),
     $         md_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(2,3,2),
     $         -qx_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qx_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(2,3,3),
     $         qy_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qy_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(1,3,1),
     $         md_ghost_cell2,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell2 failed'')'
          end if

        end function test_compute_wall_W_ghost_cell


        function test_compute_wall_E_ghost_cell(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(5)     :: x_map
          real(rkind), dimension(5)     :: y_map
          real(rkind), dimension(5,5,4) :: nodes

          real(rkind) :: velocity_x
          real(rkind) :: velocity_y

          real(rkind) :: T_guess
          real(rkind) :: md_guess

          real(rkind) :: md_ghost_cell

          logical :: test_loc

          call get_test_parameters_wall_ghost_cell(
     $         T_guess,
     $         md_guess,
     $         velocity_x,
     $         velocity_y,
     $         x_map,
     $         y_map,
     $         nodes)
          
          call compute_wall_E_ghost_cell(
     $         4,3,
     $         0.0d0,
     $         x_map,
     $         y_map,
     $         nodes,
     $         T_guess,
     $         md_guess,
     $         sd_interior_type)

          md_ghost_cell = get_wall_x_md_ghost_cell(
     $         T_guess,
     $         md_guess,
     $         [nodes(2,3,1),nodes(3,3,1)],
     $         [nodes(3,2,1),nodes(3,4,1)],
     $         velocity_x,
     $         velocity_y,
     $         nodes(3,3,4),
     $         x_map(2)-x_map(1),
     $         y_map(2)-y_map(1),
     $         90.0d0,
     $         0.005d0,
     $         sd_interior_type,
     $         sd_interior_type)

          test_validated = .true.

          test_loc = is_real_validated(
     $         nodes(4,3,1),
     $         nodes(3,3,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell : no equal reflection'')'
          end if

          test_loc = is_real_validated(
     $         nodes(4,3,1),
     $         md_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(4,3,2),
     $         md_ghost_cell*velocity_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qx_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(4,3,3),
     $         -md_ghost_cell*velocity_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''qy_ghost_cell failed'')'
          end if

          test_loc = is_real_validated(
     $         nodes(5,3,1),
     $         nodes(2,3,1)-3*nodes(3,3,1)+3*md_ghost_cell,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''md_ghost_cell2 failed'')'
          end if

        end function test_compute_wall_E_ghost_cell


        subroutine get_test_parameters_wall_ghost_cell(
     $     T_guess,
     $     md_guess,
     $     velocity_x,
     $     velocity_y,
     $     x_map,
     $     y_map,
     $     nodes)

          implicit none

          real(rkind)                  , intent(out) :: T_guess
          real(rkind)                  , intent(out) :: md_guess
          real(rkind)                  , intent(out) :: velocity_x
          real(rkind)                  , intent(out) :: velocity_y
          real(rkind), dimension(5)    , intent(out) :: x_map
          real(rkind), dimension(5)    , intent(out) :: y_map
          real(rkind), dimension(5,5,4), intent(out) :: nodes

          real(rkind) :: dx
          real(rkind) :: dy

          T_guess  = 0.915
          md_guess = 1.45d0

          velocity_x = 0.02d0
          velocity_y = 0.03d0

          dx = 0.5
          dy = 0.8

          x_map = [0.0,dx,2*dx,3*dx,4*dx]
          y_map = [0.0,dy,2*dy,3*dy,4*dy]

          nodes(2,3,1) = 0.2d0
          nodes(3,3,1) = 0.5d0

          nodes(3,2,1) = 1.2d0
          nodes(3,4,1) = 0.9d0

          nodes(3,3,2) = nodes(3,3,1)*velocity_x
          nodes(3,3,3) = nodes(3,3,1)*velocity_y
          nodes(3,3,4) = 2.26d0

        end subroutine get_test_parameters_wall_ghost_cell


        function test_get_wall_x_root_brackets(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc

          type(root_fct2) :: root_fct2_used
          type(root_fct4) :: root_fct4_used

          real(rkind) :: md_vap
          real(rkind) :: md_liq

          test_validated = .true.

          md_vap = get_mass_density_vapor(0.95d0)
          md_liq = get_mass_density_liquid(0.95d0)

          test_loc = is_real_vector_validated(
     $         get_wall_x_root_brackets(root_fct2_used,0.95d0,0.0d0),
     $         [md_vap,md_liq],
     $         detailled)
          test_validated = test_validated.and.test_loc


          test_loc = is_real_vector_validated(
     $         get_wall_x_root_brackets(root_fct4_used,0.95d0,1.5d0),
     $         [1.5d0,2.01d0],
     $         detailled)
          test_validated = test_validated.and.test_loc          

        end function test_get_wall_x_root_brackets
      

        function test_wall_x_root_fct(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(wall_x_root_fct) :: wall_x_root_fct_used

          call wall_x_root_fct_used%ini(
     $         [0.2d0,1.5d0],
     $         [1.2d0,0.9d0],
     $         0.02d0,
     $         0.03d0,
     $         2.26d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.005d0,
     $         sd_interior_type,
     $         sd_interior_type)
          
          test_validated = is_real_validated(
     $         wall_x_root_fct_used%f(0.5d0),
     $         wall_x_equilibrium_root_fct(
     $         1.5d0,
     $         [0.2d0,0.5d0],
     $         [1.2d0,0.9d0],
     $         0.02d0,
     $         0.03d0,
     $         2.26d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.005d0,
     $         sd_interior_type,
     $         sd_interior_type),
     $         detailled)

        end function test_wall_x_root_fct


        function test_wall_x_equilibrium_root_fct(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         wall_x_equilibrium_root_fct(
     $         1.5d0,
     $         [0.2d0,0.5d0],
     $         [1.2d0,0.9d0],
     $         0.02d0,
     $         0.03d0,
     $         2.26d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.005d0,
     $         sd_interior_type,
     $         sd_interior_type),
     $         -0.257602898833856d0,
     $         detailled)

        end function test_wall_x_equilibrium_root_fct


        function test_dwallInternalEnergy_dmd(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         dwallInternalEnergy_dmd(1.5d0,0.95d0,45.0d0),
     $         0.007151129d0,
     $         detailled)

        end function test_dwallInternalEnergy_dmd
        
        function test_temperature_average(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         temperature_average(1.5d0,0.2d0,5.0d0),
     $         11.5d0,
     $         detailled)

        end function test_temperature_average

        function test_md_average(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         md_average(1.5d0,0.5d0),
     $         1.0d0,
     $         detailled)

        end function test_md_average


        function test_temperature(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         temperature(1.5d0,50.0d0,-0.1d0,0.2d0,6.3d0),
     $         1.0512500000d0,
     $         detailled)

        end function test_temperature


        function test_dmddx(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         dmddx(0.5d0,[-2.0d0,-99.0d0,3.0d0],sd_interior_type),
     $         5.0d0,
     $         detailled)

        end function test_dmddx


        subroutine reflect_x(
     $     velocity_x,
     $     nodes)

          implicit none

          real(rkind)                  , intent(inout) :: velocity_x
          real(rkind), dimension(5,5,4), intent(inout) :: nodes

          real(rkind)           :: temp
          integer, dimension(4) :: prefactor
          integer               :: j,k
          
          prefactor = [1,-1,1,1]

          velocity_x = -velocity_x

          do k=1,4
             do j=1,5
                temp         = nodes(1,j,k)
                nodes(1,j,k) = prefactor(k)*nodes(5,j,k)
                nodes(5,j,k) = prefactor(k)*temp

                temp         = nodes(2,j,k)
                nodes(2,j,k) = prefactor(k)*nodes(4,j,k)
                nodes(4,j,k) = prefactor(k)*temp                
                
                nodes(3,j,k) = prefactor(k)*nodes(3,j,k)

             end do
          end do

        end subroutine reflect_x


        subroutine reflect_y(
     $     velocity_y,
     $     nodes)

          implicit none

          real(rkind)                  , intent(inout) :: velocity_y
          real(rkind), dimension(5,5,4), intent(inout) :: nodes

          real(rkind)           :: temp
          integer, dimension(4) :: prefactor
          integer               :: i,k
          
          prefactor = [1,1,-1,1]

          velocity_y = -velocity_y

          do k=1,4
             do i=1,5
                temp         = nodes(i,1,k)
                nodes(i,1,k) = prefactor(k)*nodes(i,5,k)
                nodes(i,5,k) = prefactor(k)*temp

                temp         = nodes(i,2,k)
                nodes(i,2,k) = prefactor(k)*nodes(i,4,k)
                nodes(i,4,k) = prefactor(k)*temp

                nodes(i,3,k) = prefactor(k)*nodes(i,3,k)
             end do
          end do

        end subroutine reflect_y


        subroutine transpose_nodes(
     $     velocity_x,
     $     velocity_y,
     $     nodes)

          implicit none

          real(rkind)                  , intent(inout) :: velocity_x
          real(rkind)                  , intent(inout) :: velocity_y
          real(rkind), dimension(5,5,4), intent(inout) :: nodes

          real(rkind)                 :: temp
          real(rkind), dimension(5,5) :: temp_nodes
          integer                     :: k
          
          temp       = velocity_x
          velocity_x = velocity_y
          velocity_y = temp

          do k=1,4
             temp_nodes = nodes(:,:,k)
             nodes(:,:,k) = transpose(temp_nodes)
          end do

          temp_nodes   = nodes(:,:,2)
          nodes(:,:,2) = nodes(:,:,3)
          nodes(:,:,3) = temp_nodes

        end subroutine transpose_nodes


        subroutine test_inputs()

          implicit none

          logical :: test_loc
          logical :: test_validated

          test_validated = .true.

          test_loc = is_real_validated(cv_r,2.5d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''cv_r should equal 2.5d0'')'
          end if

          test_loc = is_real_validated(We,10.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''we should equal 10.0d0'')'
          end if

          test_loc = is_real_validated(Pr,20.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''Pr should equal 20.0d0'')'
          end if

          test_loc = is_real_validated(Re,5.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''Re should equal 5.0d0'')'
          end if

          test_loc = is_real_validated(
     $         get_wall_micro_contact_angle(0.0d0,0.0d0),
     $         90.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''wall_contact_angle should equal 90.0d0'')'
          end if

          test_loc = is_real_validated(
     $         get_wall_heat_flux(0.0d0,0.0d0,0.0d0),
     $         0.005d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''wall_heat_flux should equal 1.0d0'')'
          end if

          if(.not.test_validated) then
             stop ''
          end if

        end subroutine test_inputs

      end program test_wall_xy_equilibrium
