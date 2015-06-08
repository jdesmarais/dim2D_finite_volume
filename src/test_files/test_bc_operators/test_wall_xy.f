      program test_wall_xy

        use bc_operators_wall_xy_class, only :
     $       bc_operators_wall_xy

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated

        use dim2d_parameters, only :
     $       Re,Pr,We,cv_r

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NW_corner_type,
     $       NE_corner_type

        use parameters_input, only :
     $       nx,ny,ne,
     $       wall_micro_contact_angle,
     $       wall_maximum_heat_flux

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.


        call test_inputs()


        test_loc = test_apply_bc_on_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_nodes: '',L1)', test_loc
        print '()'

        test_loc = test_apply_bc_on_fluxes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_fluxes: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated


        contains

        function test_apply_bc_on_fluxes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bc_operators_wall_xy)         :: bc_operators_used
          type(sd_operators)                 :: s
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx+1,ny,ne) :: flux_x_test
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(nx,ny+1,ne) :: flux_y_test

          integer               :: k
          integer, dimension(4) :: prefactor_x
          integer, dimension(4) :: prefactor_y

          logical        :: test_loc

          prefactor_x = [1,-1, 1,1]
          prefactor_y = [1, 1,-1,1]


          test_validated = .true.

          ! initialize the nodes
          call ini_nodes_for_apply_bc_on_fluxes(
     $         x_map, y_map, nodes)


          ! compute the fluxes for the E_edge
          call bc_operators_used%apply_bc_on_fluxes(
     $         [E_edge_type,9,3,8],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          flux_x_test = flux_x


          ! compute the fluxes for the W_edge
          call reflect_x_nodes(nodes)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [W_edge_type,1,3,8],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_x(3,3:8,k),
     $            -prefactor_x(k)*flux_x_test(9,3:8,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_W('',I2,'') failed'')',k
             end if
          end do


          ! compute the fluxes for the S_edge
          call transpose_nodes(nodes)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [S_edge_type,3,1,8],
     $         0.0d0,y_map,x_map,nodes,s,
     $         flux_x,flux_y)

          flux_y_test(3:8,3,1) = flux_x(3,3:8,1)
          flux_y_test(3:8,3,2) = flux_x(3,3:8,3)
          flux_y_test(3:8,3,3) = flux_x(3,3:8,2)
          flux_y_test(3:8,3,4) = flux_x(3,3:8,4)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_y(3:8,3,k),
     $            flux_y_test(3:8,3,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_S('',I2,'') failed'')',k
             end if
          end do


          ! compute the fluxes for the N_edge
          call reflect_y_nodes(nodes)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [N_edge_type,3,9,8],
     $         0.0d0,y_map,x_map,nodes,s,
     $         flux_x,flux_y)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_y(3:8,9,k),
     $           -prefactor_y(k)*flux_y_test(3:8,3,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_N('',I2,'') failed'')',k
             end if
          end do


          ! compute the fluxes for the NW_corner
          call bc_operators_used%apply_bc_on_fluxes(
     $         [NW_corner_type,1,9,0],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          flux_y_test = flux_y

          ! compute the fluxes for the NE corner
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(flux_y_test)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [NE_corner_type,9,9,0],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_y(9:10,9,k),
     $            flux_y_test(9:10,9,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_NE('',I2,'') failed'')',k
             end if
          end do

          ! compute the fluxes for the SE corner
          call reflect_y_nodes(nodes)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [SE_corner_type,9,1,0],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_y(9:10,3,k),
     $            -prefactor_y(k)*flux_y_test(9:10,9,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_SE('',I2,'') failed'')',k
             end if
          end do

          ! compute the fluxes for the SW corner
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(flux_y_test)

          call bc_operators_used%apply_bc_on_fluxes(
     $         [SW_corner_type,1,1,0],
     $         0.0d0,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

          do k=1,ne
             test_loc = is_real_vector_validated(
     $            flux_y(1:2,3,k),
     $           -prefactor_y(k)*flux_y_test(1:2,9,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(.not.test_loc) then
                print '(''flux_SW('',I2,'') failed'')',k
             end if
          end do


        end function test_apply_bc_on_fluxes


        function test_apply_bc_on_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bc_operators_wall_xy)       :: bc_operators_used
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind), dimension(nx,ny,ne) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: nodes

          logical                          :: test_loc
          real(rkind), dimension(nx,ny,ne) :: nodes_test

          integer(ikind) :: i,j
          integer        :: k


          test_validated = .true.


          ! initialize the nodes for the computation
          ! of the E boundary conditions
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)

          ! apply the boundary conditions on the E layer
          call bc_operators_used%apply_bc_on_nodes(
     $         [E_edge_type,9,3,8],
     $         0.0d0,x_map,y_map,nodes_tmp,p_model,
     $         nodes)

          nodes_test = nodes          

          !print *, nodes(9,3:8,1)


          ! reflect the nodes_test in the x-direction
          ! to test the W layer
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)

          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(nodes_test)


          ! apply the boundary conditions on the W layer
          call bc_operators_used%apply_bc_on_nodes(
     $         [W_edge_type,1,3,8],
     $         0.0d0,x_map,y_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(2,3:8,1),
     $         nodes_test(2,3:8,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_W(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(2,3:8,2),
     $         nodes_test(2,3:8,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_W(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(2,3:8,3),
     $         nodes_test(2,3:8,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_W(3): failed'')'
          end if


          ! transpose the nodes_test
          ! to test the S layer
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)

          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)

          call transpose_nodes(nodes_test)
          

          ! apply the boundary conditions on the S layer
          call bc_operators_used%apply_bc_on_nodes(
     $         [S_edge_type,3,1,8],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(3:8,2,1),
     $         nodes_test(3:8,2,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_S(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(3:8,2,2),
     $         nodes_test(3:8,2,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_S(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(3:8,2,3),
     $         nodes_test(3:8,2,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_S(3): failed'')'
          end if


          ! reflect the nodes_test
          ! to test the N layer
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)

          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)
          call reflect_y_nodes(nodes)

          call reflect_y_nodes(nodes_test)
          

          ! apply the boundary conditions on the N layer
          call bc_operators_used%apply_bc_on_nodes(
     $         [N_edge_type,3,9,8],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(3:8,9,1),
     $         nodes_test(3:8,9,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_N(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(3:8,9,2),
     $         nodes_test(3:8,9,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_N(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(3:8,9,3),
     $         nodes_test(3:8,9,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_N(3): failed'')'
          end if

          ! compute the SW corner boundary points
          nodes = reshape((/ (((-99.0d0,i=1,10),j=1,10),k=1,ne) /),
     $         (/10,10,ne/))

          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)

          call bc_operators_used%apply_bc_on_nodes(
     $         [SW_corner_type,1,1,0],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          nodes_test = nodes


          ! SE corner boundary points
          nodes = reshape((/ (((-99.0d0,i=1,10),j=1,10),k=1,ne) /),
     $         (/10,10,ne/))

          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)

          call reflect_x_nodes(nodes_test)

          call bc_operators_used%apply_bc_on_nodes(
     $         [SE_corner_type,9,1,0],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(9:10,2,1),
     $         nodes_test(9:10,2,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_SW(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(9:10,2,2),
     $         nodes_test(9:10,2,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_SW(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(9:10,2,3),
     $         nodes_test(9:10,2,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_SW(3): failed'')'
          end if


          ! NE corner boundary points
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)

          call reflect_y_nodes(nodes_tmp)
          call reflect_y_nodes(nodes)

          call reflect_y_nodes(nodes_test)

          call bc_operators_used%apply_bc_on_nodes(
     $         [NE_corner_type,9,9,0],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(9:10,9,1),
     $         nodes_test(9:10,9,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NE(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(9:10,9,2),
     $         nodes_test(9:10,9,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NE(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(9:10,9,3),
     $         nodes_test(9:10,9,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NE(3): failed'')'
          end if


          ! NW corner boundary points
          call ini_nodes_for_apply_bc_on_nodes(
     $         x_map,y_map,nodes_tmp,nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call transpose_nodes(nodes_tmp)
          call transpose_nodes(nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)
          call reflect_y_nodes(nodes_tmp)
          call reflect_y_nodes(nodes)
          call reflect_x_nodes(nodes_tmp)
          call reflect_x_nodes(nodes)

          call reflect_x_nodes(nodes_test)

          call bc_operators_used%apply_bc_on_nodes(
     $         [NW_corner_type,1,9,0],
     $         0.0d0,y_map,x_map,nodes_tmp,p_model,
     $         nodes)

          test_loc = is_real_vector_validated(
     $         nodes(1:2,9,1),
     $         nodes_test(1:2,9,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NW(1): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(1:2,9,2),
     $         nodes_test(1:2,9,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NW(2): failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         nodes(1:2,9,3),
     $         nodes_test(1:2,9,3),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''apply_bc_on_nodes_NW(3): failed'')'
          end if
                    

        end function test_apply_bc_on_nodes


        subroutine reflect_x_nodes(nodes)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          integer, dimension(4) :: prefactor
          integer(ikind)        :: i,j
          integer               :: k
          integer(ikind)        :: i_r
          real(rkind)           :: tmp

          prefactor = [1,-1,1,1]
          
          do k=1,ne
             do j=1,size(nodes,2)
                do i=1,size(nodes,1)/2

                   i_r = size(nodes,1)-(i-1)

                   tmp                = nodes(i,j,k)
                   nodes(i,j,k)       = prefactor(k)*nodes(i_r,j,k)
                   nodes(i_r,j,k)     = prefactor(k)*tmp

                end do
             end do
          end do

        end subroutine reflect_x_nodes


        subroutine reflect_y_nodes(nodes)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          integer, dimension(4) :: prefactor
          integer(ikind)        :: i,j
          integer               :: k
          integer(ikind)        :: j_r
          real(rkind)           :: tmp

          prefactor = [1,1,-1,1]
          
          do k=1,ne
             do j=1,size(nodes,2)/2

                j_r = size(nodes,2)-(j-1)

                do i=1,size(nodes,1)

                   tmp            = nodes(i,j,k)
                   nodes(i,j,k)   = prefactor(k)*nodes(i,j_r,k)
                   nodes(i,j_r,k) = prefactor(k)*tmp

                end do
             end do
          end do

        end subroutine reflect_y_nodes


        subroutine transpose_nodes(nodes)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer                       :: k
          real(rkind), dimension(nx,ny) :: tmp
          
          do k=1,ne
             nodes(:,:,k) = transpose(nodes(:,:,k))
          end do

          tmp = nodes(:,:,2)
          nodes(:,:,2) = nodes(:,:,3)
          nodes(:,:,3) = tmp

        end subroutine transpose_nodes


        subroutine ini_nodes_for_apply_bc_on_fluxes(
     $     x_map, y_map, nodes)

          implicit none

          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          real(rkind) :: dx
          real(rkind) :: dy

          integer(ikind) :: i,j

          dx = 0.5d0
          dy = 0.8d0

          !nodes_tmp
          nodes(7:10,1:10,1) = reshape((/
     $         0.20d0, 0.50d0, 0.46d0, 0.43d0,
     $         0.21d0, 0.52d0, 0.52d0, 0.51d0,
     $         0.20d0, 0.50d0, 0.52d0, 0.53d0,
     $         0.21d0, 0.52d0, 0.48d0, 0.46d0,
     $         0.19d0, 0.48d0, 0.46d0, 0.42d0,
     $         0.17d0, 0.51d0, 0.49d0, 0.44d0,
     $         0.22d0, 0.47d0, 0.45d0, 0.41d0,
     $         0.18d0, 0.49d0, 0.46d0, 0.43d0,
     $         0.17d0, 0.48d0, 0.44d0, 0.39d0,
     $         0.15d0, 0.42d0, 0.46d0, 0.36d0/),
     $         (/4,10/))

          nodes(7:10,1:10,2) = reshape((/
     $         0.59d0,-0.78d0, 0.76d0, 0.72d0,
     $         0.57d0, 0.55d0,-0.79d0, 0.77d0,
     $         0.22d0, 0.77d0, 0.48d0, -0.48d0,
     $        -0.20d0,-0.80d0, 0.46d0, 0.43d0,
     $         0.85d0, 0.52d0, 0.52d0, 0.55d0,
     $        -0.20d0,-0.50d0, 0.52d0, 0.53d0,
     $         0.25d0, 0.52d0, 0.48d0,-0.46d0,
     $        -0.58d0, 0.49d0,-0.46d0, 0.44d0,
     $        -0.57d0, 0.48d0,-0.44d0,-0.49d0,
     $        -0.33d0, 0.15d0,-0.53d0, 0.56d0/),
     $         (/4,10/))

          nodes(7:10,1:10,3) = reshape((/
     $         0.24d0, 0.56d0, 0.53d0, 0.53d0, 
     $         0.20d0, 0.52d0, 0.53d0, 0.50d0, 
     $         0.22d0, 0.54d0, 0.56d0, 0.52d0, 
     $         0.23d0, 0.46d0, 0.42d0, 0.44d0, 
     $         0.32d0, 0.53d0, 0.52d0, 0.53d0, 
     $         0.12d0, 0.51d0, 0.43d0, 0.44d0,
     $         0.30d0, 0.50d0, 0.46d0, 0.43d0,
     $         0.33d0, 0.47d0, 0.45d0, 0.41d0,
     $         0.17d0, 0.48d0, 0.44d0, 0.33d0,
     $         0.15d0, 0.36d9, 0.89d0, 0.75d0/),
     $         (/4,10/))

          nodes(7:10,1:10,4) = reshape((/
     $         0.20d0, 0.46d0, 0.50d0, 0.43d0,
     $         0.21d0, 0.52d0, 0.52d0, 0.51d0,
     $         0.20d0, 0.52d0, 0.50d0, 0.53d0,
     $         0.21d0, 0.48d0, 0.52d0, 0.46d0,
     $         0.19d0, 0.46d0, 0.48d0, 0.42d0,
     $         0.17d0, 0.49d0, 0.51d0, 0.44d0,
     $         0.22d0, 0.45d0, 0.47d0, 0.41d0,
     $         0.18d0, 0.46d0, 0.49d0, 0.43d0,
     $         0.17d0, 0.44d0, 0.48d0, 0.39d0,
     $         0.89d0, 0.26d0, 0.45d0, 0.45d0/),
     $         (/4,10/))

          x_map = (/ ((i-1)*dx,i=1,nx) /)
          y_map = (/ ((j-1)*dy,j=1,ny) /)

        end subroutine ini_nodes_for_apply_bc_on_fluxes


        subroutine ini_nodes_for_apply_bc_on_nodes(
     $     x_map,y_map,nodes_tmp,nodes)

          implicit none

          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          real(rkind) :: velocity_x
          real(rkind) :: velocity_y
          real(rkind) :: temperature
          
          real(rkind) :: dx
          real(rkind) :: dy

          real(rkind) :: grad_x
          real(rkind) :: grad_y

          integer(ikind) :: i,j

          velocity_x  = 0.02d0
          velocity_y  = 0.03d0
          temperature = 0.90d0

          dx = 0.5d0
          dy = 0.8d0

          !nodes_tmp
          nodes_tmp(7:9,1:8,1) = reshape((/
     $         0.20d0,0.50d0 ,0.46d0,
     $         0.21d0,0.52d0 ,0.52d0,
     $         0.20d0,0.50d0 ,0.52d0,
     $         0.21d0,0.52d0 ,0.48d0,
     $         0.19d0,0.48d0 ,0.46d0,
     $         0.17d0,0.51d0 ,0.49d0,
     $         0.22d0,0.47d0 ,0.45d0,
     $         0.18d0,0.49d0 ,0.20d0/),
     $         (/3,8/))

          nodes_tmp(8,9,1) = 0.23d0

          j=1
          do i=7,8
             nodes_tmp(i,j,2) = velocity_x*nodes_tmp(i,j,1)
             nodes_tmp(i,j,3) = velocity_y*nodes_tmp(i,j,1)

             grad_x = (nodes_tmp(i+1,j,1)-nodes_tmp(i-1,j,1))/(2*dx)
             grad_y = (nodes_tmp(i,j+1,1)-nodes_tmp(i,j,1))/(dy)

             nodes_tmp(i,j,4) =
     $            0.5d0*nodes_tmp(i,j,1)*(velocity_x**2+velocity_y**2)+
     $            nodes_tmp(i,j,1)*(8.0d0/3.0d0*cv_r*temperature-3.0d0*nodes_tmp(i,j,1))+
     $            0.5d0/we*(grad_x**2+grad_y**2)
          end do

          do j=2,8
             do i=7,8
                nodes_tmp(i,j,2) = velocity_x*nodes_tmp(i,j,1)
                nodes_tmp(i,j,3) = velocity_y*nodes_tmp(i,j,1)
                
                grad_x = (nodes_tmp(i+1,j,1)-nodes_tmp(i-1,j,1))/(2*dx)
                grad_y = (nodes_tmp(i,j+1,1)-nodes_tmp(i,j-1,1))/(2*dy)

                nodes_tmp(i,j,4) =
     $               0.5d0*nodes_tmp(i,j,1)*(velocity_x**2+velocity_y**2)+
     $               nodes_tmp(i,j,1)*(8.0d0/3.0d0*cv_r*temperature-3.0d0*nodes_tmp(i,j,1))+
     $               0.5d0/we*(grad_x**2+grad_y**2)
             end do
          end do

          !nodes
          nodes(7:8,1:8,1) = reshape((/
     $         0.16d0,0.50d0,
     $         0.22d0,0.53d0,
     $         0.19d0,0.48d0,
     $         0.21d0,0.52d0,
     $         0.19d0,0.48d0,
     $         0.17d0,0.51d0,
     $         0.22d0,0.47d0,
     $         0.18d0,0.49d0/),
     $         (/2,8/))

          do j=1,8
             do i=7,8
                nodes(i,j,2) = velocity_x*nodes(i,j,1)
                nodes(i,j,3) = velocity_y*nodes(i,j,1)
                nodes(i,j,4) =
     $               0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2)+
     $               nodes(i,j,1)*(8.0d0/3.0d0*cv_r*temperature-3.0d0*nodes(i,j,1))+
     $               0.5d0/we*((0.1d0/dx)**2+(0.1/dy)**2)
             end do
          end do

          nodes(8,9,1) = 0.23d0

          x_map = (/ ((i-1)*dx,i=1,nx) /)
          y_map = (/ ((j-1)*dy,j=1,ny) /)

        end subroutine ini_nodes_for_apply_bc_on_nodes


        subroutine test_inputs()

          implicit none

          logical :: test_loc
          logical :: test_validated

          test_validated = .true.

          test_loc = nx.eq.10
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''nx should equal 10'')'
          end if

          test_loc = ny.eq.10
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''ny should equal 10'')'
          end if

          test_loc = ne.eq.4
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''ne should equal 4: DIM physical model'')'
          end if

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
     $         wall_micro_contact_angle,
     $         90.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''wall_contact_angle should equal 90.0d0'')'
          end if

          test_loc = is_real_validated(
     $         wall_maximum_heat_flux,
     $         0.005d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''wall_heat_flux should equal 0.005d0'')'
          end if

          if(.not.test_validated) then
             stop ''
          end if

        end subroutine test_inputs

      end program test_wall_xy
