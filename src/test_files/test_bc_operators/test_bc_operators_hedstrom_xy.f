      program test_bc_operators_hedstrom_xy

        use bc_operators_hedstrom_xy_class, only :
     $       bc_operators_hedstrom_xy

        use check_data_module, only :
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SW_corner_type,
     $       SE_corner_type

        use parameters_input, only :
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice,
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        call test_inputs()

        test_loc = test_apply_bc_on_timedev(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated

        contains

        function test_apply_bc_on_timedev(detailled)
     $       result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators_hedstrom_xy)     :: bc_operators_used
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind), dimension(nx,ny,ne)   :: nodes
          type(pmodel_eq)                    :: p_model
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(nx,ny,ne)   :: timedev
          real(rkind), dimension(nx,ny,ne)   :: timedev_test

          logical :: test_loc

          test_validated = .true.


          ! initialize the x_map,y_map,nodes
          call get_nodes_for_test_apply_bc_on_timedev(
     $         x_map,
     $         y_map,
     $         nodes)


          ! compute the time dev for the E edge
          call bc_operators_used%apply_bc_on_timedev(
     $         [E_edge_type,9,3,8],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          timedev_test = timedev


          ! compare with the time dev for the W edge
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [W_edge_type,1,3,8],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(1:2,3:8,:),
     $         timedev_test(1:2,3:8,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test W_edge failed'')'
          end if


          ! compare with the time dev for the S edge
          call transpose_nodes(nodes)
          call transpose_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [S_edge_type,3,1,8],
     $         0.0,y_map,x_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(3:8,1:2,:),
     $         timedev_test(3:8,1:2,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test S_edge failed'')'
          end if


          ! compare with the time dev for the N edge
          call reflect_y_nodes(nodes)
          call reflect_y_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [N_edge_type,3,9,8],
     $         0.0,y_map,x_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(3:8,9:10,:),
     $         timedev_test(3:8,9:10,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test N_edge failed'')'
          end if


          ! compute the NW corner
          call bc_operators_used%apply_bc_on_timedev(
     $         [NW_corner_type,1,9,0],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          timedev_test = timedev


          ! compute the NE corner
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [NE_corner_type,9,9,0],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(9:10,9:10,:),
     $         timedev_test(9:10,9:10,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test NE_corner failed'')'
          end if


          ! compute the SE corner
          call reflect_y_nodes(nodes)
          call reflect_y_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [SE_corner_type,9,1,0],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(9:10,1:2,:),
     $         timedev_test(9:10,1:2,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test SE_corner failed'')'
          end if


          ! compute the SW corner
          call reflect_x_nodes(nodes)
          call reflect_x_nodes(timedev_test)

          call bc_operators_used%apply_bc_on_timedev(
     $         [SW_corner_type,1,1,0],
     $         0.0,x_map,y_map,nodes,
     $         p_model,flux_x,flux_y,
     $         timedev)

          test_loc = is_real_matrix3D_validated(
     $         timedev(1:2,1:2,:),
     $         timedev_test(1:2,1:2,:),
     $         detailled)

          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''test SW_corner failed'')'
          end if
          

        end function test_apply_bc_on_timedev


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


        subroutine get_nodes_for_test_apply_bc_on_timedev(
     $     x_map,
     $     y_map,
     $     nodes)

          implicit none

          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          integer     :: i
          real(rkind) :: dx
          real(rkind) :: dy

          dx = 0.1d0
          dy = 0.2d0

          x_map = (/ ((i-1)*dx,i=1,nx) /)
          y_map = (/ ((i-1)*dy,i=1,ny) /)


          nodes(8:10,1:10,1) = reshape((/
     $         1.46d0,  1.27d0,  1.47d0,
     $         1.23d0,  1.19d0,  1.35d0,
     $         1.26d0,  1.22d0,  1.25d0,
     $         1.33d0,  1.23d0,  1.46d0,
     $         1.44d0,  1.21d0,  1.46d0,
     $         1.27d0,  1.34d0,  1.46d0,
     $         1.23d0,  1.37d0,  1.45d0,
     $         1.44d0,  1.37d0,  1.13d0,
     $         1.46d0,  1.45d0,  1.22d0,
     $         1.41d0,  1.24d0,  1.30d0/),
     $         (/3,10/))

          nodes(8:10,1:10,2) = reshape((/
     $         0.143d0, 0.143d0, 0.125d0, 
     $         0.146d0, 0.135d0, 0.137d0, 
     $         0.142d0, 0.133d0, 0.123d0,
     $         0.135d0, 0.145d0, 0.139d0, 
     $         0.138d0, 0.145d0, 0.140d0, 
     $         0.142d0, 0.126d0, 0.141d0, 
     $         0.132d0, 0.122d0, 0.146d0, 
     $         0.141d0, 0.144d0, 0.145d0, 
     $         0.146d0, 0.143d0, 0.124d0,
     $         0.141d0, 0.147d0, 0.131d0/),
     $         (/3,10/))

          nodes(8:10,1:10,3) = reshape((/
     $         0.072d0, 0.050d0, 0.045d0, 
     $         0.012d0, 0.034d0, 0.020d0,
     $         0.002d0, 0.022d0, 0.033d0,
     $         0.013d0, 0.013d0, 0.002d0, 
     $         0.045d0, 0.042d0, 0.016d0, 
     $         0.050d0, 0.035d0, 0.004d0, 
     $         0.032d0, 0.012d0, 0.092d0,
     $         0.054d0, 0.005d0, 0.052d0,
     $         0.011d0, 0.012d0, 0.061d0,
     $         0.015d0, 0.013d0, 0.071d0/),
     $         (/3,10/))

          nodes(8:10,1:10,4) = reshape((/
     $         4.86d0,  4.86d0,  4.63d0,  
     $         4.83d0,  4.90d0,  4.78d0,
     $         4.45d0,  4.86d0,  4.89d0,  
     $         4.45d0,  4.86d0,  4.87d0,
     $         4.89d0,  4.86d0,  4.74d0,  
     $         4.81d0,  4.56d0,  4.78d0,  
     $         4.76d0,  4.83d0,  4.85d0,  
     $         4.86d0,  4.82d0,  4.89d0,  
     $         4.52d0,  4.83d0,  4.91d0,  
     $         4.81d0,  4.82d0,  4.83d0/),
     $         (/3,10/))

        end subroutine get_nodes_for_test_apply_bc_on_timedev


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

          if(.not.test_validated) then
             stop ''
          end if

        end subroutine test_inputs

      end program test_bc_operators_hedstrom_xy
