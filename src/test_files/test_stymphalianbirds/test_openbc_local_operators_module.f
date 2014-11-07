      ! module to test the application of the boundary conditions
      ! on the time derivatives using the edge and corner procedures:
      !
      ! functions tested:
      ! - bc_operators%apply_bc_on_timedev_x_edge
      ! - bc_operators%apply_bc_on_timedev_y_edge
      ! - bc_operators%apply_bc_on_timedev_xy_corner
      !---------------------------------------------------------
      module test_openbc_local_operators_module

        !use ifport

        use bc_operators_class, only :
     $       bc_operators

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       left,
     $       right

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use pmodel_eq_class , only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0

        implicit none

        contains

        subroutine test_openbc_local_operators(detailled)

          implicit none

          logical, intent(in) :: detailled

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind)                      :: dx, dy
          type(pmodel_eq)                  :: p_model
          type(bc_operators)               :: bc_used

          if((nx.ne.7).or.(ny.ne.5)) then
             stop 'the test requires (nx,ny)=(7,5)'
          end if


          !initialization of the nodes
          call initialize_nodes(nodes,dx,dy)
  
  
          !test apply_bc_on_timedev_x_edge
          !test apply_bc_on_timedev_y_edge
          !test apply_bc_on_timedev_xy_corner
          print '(''----------------------------------'')'
          print '(''test_apply_bc_on_timedev_local'')'
          print '(''----------------------------------'')'
          call test_apply_bc_on_timedev_local(
     $         nodes,dx,dy,
     $         p_model,
     $         bc_used,
     $         detailled)
          print '()'

        end subroutine test_openbc_local_operators        

        
        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated        


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          
          dx=0.5
          dy=0.6

          !position
          nodes(1,1,1)=0.5
          nodes(2,1,1)=0.2
          nodes(3,1,1)=1.2
          nodes(4,1,1)=5.0
          nodes(5,1,1)=0.6
          nodes(6,1,1)=-3.6
          nodes(7,1,1)=-6.52

          nodes(1,2,1)=3.0
          nodes(2,2,1)=4.2
          nodes(3,2,1)=11.0
          nodes(4,2,1)=10.6
          nodes(5,2,1)=5.2
          nodes(6,2,1)=1.2
          nodes(7,2,1)=7.89

          nodes(1,3,1)=-14.2
          nodes(2,3,1)=23
          nodes(3,3,1)=9.8
          nodes(4,3,1)=3.4
          nodes(5,3,1)=9.1
          nodes(6,3,1)=6.7
          nodes(7,3,1)=4.12

          nodes(1,4,1)=2.45
          nodes(2,4,1)=0.2
          nodes(3,4,1)=9.0
          nodes(4,4,1)=5.4
          nodes(5,4,1)=-2.3
          nodes(6,4,1)=1.0
          nodes(7,4,1)=-5.62

          nodes(1,5,1)=3.6
          nodes(2,5,1)=0.1
          nodes(3,5,1)=6.3
          nodes(4,5,1)=8.9
          nodes(5,5,1)=-4.23
          nodes(6,5,1)=8.9
          nodes(7,5,1)=8.95


          !velocity_x
          nodes(1,1,2)=7.012
          nodes(2,1,2)=-6.323
          nodes(3,1,2)=3.012
          nodes(4,1,2)=4.5
          nodes(5,1,2)=9.6
          nodes(6,1,2)=9.57
          nodes(7,1,2)=6.012
                    
          nodes(1,2,2)=4.26
          nodes(2,2,2)=4.23
          nodes(3,2,2)=4.5
          nodes(4,2,2)=7.56
          nodes(5,2,2)=7.21
          nodes(6,2,2)=5.62
          nodes(7,2,2)=3.62
                    
          nodes(1,3,2)=0.23
          nodes(2,3,2)=7.23
          nodes(3,3,2)=3.1
          nodes(4,3,2)=8.9
          nodes(5,3,2)=9.3
          nodes(6,3,2)=1.29
          nodes(7,3,2)=7.26
                    
          nodes(1,4,2)=8.23
          nodes(2,4,2)=-3.1
          nodes(3,4,2)=6.03
          nodes(4,4,2)=6.25
          nodes(5,4,2)=5.12
          nodes(6,4,2)=0.36
          nodes(7,4,2)=6.89
                    
          nodes(1,5,2)=3.2
          nodes(2,5,2)=8.12
          nodes(3,5,2)=8.9
          nodes(4,5,2)=4.2
          nodes(5,5,2)=7.8
          nodes(6,5,2)=6.3
          nodes(7,5,2)=1.2          

        end subroutine initialize_nodes


        subroutine test_apply_bc_on_timedev_local(
     $     nodes,dx,dy,
     $     p_model,
     $     bc_used,
     $     detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          type(pmodel_eq)                 , intent(in)    :: p_model
          type(bc_operators)              , intent(inout) :: bc_used
          logical                         , intent(in)    :: detailled

          
          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_validated
          logical           :: loc

          real(rkind), dimension(nx) :: x_map
          real(rkind), dimension(ny) :: y_map
          real(rkind)                :: t

          real(rkind), dimension(ne) :: timedev_local
          logical                    :: side_x
          logical                    :: side_y


          do i=1,nx
             x_map(i) = (i-1)*dx
          end do

          do j=1,ny
             y_map(j) = (j-1)*dy
          end do


          !initialize the y flux at zero
          do k=1,ne
             do j=1,ny+1
                do i=1,nx
                   flux_y(i,j,k)=0.0d0
                end do
             end do
          end do


          !initialize the flux_x randomly
          call srand(10)

          do k=1,ne
             do j=1,ny
                do i=1,nx+1
                   flux_x(i,j,k) = RAND()
                end do
             end do
          end do

          
          !initialize the boundary conditions
          call bc_used%ini(p_model)


          !computation of the time derivatives using open b.c.
          !with the function that has been validated previously
          call bc_used%apply_bc_on_timedev(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         time_dev)


          !computation of N_edge, S_edge, E_edge, W_edge
          !and NE_corner,SE_corner,NW_corner,SW_corner
          !using the local functions
          
          !N_edge: ny-1
          i      = 3
          j      = ny-1
          side_y = right

          timedev_local =
     $         bc_used%apply_bc_on_timedev_y_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_x,
     $         side_y,
     $         gradient_y_y_oneside_R1)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated N_edge(1/2): '',L1)', test_validated

          !N_edge: ny
          i      = 3
          j      = ny
          side_y = right

          timedev_local =
     $         bc_used%apply_bc_on_timedev_y_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_x,
     $         side_y,
     $         gradient_y_y_oneside_R0)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated N_edge(2/2): '',L1)', test_validated


          !S_edge: 1
          i      = 3
          j      = 1
          side_y = left

          timedev_local =
     $         bc_used%apply_bc_on_timedev_y_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_x,
     $         side_y,
     $         gradient_y_y_oneside_L0)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated S_edge(1/2): '',L1)', test_validated

          !S_edge: 2
          i      = 3
          j      = 2
          side_y = left

          timedev_local =
     $         bc_used%apply_bc_on_timedev_y_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_x,
     $         side_y,
     $         gradient_y_y_oneside_L1)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated S_edge(2/2): '',L1)', test_validated


          !E_edge: nx-1
          i      = nx-1
          j      = 3
          side_x = right

          timedev_local =
     $         bc_used%apply_bc_on_timedev_x_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_y,
     $         side_x,
     $         gradient_x_x_oneside_R1)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated E_edge(1/2): '',L1)', test_validated

          !E_edge: nx
          i      = nx
          j      = 3
          side_x = right

          timedev_local =
     $         bc_used%apply_bc_on_timedev_x_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_y,
     $         side_x,
     $         gradient_x_x_oneside_R0)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated E_edge(2/2): '',L1)', test_validated

          
          !W_edge: 1
          i      = 1
          j      = 3
          side_x = left

          timedev_local =
     $         bc_used%apply_bc_on_timedev_x_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_y,
     $         side_x,
     $         gradient_x_x_oneside_L0)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated W_edge(1/2): '',L1)', test_validated

          !E_edge: 2
          i      = 2
          j      = 3
          side_x = left

          timedev_local =
     $         bc_used%apply_bc_on_timedev_x_edge(
     $         p_model, t,
     $         nodes,x_map,y_map,i,j,
     $         flux_y,
     $         side_x,
     $         gradient_x_x_oneside_L1)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(''test_validated W_edge(2/2): '',L1)', test_validated

          print '()'

          !corner SW
          i=1
          j=1
          side_x = left
          side_y = left

          call test_corner(
     $         'test_validated SW_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L0,
     $         detailled)

          i=2
          j=1
          side_x = left
          side_y = left

          call test_corner(
     $         'test_validated SW_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L0,
     $         detailled)

          i=1
          j=2
          side_x = left
          side_y = left

          call test_corner(
     $         'test_validated SW_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L1,
     $         detailled)

          i=2
          j=2
          side_x = left
          side_y = left

          call test_corner(
     $         'test_validated SW_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L1,
     $         detailled)

          !corner NE
          i=nx-1
          j=ny-1
          side_x = right
          side_y = right

          call test_corner(
     $         'test_validated NE_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R1,
     $         detailled)

          i=nx
          j=ny-1
          side_x = right
          side_y = right

          call test_corner(
     $         'test_validated NE_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R1,
     $         detailled)

          i=nx-1
          j=ny
          side_x = right
          side_y = right

          call test_corner(
     $         'test_validated NE_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R0,
     $         detailled)

          i=nx
          j=ny
          side_x = right
          side_y = right

          call test_corner(
     $         'test_validated NE_corner: ',
     $         time_dev,
     $         bc_used,
     $         p_model,t,
     $         nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0,
     $         detailled)

        end subroutine test_apply_bc_on_timedev_local


        subroutine test_corner(
     $     test_name,
     $     time_dev,
     $     bc_used,
     $     p_model,t,
     $     nodes,x_map,y_map,i,j,
     $     side_x,side_y,
     $     gradient_x,
     $     gradient_y,
     $     detailled)

          implicit none

          character(len=26)               , intent(in) :: test_name
          real(rkind), dimension(nx,ny,ne), intent(in) :: time_dev
          type(bc_operators)              , intent(in) :: bc_used
          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: t
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          logical                         , intent(in) :: side_x
          logical                         , intent(in) :: side_y
          procedure(gradient_x_proc)                   :: gradient_x
          procedure(gradient_y_proc)                   :: gradient_y
          logical                         , intent(in) :: detailled

          real(rkind), dimension(ne) :: timedev_local
          logical                    :: test_validated
          integer                    :: k
          logical                    :: loc
          
          timedev_local =
     $         bc_used%apply_bc_on_timedev_xy_corner(
     $         p_model, t,
     $         nodes, x_map, y_map, i,j,
     $         side_x,
     $         side_y,
     $         gradient_x,
     $         gradient_y)

          test_validated = .true.

          do k=1,ne
             loc = is_test_validated(
     $            timedev_local(k), time_dev(i,j,k), detailled)
             test_validated = test_validated.and.loc
          end do

          print '(A26,L1)', test_name, test_validated

        end subroutine test_corner

      end module test_openbc_local_operators_module
