      ! module to test the application of the boundary conditions
      ! on the time derivatives using the edge and corner procedures:
      !
      ! functions tested:
      ! - bc_operators%apply_bc_on_timedev_x_edge
      ! - bc_operators%apply_bc_on_timedev_y_edge
      ! - bc_operators%apply_bc_on_timedev_xy_corner
      !---------------------------------------------------------
      module test_openbc_local_operators_module

        use ifport

        use bc_operators_class, only :
     $       bc_operators

        use bf_layer_bc_procedure_module, only :
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       N_edge_type

c$$$        use interface_primary, only :
c$$$     $       gradient_x_proc,
c$$$     $       gradient_y_proc
c$$$
c$$$        use openbc_operators_module, only :
c$$$     $       incoming_left,
c$$$     $       incoming_right
c$$$
c$$$        use parameters_constant, only :
c$$$     $       left,
c$$$     $       right

        use parameters_constant, only :
     $       peak

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size,
     $       ic_choice

        use pmodel_eq_class , only :
     $       pmodel_eq

c$$$        use sd_operators_fd_module, only :
c$$$     $       gradient_x_x_oneside_L0,
c$$$     $       gradient_x_x_oneside_L1,
c$$$     $       gradient_x_x_oneside_R1,
c$$$     $       gradient_x_x_oneside_R0,
c$$$     $       gradient_y_y_oneside_L0,
c$$$     $       gradient_y_y_oneside_L1,
c$$$     $       gradient_y_y_oneside_R1,
c$$$     $       gradient_y_y_oneside_R0

        implicit none

        contains

        subroutine test_openbc_local_operators(detailled)

          implicit none

          logical, intent(in) :: detailled


          type(bc_operators)  :: bc_used
          logical             :: test_validated

          if((nx.ne.10).or.(ny.ne.10)) then
             stop 'the test requires (nx,ny)=(10,10)'
          end if

          if(ic_choice.ne.peak) then
             stop 'the test requires NS equations + peak ic'
          end if


          !test apply_bc_on_timedev_x_edge
          !test apply_bc_on_timedev_y_edge
          !test apply_bc_on_timedev_xy_corner
          print '(''----------------------------------'')'
          print '(''test_apply_bc_on_timedev_local'')'
          print '(''----------------------------------'')'
          test_validated = test_apply_bc_on_timedev_local(
     $         bc_used,
     $         detailled)
          print '(''test_validated: '',L1)', test_validated
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


        

        function test_apply_bc_on_timedev_local(
     $     bc_used,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(bc_operators), intent(inout) :: bc_used
          logical           , intent(in)    :: detailled
          logical                           :: test_validated

          real(rkind)                        :: t
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          type(pmodel_eq)                    :: p_model

          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx,ny,ne)   :: time_dev_nopt
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)                     :: i,j
          integer                            :: k        

          integer, dimension(:,:), allocatable :: bc_sections
          logical                              :: test_edge_N
          logical                              :: test_edge_S
          logical                              :: test_edge_E
          logical                              :: test_edge_W
          logical                              :: test_loc


          !the x_map and y_map are initialized such that the
          !peak initial conditions for the NS equations cover
          !the computational domain entirely
          dx = 0.2/nx
          dy = 0.2/ny

          do i=1,nx
             x_map(i) = -0.1 + (i-1)*dx
          end do

          do j=1,ny
             y_map(j) = -0.1 + (j-1)*dy
          end do


          !initialize the nodes
          call p_model%apply_ic(nodes,x_map,y_map)

          print '(''mass: '')'
          print '(''------'')'
          do j=1, ny
             print '(10F8.3)', nodes(:,j,1)
          end do
          print '()'

          print '(''momentum_x: '')'
          print '(''------------'')'
          do j=1, ny
             print '(10F8.3)', nodes(:,j,2)
          end do
          print '()'

          print '(''momentum_y: '')'
          print '(''------------'')'
          do j=1, ny
             print '(10F8.3)', nodes(:,j,3)
          end do
          print '()'

          print '(''total_energy: '')'
          print '(''--------------'')'
          do j=1, ny
             print '(10F8.3)', nodes(:,j,4)
          end do
          print '()'


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

          !1) definition of the bc_sections
          allocate(bc_sections(4,8))

          bc_sections(:,1) = [SW_corner_type,1,1,0]
          bc_sections(:,2) = [S_edge_type,bc_size+1,1,nx-bc_size]
          bc_sections(:,3) = [SE_corner_type,nx-bc_size+1,1,0]
          bc_sections(:,4) = [W_edge_type,1,bc_size+1,ny-bc_size]
          bc_sections(:,5) = [E_edge_type,nx-bc_size+1,bc_size+1,ny-bc_size]
          bc_sections(:,6) = [NW_corner_type,1,ny-bc_size+1,0]
          bc_sections(:,7) = [N_edge_type,bc_size+1,ny-bc_size+1,nx-bc_size]
          bc_sections(:,8) = [NE_corner_type,nx-bc_size+1,ny-bc_size+1,0]


          !2) computation of the time derivatives using the
          !   non-optimized function for the interior
          call bc_used%apply_bc_on_timedev_nopt(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x, flux_y,
     $         time_dev_nopt,
     $         bc_sections)

          if(detailled) then
             print '(''timedev mass: '')'
             print '(''--------------'')'
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev(:,j,1)
             end do
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev_nopt(:,j,1)
             end do
             print '()'

             print '(''timedev momentum_x: '')'
             print '(''--------------------'')'
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev(:,j,2)
             end do
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev_nopt(:,j,2)
             end do
             print '()'

             print '(''timedev momentum_y: '')'
             print '(''--------------------'')'
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev(:,j,3)
             end do
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev_nopt(:,j,3)
             end do
             print '()'

             print '(''timedev total_energy: '')'
             print '(''----------------------'')'
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev(:,j,4)
             end do
             print '()'
             do j=1, ny
                print '(10F8.3)', time_dev_nopt(:,j,4)
             end do
             print '()'
          end if


          !3) comparison of the results from the two functions
          !comparison S_edge
          test_validated = .true.
          do k=1,ne
             do j=1,bc_size
                do i=1, nx
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k),
     $                  time_dev_nopt(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do
          test_edge_S = test_validated
          if(.not.detailled) then
             print '(''test (SW_corner)+(edge_S)+(SE_corner): '',L1)',
     $            test_edge_S
          end if

          !comparison W_edge
          test_validated = .true.
          do k=1,ne
             do j=bc_size+1,ny-bc_size
                do i=1, bc_size
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k),
     $                  time_dev_nopt(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do
          test_edge_W = test_validated
          if(.not.detailled) then
             print '(''test_edge_W: '',L1)', test_edge_W
          end if

          !comparison E_edge
          test_validated = .true.
          do k=1,ne
             do j=bc_size+1,ny-bc_size
                do i=nx-bc_size+1, nx
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k),
     $                  time_dev_nopt(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do
          test_edge_E = test_validated
          if(.not.detailled) then
             print '(''test_edge_E: '',L1)', test_edge_E
          end if

          !comparison N_edge
          test_validated = .true.
          do k=1,ne
             do j=ny-bc_size+1,ny
                do i=1, nx
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k),
     $                  time_dev_nopt(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do
          test_edge_N = test_validated
          if(.not.detailled) then
             print '(''test (NW_corner)+(edge_N)+(NE_corner): '',L1)',
     $            test_edge_N
          end if


c$$$          !N_edge: ny-1
c$$$          i      = 3
c$$$          j      = ny-1
c$$$          side_y = right
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_y_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_x,
c$$$     $         side_y,
c$$$     $         gradient_y_y_oneside_R1)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated N_edge(1/2): '',L1)', test_validated
c$$$
c$$$          !N_edge: ny
c$$$          i      = 3
c$$$          j      = ny
c$$$          side_y = right
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_y_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_x,
c$$$     $         side_y,
c$$$     $         gradient_y_y_oneside_R0)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated N_edge(2/2): '',L1)', test_validated
c$$$
c$$$
c$$$          !S_edge: 1
c$$$          i      = 3
c$$$          j      = 1
c$$$          side_y = left
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_y_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_x,
c$$$     $         side_y,
c$$$     $         gradient_y_y_oneside_L0)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated S_edge(1/2): '',L1)', test_validated
c$$$
c$$$          !S_edge: 2
c$$$          i      = 3
c$$$          j      = 2
c$$$          side_y = left
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_y_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_x,
c$$$     $         side_y,
c$$$     $         gradient_y_y_oneside_L1)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated S_edge(2/2): '',L1)', test_validated
c$$$
c$$$
c$$$          !E_edge: nx-1
c$$$          i      = nx-1
c$$$          j      = 3
c$$$          side_x = right
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_x_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_y,
c$$$     $         side_x,
c$$$     $         gradient_x_x_oneside_R1)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated E_edge(1/2): '',L1)', test_validated
c$$$
c$$$          !E_edge: nx
c$$$          i      = nx
c$$$          j      = 3
c$$$          side_x = right
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_x_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_y,
c$$$     $         side_x,
c$$$     $         gradient_x_x_oneside_R0)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated E_edge(2/2): '',L1)', test_validated
c$$$
c$$$          
c$$$          !W_edge: 1
c$$$          i      = 1
c$$$          j      = 3
c$$$          side_x = left
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_x_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_y,
c$$$     $         side_x,
c$$$     $         gradient_x_x_oneside_L0)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated W_edge(1/2): '',L1)', test_validated
c$$$
c$$$          !E_edge: 2
c$$$          i      = 2
c$$$          j      = 3
c$$$          side_x = left
c$$$
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_x_edge(
c$$$     $         p_model, t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         flux_y,
c$$$     $         side_x,
c$$$     $         gradient_x_x_oneside_L1)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(''test_validated W_edge(2/2): '',L1)', test_validated
c$$$
c$$$          print '()'
c$$$
c$$$          !corner SW
c$$$          i=1
c$$$          j=1
c$$$          side_x = left
c$$$          side_y = left
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated SW_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_L0,
c$$$     $         gradient_y_y_oneside_L0,
c$$$     $         detailled)
c$$$
c$$$          i=2
c$$$          j=1
c$$$          side_x = left
c$$$          side_y = left
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated SW_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_L1,
c$$$     $         gradient_y_y_oneside_L0,
c$$$     $         detailled)
c$$$
c$$$          i=1
c$$$          j=2
c$$$          side_x = left
c$$$          side_y = left
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated SW_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_L0,
c$$$     $         gradient_y_y_oneside_L1,
c$$$     $         detailled)
c$$$
c$$$          i=2
c$$$          j=2
c$$$          side_x = left
c$$$          side_y = left
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated SW_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_L1,
c$$$     $         gradient_y_y_oneside_L1,
c$$$     $         detailled)
c$$$
c$$$          !corner NE
c$$$          i=nx-1
c$$$          j=ny-1
c$$$          side_x = right
c$$$          side_y = right
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated NE_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_R1,
c$$$     $         gradient_y_y_oneside_R1,
c$$$     $         detailled)
c$$$
c$$$          i=nx
c$$$          j=ny-1
c$$$          side_x = right
c$$$          side_y = right
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated NE_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_R0,
c$$$     $         gradient_y_y_oneside_R1,
c$$$     $         detailled)
c$$$
c$$$          i=nx-1
c$$$          j=ny
c$$$          side_x = right
c$$$          side_y = right
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated NE_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_R1,
c$$$     $         gradient_y_y_oneside_R0,
c$$$     $         detailled)
c$$$
c$$$          i=nx
c$$$          j=ny
c$$$          side_x = right
c$$$          side_y = right
c$$$
c$$$          call test_corner(
c$$$     $         'test_validated NE_corner: ',
c$$$     $         time_dev,
c$$$     $         bc_used,
c$$$     $         p_model,t,
c$$$     $         nodes,x_map,y_map,i,j,
c$$$     $         side_x,side_y,
c$$$     $         gradient_x_x_oneside_R0,
c$$$     $         gradient_y_y_oneside_R0,
c$$$     $         detailled)

        end function test_apply_bc_on_timedev_local


c$$$        subroutine test_corner(
c$$$     $     test_name,
c$$$     $     time_dev,
c$$$     $     bc_used,
c$$$     $     p_model,t,
c$$$     $     nodes,x_map,y_map,i,j,
c$$$     $     side_x,side_y,
c$$$     $     gradient_x,
c$$$     $     gradient_y,
c$$$     $     detailled)
c$$$
c$$$          implicit none
c$$$
c$$$          character(len=26)               , intent(in) :: test_name
c$$$          real(rkind), dimension(nx,ny,ne), intent(in) :: time_dev
c$$$          type(bc_operators)              , intent(in) :: bc_used
c$$$          type(pmodel_eq)                 , intent(in) :: p_model
c$$$          real(rkind)                     , intent(in) :: t
c$$$          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
c$$$          real(rkind), dimension(nx)      , intent(in) :: x_map
c$$$          real(rkind), dimension(ny)      , intent(in) :: y_map
c$$$          integer(ikind)                  , intent(in) :: i
c$$$          integer(ikind)                  , intent(in) :: j
c$$$          logical                         , intent(in) :: side_x
c$$$          logical                         , intent(in) :: side_y
c$$$          procedure(gradient_x_proc)                   :: gradient_x
c$$$          procedure(gradient_y_proc)                   :: gradient_y
c$$$          logical                         , intent(in) :: detailled
c$$$
c$$$          real(rkind), dimension(ne) :: timedev_local
c$$$          logical                    :: test_validated
c$$$          integer                    :: k
c$$$          logical                    :: loc
c$$$          
c$$$          timedev_local =
c$$$     $         bc_used%apply_bc_on_timedev_xy_corner(
c$$$     $         p_model, t,
c$$$     $         nodes, x_map, y_map, i,j,
c$$$     $         side_x,
c$$$     $         side_y,
c$$$     $         gradient_x,
c$$$     $         gradient_y)
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          do k=1,ne
c$$$             loc = is_test_validated(
c$$$     $            timedev_local(k), time_dev(i,j,k), detailled)
c$$$             test_validated = test_validated.and.loc
c$$$          end do
c$$$
c$$$          print '(A26,L1)', test_name, test_validated
c$$$
c$$$        end subroutine test_corner

c$$$        subroutine initialize_nodes(nodes,dx,dy)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
c$$$          real(rkind)                     , intent(out) :: dx
c$$$          real(rkind)                     , intent(out) :: dy
c$$$          
c$$$          dx=0.5
c$$$          dy=0.6
c$$$
c$$$          !position
c$$$          nodes(1,1,1)=0.5
c$$$          nodes(2,1,1)=0.2
c$$$          nodes(3,1,1)=1.2
c$$$          nodes(4,1,1)=5.0
c$$$          nodes(5,1,1)=0.6
c$$$          nodes(6,1,1)=-3.6
c$$$          nodes(7,1,1)=-6.52
c$$$
c$$$          nodes(1,2,1)=3.0
c$$$          nodes(2,2,1)=4.2
c$$$          nodes(3,2,1)=11.0
c$$$          nodes(4,2,1)=10.6
c$$$          nodes(5,2,1)=5.2
c$$$          nodes(6,2,1)=1.2
c$$$          nodes(7,2,1)=7.89
c$$$
c$$$          nodes(1,3,1)=-14.2
c$$$          nodes(2,3,1)=23
c$$$          nodes(3,3,1)=9.8
c$$$          nodes(4,3,1)=3.4
c$$$          nodes(5,3,1)=9.1
c$$$          nodes(6,3,1)=6.7
c$$$          nodes(7,3,1)=4.12
c$$$
c$$$          nodes(1,4,1)=2.45
c$$$          nodes(2,4,1)=0.2
c$$$          nodes(3,4,1)=9.0
c$$$          nodes(4,4,1)=5.4
c$$$          nodes(5,4,1)=-2.3
c$$$          nodes(6,4,1)=1.0
c$$$          nodes(7,4,1)=-5.62
c$$$
c$$$          nodes(1,5,1)=3.6
c$$$          nodes(2,5,1)=0.1
c$$$          nodes(3,5,1)=6.3
c$$$          nodes(4,5,1)=8.9
c$$$          nodes(5,5,1)=-4.23
c$$$          nodes(6,5,1)=8.9
c$$$          nodes(7,5,1)=8.95
c$$$
c$$$
c$$$          !velocity_x
c$$$          nodes(1,1,2)=7.012
c$$$          nodes(2,1,2)=-6.323
c$$$          nodes(3,1,2)=3.012
c$$$          nodes(4,1,2)=4.5
c$$$          nodes(5,1,2)=9.6
c$$$          nodes(6,1,2)=9.57
c$$$          nodes(7,1,2)=6.012
c$$$                    
c$$$          nodes(1,2,2)=4.26
c$$$          nodes(2,2,2)=4.23
c$$$          nodes(3,2,2)=4.5
c$$$          nodes(4,2,2)=7.56
c$$$          nodes(5,2,2)=7.21
c$$$          nodes(6,2,2)=5.62
c$$$          nodes(7,2,2)=3.62
c$$$                    
c$$$          nodes(1,3,2)=0.23
c$$$          nodes(2,3,2)=7.23
c$$$          nodes(3,3,2)=3.1
c$$$          nodes(4,3,2)=8.9
c$$$          nodes(5,3,2)=9.3
c$$$          nodes(6,3,2)=1.29
c$$$          nodes(7,3,2)=7.26
c$$$                    
c$$$          nodes(1,4,2)=8.23
c$$$          nodes(2,4,2)=-3.1
c$$$          nodes(3,4,2)=6.03
c$$$          nodes(4,4,2)=6.25
c$$$          nodes(5,4,2)=5.12
c$$$          nodes(6,4,2)=0.36
c$$$          nodes(7,4,2)=6.89
c$$$                    
c$$$          nodes(1,5,2)=3.2
c$$$          nodes(2,5,2)=8.12
c$$$          nodes(3,5,2)=8.9
c$$$          nodes(4,5,2)=4.2
c$$$          nodes(5,5,2)=7.8
c$$$          nodes(6,5,2)=6.3
c$$$          nodes(7,5,2)=1.2          
c$$$
c$$$        end subroutine initialize_nodes

      end module test_openbc_local_operators_module
