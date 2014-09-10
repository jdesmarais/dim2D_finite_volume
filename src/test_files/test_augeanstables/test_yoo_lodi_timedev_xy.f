      program test_yoo_lodi_timedev_xy

        use lodi_corner_inflow_inflow_class, only :
     $     lodi_corner_inflow_inflow

        use lodi_corner_inflow_outflow_class, only :
     $       lodi_corner_inflow_outflow

        use lodi_corner_outflow_outflow_class, only :
     $       lodi_corner_outflow_outflow

        use lodi_edge_inflow_class, only :
     $       lodi_edge_inflow

        use lodi_edge_outflow_class, only :
     $       lodi_edge_outflow

        use lodi_timedev_xy_module, only :
     $       get_flow_config,
     $       compute_timedev_x_edge,
     $       compute_timedev_corner

        use ns2d_parameters, only :
     $       gamma,
     $       mach_infty

        use parameters_constant, only :
     $       left, right,
     $       inflow_type, outflow_type,
     $       ask_flow, always_inflow, always_outflow,
     $       vector_x,x_direction

        use parameters_input, only :
     $       sigma_P,
     $       flow_direction,
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
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

        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind), dimension(nx)         :: x_map
        real(rkind), dimension(ny)         :: y_map
        real(rkind)                        :: t
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        type(pmodel_eq)                    :: p_model

        logical :: test_validated
        logical :: detailled
        

        !check the parameters
        if((nx.ne.5).or.(ny.ne.5).or.(ne.ne.4)) then
           print '(''test designed for:'')'
           print '(''nx=5'')'
           print '(''ny=5'')'
           print '(''pm_model=ns2d'')'
           stop 'change inputs'
        end if

        
        detailled = .false.
        if(
     $       (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(mach_infty,0.2d0,detailled)).or.
     $       (.not.is_test_validated(sigma_P,0.25d0,detailled)).or.
     $       (.not.(flow_direction.eq.x_direction))) then

           print '(''the test requires: '')'
           print '(''gamma=5/3'')'
           print '(''mach_infty=0.2'')'
           print '(''sigma_P=0.25'')'
           print '(''flow_direction=x-direction'')'
           stop ''

        end if

        
        !test get_flow_config
        print '(''test get_flow_config'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        test_validated = test_get_flow_config(detailled)

        print '()'


        !test compute_timedev_x_edge
        print '(''test compute_timedev_x_edge'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

        detailled = .false.

        test_validated = test_compute_timedev_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       detailled)

        print '()'


        !test compute_timedev_corner
        print '(''test compute_timedev_corner'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

        detailled = .false.

        test_validated = test_compute_timedev_corner(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       detailled)

        print '()'
        
        

        contains

        !initialize the nodes
        subroutine initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

          implicit none

          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy

          integer, dimension(ne) :: var_type
          integer(ikind) :: i,j
          integer        :: k


          !fill the nodes (i in [1,3])x(j in [1,5])
          !mass----------------
          nodes(1,1,1) =  2.3d0
          nodes(2,1,1) =  1.02d0
          nodes(3,1,1) =  3.2d0
                          
          nodes(1,2,1) =  1.2d0
          nodes(2,2,1) =  8.6d0
          nodes(3,2,1) =  6.13d0
                          
          nodes(1,3,1) =  0.23d0
          nodes(2,3,1) =  4.5d0
          nodes(3,3,1) =  7.13d0
                          
          nodes(1,4,1) =  8.5d0
          nodes(2,4,1) =  7.8d0
          nodes(3,4,1) =  1.5d0
                          
          nodes(1,5,1) =  0.2d0
          nodes(2,5,1) =  3.6d0
          nodes(3,5,1) =  9.23d0

          !momentum-x----------
          nodes(1,1,2) = -6.045d0
          nodes(2,1,2) =  8.125d0
          nodes(3,1,2) =  0.0d0
                    
          nodes(1,2,2) = -6.3d0
          nodes(2,2,2) =  7.98d0
          nodes(3,2,2) =  0.0d0
                    
          nodes(1,3,2) = -0.15d0
          nodes(2,3,2) = -6.213d0
          nodes(3,3,2) =  0.0d0
                    
          nodes(1,4,2) =  8.23d0
          nodes(2,4,2) =  3.012d0
          nodes(3,4,2) =  0.0d0
                    
          nodes(1,5,2) = -1.23d0
          nodes(2,5,2) =  7.8d0
          nodes(3,5,2) =  0.0d0

          !momentum-y----------
          nodes(1,1,3) =  2.01d0
          nodes(2,1,3) =  3.25d0
          nodes(3,1,3) =  6.2d0
                    
          nodes(1,2,3) =  7.135d0
          nodes(2,2,3) = -2.01d0
          nodes(3,2,3) =  3.06d0
                    
          nodes(1,3,3) =  9.46d0
          nodes(2,3,3) =  9.16d0
          nodes(3,3,3) =  4.12d0
                    
          nodes(1,4,3) =  2.13d0
          nodes(2,4,3) = -2.15d0
          nodes(3,4,3) = -3.25d0
                    
          nodes(1,5,3) =  6.1023d0
          nodes(2,5,3) =  5.23d0
          nodes(3,5,3) =  1.12d0

          !total_energy--------
          nodes(1,1,4) =  20.1d0
          nodes(2,1,4) =  895.26d0
          nodes(3,1,4) =  961.23d0

          nodes(1,2,4) =  78.256d0
          nodes(2,2,4) =  8.45d0
          nodes(3,2,4) =  7.4d0

          nodes(1,3,4) =  256.12d0
          nodes(2,3,4) =  163.48d0
          nodes(3,3,4) =  9.56d0
                    
          nodes(1,4,4) =  56.12d0
          nodes(2,4,4) =  7.89d0
          nodes(3,4,4) =  629.12d0
                    
          nodes(1,5,4) =  102.3d0
          nodes(2,5,4) =  231.02d0
          nodes(3,5,4) =  7.123d0

          
          !fill the nodes (i in [4,5])x(j in [1,5]) by using
          !the symmetry along the x-axis
          var_type = p_model%get_var_type()

          do k=1, ne
             do j=1,5
                do i=1,2
                   if(var_type(k).eq.vector_x) then
                      nodes(6-i,j,k) = - nodes(i,j,k)
                   else
                      nodes(6-i,j,k) =   nodes(i,j,k)
                   end if
                end do
             end do
          end do


          !set to zero the nodes (i in [3])x(j in [1,5]) for
          !variables of type vector_x
          i=3
          do k=1,ne
             if(var_type(k).eq.vector_x) then
                do j=1,5
                   nodes(i,j,k)=0.0
                end do
             end if
          end do


          !initialize dx
          dx = 0.4d0

          !intiialize dy
          dy = 0.6d0

          !initialize the x_map
          do i=1,5
             x_map(i) = (i-1)*dx
          end do

          !initialize the y_map
          do i=1,5
             y_map(i) = (i-1)*dy
          end do
          
        end subroutine initialize_nodes


        !check the data
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


        !check vector data
        function is_vector_validated(matrix,test_data,detailled)
     $     result(test_validated)
        
          implicit none

          real(rkind), dimension(:), intent(in) :: matrix
          real(rkind), dimension(:), intent(in) :: test_data
          logical                  , intent(in) :: detailled
          logical                               :: test_validated

          integer :: i
          logical :: test_loc

          test_validated = .true.

          do i=1, size(matrix,1)
                
             test_loc = is_test_validated(
     $            matrix(i),
     $            test_data(i),
     $            detailled)
             
             test_validated = test_validated.and.test_loc
             if(detailled) then
                print '(''('',I2,''):'',L2)',i,test_loc
             end if
                
          end do

        end function is_vector_validated


        function is_matrix_validated(matrix,test_data,detailled)
     $     result(test_validated)
        
          implicit none

          real(rkind), dimension(:,:), intent(in) :: matrix
          real(rkind), dimension(:,:), intent(in) :: test_data
          logical                    , intent(in) :: detailled
          logical                                 :: test_validated

          integer :: i,j
          logical :: test_loc

          test_validated = .true.

          do j=1, size(matrix,2)
             do i=1, size(matrix,1)
                
                test_loc = is_test_validated(
     $               matrix(i,j),
     $               test_data(i,j),
     $               detailled)

                test_validated = test_validated.and.test_loc
                if(detailled) then
                   print '(''('',2I2,''):'',L2)',i,j,test_loc
                end if

             end do
          end do

        end function is_matrix_validated


        !get_flow_config
        function test_get_flow_config(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer    , parameter           :: nb_tests = 8
          real(rkind), dimension(nb_tests) :: vector
          logical    , dimension(nb_tests) :: side
          integer    , dimension(nb_tests) :: usr_config
          logical    , dimension(nb_tests) :: test_data

          integer :: i
          logical :: test_loc


          test_validated = .true.


          !fill the test data
          vector(1)     = 1.0d0
          side(1)       = left
          usr_config(1) = ask_flow
          test_data(1)  = inflow_type

          vector(2)     = -1.0d0
          side(2)       = left
          usr_config(2) = ask_flow
          test_data(2)  = outflow_type

          vector(3)     = 1.0d0
          side(3)       = right
          usr_config(3) = ask_flow
          test_data(3)  = outflow_type

          vector(4)     = -1.0d0
          side(4)       = right
          usr_config(4) = ask_flow
          test_data(4)  = inflow_type

          vector(5)     = 1.0d0
          side(5)       = left
          usr_config(5) = always_outflow
          test_data(5)  = outflow_type

          vector(6)     = -1.0d0
          side(6)       = left
          usr_config(6) = always_inflow
          test_data(6)  = inflow_type

          vector(7)     = 1.0d0
          side(7)       = right
          usr_config(7) = always_inflow
          test_data(7)  = inflow_type

          vector(8)     = -1.0d0
          side(8)       = right
          usr_config(8) = always_outflow
          test_data(8)  = outflow_type


          !test the data
          do i=1, nb_tests

             test_loc = (
     $            get_flow_config(
     $            vector(i),
     $            side(i),
     $            usr_config(i)).eqv.test_data(i))

             test_validated = test_validated.and.test_loc

             if(detailled) then
                print '(''test('',I1,''): '',L1)', i, test_validated
             end if

          end do

          if(.not.detailled) then
             print '(''test_get_flow_config: '',L1)', test_validated
          end if

        end function test_get_flow_config


        !compute_timedev_x_edge
        function test_compute_timedev_edge(
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: t
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          real(rkind), dimension(ne)         :: transverse_lodi
          real(rkind), dimension(ne)         :: viscous_lodi
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(4,ne)       :: test_data
          type(lodi_edge_inflow)             :: inflow_bc
          type(lodi_edge_outflow)            :: outflow_bc

          integer :: k
          integer :: flow_x_user


          !initialize the flux_y, transverse_lodi and viscous_lodi
          do k=1, size(transverse_lodi,1)
             transverse_lodi(k) = 0.1d0
             viscous_lodi(k)    = 0.2d0
          end do

          do k=1, ne
             flux_y(1,3,k) = 2.3d0
             flux_y(1,4,k) = 6.25d0

             flux_y(2,3,k) = 7.89d0
             flux_y(2,4,k) = 9.45d0

             flux_y(3,3,k) =-4.1263d0
             flux_y(3,4,k) =-3.6d0

             flux_y(4,3,k) =  3.6d0
             flux_y(4,4,k) = -5.6d0

             flux_y(5,3,k) = 1.23d0
             flux_y(5,4,k) = -8.9d0
          end do          


          !initialize the test data for the inflow b.c.
          !------------------------------------------------
          flow_x_user = always_inflow

          !mass density test data
          test_data(1,1) = -6.02155766d0
          test_data(2,1) = -9.45357331d0
          test_data(3,1) = 8.378420918d0
          test_data(4,1) = 17.43249201d0       
          
          !momentum_x test data
          test_data(1,2) = -86.92413409d0
          test_data(2,2) =  48.43983127d0
          test_data(3,2) = -35.22995451d0
          test_data(4,2) =  97.43342286d0 

          !momentum_y test data
          test_data(1,3) = -54.02910608d0
          test_data(2,3) = -29.2929531d0 
          test_data(3,3) = -11.56590113d0
          test_data(4,3) = -31.08138214d0 

          !total energy test data
          test_data(1,4) = -303.6119102d0
          test_data(2,4) = -473.0514349d0
          test_data(3,4) = -460.1985142d0
          test_data(4,4) = -296.3032376d0


          !compare the test data with the timedev 
          !computed
          test_validated = test_compute_timedev_x_edge(
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         flux_y,
     $         transverse_lodi, viscous_lodi,
     $         inflow_bc, outflow_bc,
     $         flow_x_user,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_edge_inflow: '',L1)', test_validated
          end if


          !initialize the test data for the outflow b.c.
          !------------------------------------------------
          flow_x_user = always_outflow

          !mass density test data
          test_data(1,1) = 4.671816424d0
          test_data(2,1) = 3.379210546d0
          test_data(3,1) = 21.31254388d0
          test_data(4,1) = 28.13848309d0
          
          !momentum_x test data
          test_data(1,2) = -94.64365664d0
          test_data(2,2) =  24.75376184d0
          test_data(3,2) =  -12.0204285d0
          test_data(4,2) =  104.9436566d0

          !momentum_y test data
          test_data(1,3) =    441.68529d0
          test_data(2,3) = -305.3705731d0
          test_data(3,3) = -287.4372398d0
          test_data(4,3) =  465.1519567d0

          !total energy test data
          test_data(1,4) =  11023.75812d0
          test_data(2,4) = -1095.596288d0
          test_data(3,4) = -1077.662955d0
          test_data(4,4) =  11047.22478d0


          !compare the test data with the timedev 
          !computed
          test_validated = test_compute_timedev_x_edge(
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         flux_y,
     $         transverse_lodi, viscous_lodi,
     $         inflow_bc, outflow_bc,
     $         flow_x_user,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_edge_outflow: '',L1)', test_validated
          end if

        end function test_compute_timedev_edge


        function test_compute_timedev_x_edge(
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     flux_y,
     $     transverse_lodi, viscous_lodi,
     $     inflow_bc, outflow_bc,
     $     flow_x_user,
     $     test_data,
     $     detailled)
     $     result(test_validated)


          implicit none

          type(pmodel_eq)                   , intent(in) :: p_model
          real(rkind)                       , intent(in) :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in) :: nodes
          real(rkind), dimension(nx)        , intent(in) :: x_map
          real(rkind), dimension(ny)        , intent(in) :: y_map
          real(rkind), dimension(nx,ny+1,ne), intent(in) :: flux_y
          real(rkind), dimension(ne)        , intent(in) :: transverse_lodi
          real(rkind), dimension(ne)        , intent(in) :: viscous_lodi
          type(lodi_edge_inflow)            , intent(in) :: inflow_bc
          type(lodi_edge_outflow)           , intent(in) :: outflow_bc
          integer                           , intent(in) :: flow_x_user
          real(rkind), dimension(4,ne)      , intent(in) :: test_data
          logical                           , intent(in) :: detailled
          logical                                        :: test_validated


          integer(ikind)                   :: i,j
          logical                          :: side_x
          integer                          :: k
          real(rkind), dimension(nx,ny,ne) :: timedev
          real(rkind), dimension(4,ne)     :: timedev_output



          !compute the timedev for the test
          j=3

          side_x = left

          i=1
          call compute_timedev_x_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,i,j,
     $       flux_y,
     $       transverse_lodi, viscous_lodi,
     $       side_x,
     $       gradient_x_x_oneside_L0,
     $       inflow_bc,
     $       outflow_bc,
     $       flow_x_user,
     $       timedev)

          i=2
          call compute_timedev_x_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,i,j,
     $       flux_y,
     $       transverse_lodi, viscous_lodi,
     $       side_x,
     $       gradient_x_x_oneside_L1,
     $       inflow_bc,
     $       outflow_bc,
     $       flow_x_user,
     $       timedev)


          side_x = right

          i=4
          call compute_timedev_x_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,i,j,
     $       flux_y,
     $       transverse_lodi, viscous_lodi,
     $       side_x,
     $       gradient_x_x_oneside_R1,
     $       inflow_bc,
     $       outflow_bc,
     $       flow_x_user,
     $       timedev)

          i=5
          call compute_timedev_x_edge(
     $       p_model,
     $       t,nodes,x_map,y_map,i,j,
     $       flux_y,
     $       transverse_lodi, viscous_lodi,
     $       side_x,
     $       gradient_x_x_oneside_R0,
     $       inflow_bc,
     $       outflow_bc,
     $       flow_x_user,
     $       timedev)

          !compare the timedev of the test with the
          !test data
          do k=1, ne
             timedev_output(1,k) = timedev(1,j,k)
             timedev_output(2,k) = timedev(2,j,k)
             timedev_output(3,k) = timedev(4,j,k)
             timedev_output(4,k) = timedev(5,j,k)
          end do

          test_validated = is_matrix_validated(
     $         timedev_output,
     $         test_data,
     $         detailled)

        end function test_compute_timedev_x_edge


        !compute_timedev_corner
        function test_compute_timedev_corner(
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)                 , intent(in) :: p_model
          real(rkind)                     , intent(in) :: t
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          integer(ikind) :: i,j
          integer :: flow_x_user
          integer :: flow_y_user
          real(rkind), dimension(ne) :: test_data

          
          !test NE corner
          i = 2
          j = 4

          !inflow/inflow
          flow_x_user = always_inflow
          flow_y_user = always_inflow

          test_data(1) =  112.3433452d0
          test_data(2) = -26.69167932d0
          test_data(3) = -14.63876117d0
          test_data(4) =  114.3962942d0

          test_validated = test_compute_timedev_NW_corner(
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         i,j,
     $         flow_x_user, flow_y_user,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_corner_inflow_inflow: '',L1)', test_validated
          end if

          !inflow/outflow
          flow_x_user = always_outflow
          flow_y_user = always_inflow

          test_data(1) =  203.855464d0 
          test_data(2) =  12.4456627d0
          test_data(3) =-22.68219185d0
          test_data(4) = 138.3270377d0

          test_validated = test_compute_timedev_NW_corner(
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         i,j,
     $         flow_x_user, flow_y_user,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_corner_inflow_outflow: '',L1)', test_validated
          end if

          !outflow/outflow
          flow_x_user = always_outflow
          flow_y_user = always_outflow

          test_data(1) =  275.6367306d0
          test_data(2) = -74.26416037d0
          test_data(3) = -77.88518611d0
          test_data(4) =  108.9000925d0

          test_validated = test_compute_timedev_NW_corner(
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         i,j,
     $         flow_x_user, flow_y_user,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_corner_inflow_outflow: '',L1)', test_validated
          end if

        end function test_compute_timedev_corner


        function test_compute_timedev_NW_corner(
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     i,j,
     $     flow_x_user, flow_y_user,
     $     test_data,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)                  , intent(in) :: p_model
          real(rkind)                      , intent(in) :: t
          real(rkind), dimension(nx,ny,ne) , intent(in) :: nodes
          real(rkind), dimension(nx)       , intent(in) :: x_map
          real(rkind), dimension(ny)       , intent(in) :: y_map
          integer(ikind)                   , intent(in) :: i,j
          integer                          , intent(in) :: flow_x_user
          integer                          , intent(in) :: flow_y_user
          real(rkind), dimension(ne)       , intent(in) :: test_data
          logical                          , intent(in) :: detailled
          logical                                       :: test_validated

          
          type(lodi_corner_inflow_inflow)   :: corner_inflow_inflow
          type(lodi_corner_inflow_outflow)  :: corner_inflow_outflow
          type(lodi_corner_outflow_outflow) :: corner_outflow_outflow

          real(rkind), dimension(nx,ny,ne) :: timedev


          !compute the timedev at NE corner
          call compute_timedev_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         left, right,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         corner_inflow_inflow,
     $         corner_inflow_outflow,
     $         corner_outflow_outflow,
     $         flow_x_user, flow_y_user,
     $         timedev)

          !compare with the test data
          test_validated = is_vector_validated(
     $         timedev(i,j,:),
     $         test_data,
     $         detailled)

        end function test_compute_timedev_NW_corner

      end program test_yoo_lodi_timedev_xy
