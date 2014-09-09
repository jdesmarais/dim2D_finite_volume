      program test_yoo_corners

        use lodi_corner_inflow_inflow_class, only :
     $       get_lodi_A_inflow_inflow,
     $       lodi_corner_inflow_inflow

        use lodi_corner_inflow_outflow_class, only :
     $       lodi_corner_inflow_outflow
        
        use lodi_corner_outflow_outflow_class, only :
     $       lodi_corner_outflow_outflow

        use ns2d_parameters, only :
     $       gamma,
     $       mach_infty

        use parameters_constant, only :
     $       left,
     $       right,
     $       inflow_type,
     $       outflow_type,
     $       vector_x,
     $       x_direction

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
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1

        implicit none

        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind), dimension(nx)         :: x_map
        real(rkind), dimension(ny)         :: y_map
        real(rkind)                        :: t
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        type(pmodel_eq)                    :: p_model
        

        !character(*), parameter :: FMT='(5F14.5)'

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


        !test get_lodi_A
        print '(''test get_lodi_A'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        test_validated = test_get_lodi_A(detailled)

        print '()'


        !test corner inflow/inflow
        print '(''test corner_inflow_inflow'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

        test_validated = test_corner_inflow_inflow(
     $       p_model,
     $       nodes, x_map, y_map,
     $       detailled)

        if(.not.detailled) then
           print '(''---------------------------------------'')'
           print '(''test_validated:'',L2)', test_validated
        end if
        print '()'
        

        !test corner inflow/outflow
        print '(''test corner_inflow_outflow'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

        !detailled = .true.

        test_validated = test_corner_inflow_outflow(
     $       p_model,
     $       nodes, x_map, y_map,
     $       detailled)

        if(.not.detailled) then
           print '(''---------------------------------------'')'
           print '(''test_validated:'',L2)', test_validated
        end if
        print '()'


        !test corner outflow/outflow
        print '(''test corner_outflow_outflow'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

        !detailled = .true.

        test_validated = test_corner_outflow_outflow(
     $       p_model,
     $       nodes, x_map, y_map,
     $       detailled)

        if(.not.detailled) then
           print '(''---------------------------------------'')'
           print '(''test_validated:'',L2)', test_validated
        end if
        print '()'

        contains

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


        function test_get_lodi_A(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6,6) :: lodi_A

          real(rkind), dimension(6,6) :: test_data
          real(rkind)                 :: md
          real(rkind)                 :: c
          real(rkind)                 :: tau
          integer                     :: sign_x
          integer                     :: sign_y

          !inputs
          md     =  0.3d0
          c      =  0.5d0
          tau    =  0.2d0
          sign_x =  1
          sign_y = -1

          !fct tested
          lodi_A = get_lodi_A_inflow_inflow(
     $         md,c,sign_x,sign_y,tau)
          
          !data output expected
          test_data(1,1) =  1.71957672d0
          test_data(2,1) =  0.0d0
          test_data(3,1) = -3.527336861d0
          test_data(4,1) =  0.423280423d0
          test_data(5,1) =  0.0d0
          test_data(6,1) =  5.996472663d0

          test_data(1,2) =  0.0d0
          test_data(2,2) =  2.777777778d0
          test_data(3,2) =  0.0d0
          test_data(4,2) =  0.0d0
          test_data(5,2) =  -2.22222222d0
          test_data(6,2) =  0.0d0

          test_data(1,3) = -0.158730159d0
          test_data(2,3) =  0.0d0
          test_data(3,3) =  2.248677249d0
          test_data(4,3) = -0.26984127d0
          test_data(5,3) =  0.0d0
          test_data(6,3) = -1.322751323d0

          test_data(1,4) =  0.423280423d0
          test_data(2,4) =  0.0d0
          test_data(3,4) = -5.996472663d0
          test_data(4,4) =  1.71957672d0
          test_data(5,4) =  0.0d0
          test_data(6,4) =  3.527336861d0

          test_data(1,5) =  0.0d0
          test_data(2,5) = -2.222222222d0
          test_data(3,5) =  0.0d0
          test_data(4,5) =  0.0d0
          test_data(5,5) =  2.777777778d0
          test_data(6,5) =  0.0d0

          test_data(1,6) =  0.26984127d0
          test_data(2,6) =  0.0d0
          test_data(3,6) = -1.322751323d0
          test_data(4,6) =  0.158730159d0
          test_data(5,6) =  0.0d0
          test_data(6,6) =  2.248677249d0

          !test the output
          test_validated = is_matrix_validated(
     $         lodi_A,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_get_lodi_A: '',L2)', test_validated
          end if

        end function test_get_lodi_A


        function test_corner_inflow_inflow(
     $     p_model,
     $     nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          integer(ikind)             :: i
          integer(ikind)             :: j
          logical                    :: side_x
          logical                    :: side_y
          real(rkind), dimension(ne) :: lodi_x
          real(rkind), dimension(ne) :: lodi_y
          real(rkind), dimension(ne) :: test_data_lodi_x
          real(rkind), dimension(ne) :: test_data_lodi_y
          logical                    :: test_loc
          real(rkind), dimension(ne) :: timedev
          real(rkind), dimension(ne) :: test_data_timedev

          type(lodi_corner_inflow_inflow) :: corner_i_i

          
          test_validated = .true.


          !NW corner: (side_x=left,side_y=right)
          i = 2
          j = 4          
          side_x = left
          side_y = right


          !test the computation of the x and y LODI vectors
          call corner_i_i%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         left, right,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) =    1.454995493d0
          test_data_lodi_x(2) =  -19.17172083d0
          test_data_lodi_x(3) = -299.0036936d0
          test_data_lodi_x(4) =   15.99921919d0

          test_data_lodi_y(1) = -11.22002954d0
          test_data_lodi_y(2) =  -4.17102088d0
          test_data_lodi_y(3) =  80.28858505d0
          test_data_lodi_y(4) =  24.96631204d0          

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NW):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NW):'',L2)', test_loc


          !test the computation of the time derivatives from
          !the x and y LODI vectors
          timedev = corner_i_i%compute_x_and_y_timedev(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1)
          
          test_data_timedev(1) =  112.3433452d0
          test_data_timedev(2) = -26.69167932d0
          test_data_timedev(3) = -14.63876117d0
          test_data_timedev(4) =  114.3962942d0

          test_loc = is_vector_validated(
     $         timedev,
     $         test_data_timedev,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''timedev(NW):'',L2)', test_loc

          
          !NE corner: (side_x=right,side_y=right)
          i = 4
          j = 4          
          side_x = right
          side_y = right


          !test the computation of the x and y LODI vectors
          call corner_i_i%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) =  1.462867208d0
          test_data_lodi_x(2) = -19.17172083d0
          test_data_lodi_x(3) =  15.23239041d0
          test_data_lodi_x(4) = -299.0036936d0

          test_data_lodi_y(1) =  9.216255037d0
          test_data_lodi_y(2) =  -4.17102088d0
          test_data_lodi_y(3) =   80.5222253d0
          test_data_lodi_y(4) =  24.96631204d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NE):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NE):'',L2)', test_loc

          
          !SE corner: (side_x=right,side_y=left)
          i = 4
          j = 2          
          side_x = right
          side_y = left


          !test the computation of the x and y LODI vectors
          call corner_i_i%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -8.494118919d0
          test_data_lodi_x(2) = -71.64831958d0
          test_data_lodi_x(3) =  66.02236803d0
          test_data_lodi_x(4) = -11.72021033d0

          test_data_lodi_y(1) = -6.5319714d0
          test_data_lodi_y(2) = -56.38772216d0
          test_data_lodi_y(3) =  385.9079204d0
          test_data_lodi_y(4) = -7.247702158d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(SE):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(SE):'',L2)', test_loc

          
          !SE corner: (side_x=right,side_y=left)
          i = 4
          j = 2          
          side_x = right
          side_y = left


          !test the computation of the x and y LODI vectors
          call corner_i_i%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -8.494118919d0
          test_data_lodi_x(2) = -71.64831958d0
          test_data_lodi_x(3) =  66.02236803d0
          test_data_lodi_x(4) = -11.72021033d0

          test_data_lodi_y(1) = -6.5319714d0
          test_data_lodi_y(2) = -56.38772216d0
          test_data_lodi_y(3) =  385.9079204d0
          test_data_lodi_y(4) = -7.247702158d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(SE):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(SE):'',L2)', test_loc


          !SW corner: (side_x=right,side_y=left)
          i = 2
          j = 2          
          side_x = left
          side_y = left


          !test the computation of the x and y LODI vectors
          call corner_i_i%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         inflow_type, inflow_type,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -8.536994696d0
          test_data_lodi_x(2) = -71.64831958d0
          test_data_lodi_x(3) = -11.72021033d0
          test_data_lodi_x(4) =  49.01483604d0

          test_data_lodi_y(1) =   0.98256584d0
          test_data_lodi_y(2) = -56.38772216d0
          test_data_lodi_y(3) =  385.9079204d0
          test_data_lodi_y(4) = -9.471904869d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(SW):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(SW):'',L2)', test_loc

        end function test_corner_inflow_inflow


        function test_corner_inflow_outflow(
     $     p_model,
     $     nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          integer(ikind)             :: i
          integer(ikind)             :: j
          logical                    :: side_x
          logical                    :: side_y
          logical                    :: flow_x
          logical                    :: flow_y
          real(rkind), dimension(ne) :: lodi_x
          real(rkind), dimension(ne) :: lodi_y
          real(rkind), dimension(ne) :: test_data_lodi_x
          real(rkind), dimension(ne) :: test_data_lodi_y
          logical                    :: test_loc
c$$$          real(rkind), dimension(ne) :: timedev
c$$$          real(rkind), dimension(ne) :: test_data_timedev

          type(lodi_corner_inflow_outflow) :: corner_i_o

          
          test_validated = .true.


          !NW corner: (side_x=left,side_y=right)
          i = 2
          j = 4          
          side_x = left
          side_y = right


          !x:inflow/y:outflow
          print '(''x:inflow/y:outflow'')'
          print '(''=================='')'
          flow_x = inflow_type
          flow_y = outflow_type

          !test the computation of the x and y LODI vectors
          call corner_i_o%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         left, right,
     $         flow_x, flow_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -1.250055995d0
          test_data_lodi_x(2) =  -27.0137852d0
          test_data_lodi_x(3) = -299.0036936d0
          test_data_lodi_x(4) = -5.522545252d0

          test_data_lodi_y(1) = -0.814825499d0
          test_data_lodi_y(2) =  10.75786491d0
          test_data_lodi_y(3) =  0.0d0
          test_data_lodi_y(4) =  24.96631204d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NW):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NW):'',L2)', test_loc
          print '()'

          !x:inflow/y:outflow
          print '(''x:outflow/y:inflow'')'
          print '(''=================='')'
          flow_x = outflow_type
          flow_y = inflow_type

          !test the computation of the x and y LODI vectors
          call corner_i_o%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         left, right,
     $         flow_x, flow_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -1.166790347d0
          test_data_lodi_x(2) = -187.9985406d0
          test_data_lodi_x(3) = -299.0036936d0
          test_data_lodi_x(4) = 0.0d0         

          test_data_lodi_y(1) = -10.680990d0
          test_data_lodi_y(2) = 84.5128104d0
          test_data_lodi_y(3) = 73.7545150d0
          test_data_lodi_y(4) = 24.9663120d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NW):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NW):'',L2)', test_loc
          print '()'

        end function test_corner_inflow_outflow


        function test_corner_outflow_outflow(
     $     p_model,
     $     nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          integer(ikind)             :: i
          integer(ikind)             :: j
          logical                    :: side_x
          logical                    :: side_y
          logical                    :: flow_x
          logical                    :: flow_y
          real(rkind), dimension(ne) :: lodi_x
          real(rkind), dimension(ne) :: lodi_y
          real(rkind), dimension(ne) :: test_data_lodi_x
          real(rkind), dimension(ne) :: test_data_lodi_y
          logical                    :: test_loc

          type(lodi_corner_outflow_outflow) :: corner_o_o

          
          test_validated = .true.


          !NW corner: (side_x=left,side_y=right)
          i = 2
          j = 4          
          side_x = left
          side_y = right


          !x:outflow/y:outflow
          flow_x = outflow_type
          flow_y = outflow_type

          !test the computation of the x and y LODI vectors
          call corner_o_o%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         flow_x, flow_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -1.166790347d0
          test_data_lodi_x(2) = -187.9985406d0
          test_data_lodi_x(3) = -299.0036936d0
          test_data_lodi_x(4) =  74.90272475d0

          test_data_lodi_y(1) = -0.814825499d0
          test_data_lodi_y(2) =  10.75786491d0
          test_data_lodi_y(3) =  2.959925659d0
          test_data_lodi_y(4) =  24.96631204d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NW):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NW):'',L2)', test_loc


          !NW corner: (side_x=left,side_y=right)
          i = 4
          j = 4          
          side_x = right
          side_y = right


          !x:outflow/y:outflow
          flow_x = outflow_type
          flow_y = outflow_type

          !test the computation of the x and y LODI vectors
          call corner_o_o%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         flow_x, flow_y,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R1,
     $         lodi_x, lodi_y)

          test_data_lodi_x(1) = -1.166790347d0
          test_data_lodi_x(2) = -187.9985406d0
          test_data_lodi_x(3) =  74.90272475d0
          test_data_lodi_x(4) = -299.0036936d0

          test_data_lodi_y(1) =  0.814825499d0
          test_data_lodi_y(2) =  10.75786491d0
          test_data_lodi_y(3) =  2.959925659d0
          test_data_lodi_y(4) =  24.96631204d0

          test_loc = is_vector_validated(
     $         lodi_x,
     $         test_data_lodi_x,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_x(NE):'',L2)', test_loc

          test_loc = is_vector_validated(
     $         lodi_y,
     $         test_data_lodi_y,
     $         detailled)
          test_validated = test_validated.and.test_loc
          print '(''lodi_y(NE):'',L2)', test_loc

        end function test_corner_outflow_outflow

      end program test_yoo_corners
