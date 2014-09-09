      program test_yoo_ns2d_edge

        use ns2d_parameters, only :
     $       gamma,
     $       mach_infty

        use lodi_edge_abstract_class, only :
     $       lodi_edge_abstract

        use lodi_edge_inflow_class, only :
     $       lodi_edge_inflow

        use lodi_edge_outflow_class, only :
     $       lodi_edge_outflow

        use parameters_constant, only :
     $       x_direction,
     $       vector_x,
     $       left,
     $       right

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
        type(lodi_edge_outflow)            :: outflow_bc
        type(lodi_edge_inflow)             :: inflow_bc

        character(*), parameter :: FMT='(5F14.5)'
        
        real(rkind), dimension(nx,ny,ne) :: test_data
        logical                          :: detailled
        logical                          :: test_validated

        real(rkind), dimension(nx,ny,ne) :: transverse_lodi
        real(rkind), dimension(nx,ny,ne) :: viscous_lodi


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


        !compute the lodi vector from the lodi outflow x
        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)        
        call print_nodes(nodes,x_map,y_map)
        call get_test_data_for_lodi_outflow_x(test_data)

        print '(''test lodi for outflow x'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        call outflow_bc%ini()

        test_validated = test_lodi_x(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       transverse_lodi, viscous_lodi,
     $       detailled)

        print '()'


        !test the computation of the time derivatives for the lodi
        !outflow x

        print '(''test time_dev for outflow x'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)        
        call get_test_data_for_lodi_outflow_timedevx(test_data)

        detailled = .false.

        call outflow_bc%ini()

        test_validated = test_lodi_timedev_x(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       transverse_lodi, viscous_lodi,
     $       detailled)

        print '()'


        !compute the lodi vector from the lodi inflow x
        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)
        call get_test_data_for_lodi_inflow_x(test_data)

        print '(''test lodi for inflow x'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        call inflow_bc%ini()

        test_validated = test_lodi_x(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       transverse_lodi, viscous_lodi,
     $       detailled)

        print '()'



        !test the computation of the time derivatives for the lodi
        !inflow x

        print '(''test time_dev for inflow x'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call get_test_data_for_lodi_inflow_timedevx(test_data)

        detailled = .false.

        call inflow_bc%ini()

        test_validated = test_lodi_timedev_x(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       transverse_lodi, viscous_lodi,
     $       detailled)

        print '()'


        contains

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


       subroutine initialize_lodi_intermediate(transverse_lodi, viscous_lodi)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: transverse_lodi
          real(rkind), dimension(nx,ny,ne), intent(out) :: viscous_lodi


          integer(ikind) :: i,j
          integer        :: k

          
          do k=1, ne
             do j=1, ny
                do i=1, nx
                   transverse_lodi(i,j,k) = 0.1d0
                   viscous_lodi(i,j,k) = 0.2d0
                end do
             end do
          end do

       end subroutine initialize_lodi_intermediate


       subroutine print_nodes(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map

          integer(ikind) :: j


          print '(''x_map'')'
          print FMT, x_map
          print '()'

          print '(''y_map'')'
          print FMT, y_map
          print '()'

          print '()'
          print '(''mass_density'')'
          do j=1,5
             print FMT, nodes(1:5,6-j,1)
          end do
          print '()'

          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print FMT, nodes(1:5,6-j,2)
          end do
          print '()'

          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print FMT, nodes(1:5,6-j,3)
          end do
          print '()'

          print '()'
          print '(''total energy'')'
          do j=1,5
             print FMT, nodes(1:5,6-j,4)
          end do
          print '()'
          print '()'

        end subroutine print_nodes


        subroutine print_timedev(timedev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: timedev

          integer(ikind) :: j


          print '(''time derivatives of governing variables'')'
          print '(''---------------------------------------'')'
          
          print '()'
          print '(''mass_density'')'
          do j=1,5
             print FMT, timedev(1:5,6-j,1)
          end do
          print '()'
          
          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print FMT, timedev(1:5,6-j,2)
          end do
          print '()'
          
          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print FMT, timedev(1:5,6-j,3)
          end do
          print '()'
          
          print '()'
          print '(''total energy'')'
          do j=1,5
             print FMT, timedev(1:5,6-j,4)
          end do
          print '()'
          print '()'

        end subroutine print_timedev


        function test_lodi_x(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     transverse_lodi,
     $     viscous_lodi,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_edge_abstract)       , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)  :: transverse_lodi
          real(rkind), dimension(nx,ny,ne), intent(in)  :: viscous_lodi
          logical                         , intent(in)  :: detailled
          logical                                       :: test_validated


          real(rkind), dimension(nx,ny,ne) :: lodi
          logical                          :: loc
          logical                          :: test_lodi_validated
          logical, dimension(ne)           :: detailled_loc

          integer(ikind) :: i,j
          integer        :: k

          
          do j=1,5
             i=1
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            left,
     $            gradient_x_x_oneside_L0)
             
             i=2
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            left,
     $            gradient_x_x_oneside_L1)
             
             i=4
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            right,
     $            gradient_x_x_oneside_R1)

             i=5
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            right,
     $            gradient_x_x_oneside_R0)

          end do


          test_validated = .true.
          detailled_loc = [.false.,.false.,.false.,.false.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,5
                do i=1,2
                   loc = is_test_validated(lodi(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''lodi('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             
                do i=4,5
                   loc = is_test_validated(lodi(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''lodi('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             end do
             if(.not.detailled_loc(k)) then
                print '(''lodi('',I1,''):'',L3)', k, test_lodi_validated
             end if
          end do

          if(.not.detailled) print '(''test_validated: '',L3)', test_validated

        end function test_lodi_x


        function test_lodi_timedev_x(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     transverse_lodi, viscous_lodi,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_edge_abstract)       , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
          real(rkind), dimension(ne)      , intent(in)  :: transverse_lodi
          real(rkind), dimension(ne)      , intent(in)  :: viscous_lodi
          logical                         , intent(in)  :: detailled
          logical                                       :: test_validated


          real(rkind), dimension(nx,ny,ne) :: timedev
          logical                          :: loc
          logical                          :: test_lodi_validated
          logical, dimension(ne)           :: detailled_loc

          integer(ikind) :: i,j
          integer        :: k

          
          do j=1,5
             i=1
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            left,
     $            gradient_x_x_oneside_L0)
             
             i=2
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            left,
     $            gradient_x_x_oneside_L1)
             
             i=4
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            right,
     $            gradient_x_x_oneside_R1)

             i=5
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            right,
     $            gradient_x_x_oneside_R0)

          end do


          test_validated = .true.
          detailled_loc = [.false.,.false.,.false.,.false.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,5
                do i=1,2
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''timedev('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             
                do i=4,5
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''timedev('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             end do
             if(.not.detailled_loc(k)) then
                print '(''timedev('',I1,''):'',L3)', k, test_lodi_validated
             end if
          end do

          if(.not.detailled) print '(''test_validated: '',L3)', test_validated

        end function test_lodi_timedev_x


        subroutine get_test_data_for_lodi_outflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !L1
          test_data(1,5,1) =  446.77785d0
          test_data(1,4,1) = -1.2737843d0
          test_data(1,3,1) = 63.7416509d0
          test_data(1,2,1) = 81.1066497d0
          test_data(1,1,1) = -15.193722d0
                                       
          test_data(2,5,1) = -82.306674d0
          test_data(2,4,1) = -1.1667903d0
          test_data(2,3,1) = 69.9870196d0
          test_data(2,2,1) = -6.3174795d0
          test_data(2,1,1) = 10.5902500d0
                                       
          test_data(4,5,1) = -82.306674d0
          test_data(4,4,1) = -1.1667903d0
          test_data(4,3,1) = 69.9870196d0
          test_data(4,2,1) = -6.3174795d0
          test_data(4,1,1) = 10.5902500d0
                                       
          test_data(5,5,1) = 446.777854d0
          test_data(5,4,1) = -1.2737843d0
          test_data(5,3,1) = 63.7416509d0
          test_data(5,2,1) = 81.1066497d0
          test_data(5,1,1) = -15.193722d0
                                       
          !L2                          
          test_data(1,5,2) = 612.011517d0
          test_data(1,4,2) = 60.8978772d0
          test_data(1,3,2) = -1973.1928d0
          test_data(1,2,2) = -3957.7409d0
          test_data(1,1,2) = 3753.61415d0
                                       
          test_data(2,5,2) = 1648.38519d0
          test_data(2,4,2) = -187.99854d0
          test_data(2,3,2) = -501.81569d0
          test_data(2,2,2) = 29.5245331d0
          test_data(2,1,2) = 2106.98704d0
                                       
          test_data(4,5,2) = 1648.38519d0
          test_data(4,4,2) = -187.99854d0
          test_data(4,3,2) = -501.81569d0
          test_data(4,2,2) = 29.5245331d0
          test_data(4,1,2) = 2106.98704d0
                                       
          test_data(5,5,2) = 612.011517d0
          test_data(5,4,2) = 60.8978772d0
          test_data(5,3,2) = -1973.1928d0
          test_data(5,2,2) = -3957.7409d0
          test_data(5,1,2) = 3753.61415d0
                                       
          !L3                      
          test_data(1,5,3) = -3872.8478d0
          test_data(1,4,3) = 69.6020753d0
          test_data(1,3,3) = -2763.6543d0
          test_data(1,2,3) = 1973.37048d0
          test_data(1,1,3) = -6295.0801d0
                                       
          test_data(2,5,3) = 1367.715676d0
          test_data(2,4,3) = -299.003693d0
          test_data(2,3,3) = 497.1727008d0
          test_data(2,2,3) = -11.7202103d0
          test_data(2,1,3) = -15463.5656d0
                       
          test_data(4,5,3) = 98.40374684d0
          test_data(4,4,3) = -7.49140142d0
          test_data(4,3,3) = 63.94362383d0
          test_data(4,2,3) = -8.81893645d0
          test_data(4,1,3) = 417.8827844d0
                             
          test_data(5,5,3) = -8.8057865d0
          test_data(5,4,3) = 14.9460083d0
          test_data(5,3,3) = 19.5734597d0
          test_data(5,2,3) = 9.17385277d0
          test_data(5,1,3) = -5.4297627d0
                                       
          !L4              
          test_data(1,5,4) = -8.8057865d0
          test_data(1,4,4) = 14.9460083d0
          test_data(1,3,4) = 19.5734597d0
          test_data(1,2,4) = 9.17385277d0
          test_data(1,1,4) = -5.4297627d0
                                       
          test_data(2,5,4) = 98.40374684d0
          test_data(2,4,4) = -7.49140142d0
          test_data(2,3,4) = 63.94362383d0
          test_data(2,2,4) = -8.81893645d0
          test_data(2,1,4) = 417.8827844d0
                                       
          test_data(4,5,4) = 1367.715676d0
          test_data(4,4,4) = -299.003693d0
          test_data(4,3,4) = 497.1727008d0
          test_data(4,2,4) = -11.7202103d0
          test_data(4,1,4) = -15463.5656d0
                       
          test_data(5,5,4) = -3872.8478d0
          test_data(5,4,4) = 69.6020753d0
          test_data(5,3,4) = -2763.6543d0
          test_data(5,2,4) = 1973.37048d0
          test_data(5,1,4) = -6295.0801d0

        end subroutine get_test_data_for_lodi_outflow_x


        subroutine get_test_data_for_lodi_outflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data


          !mass
          test_data(1,5,1) = 44.10934585d0
          test_data(1,4,1) = -15.2165587d0
          test_data(1,3,1) = 11.25514976d0
          test_data(1,2,1) = 79.09267117d0
          test_data(1,1,1) = -110.744439d0

          test_data(2,5,1) = -35.2692175d0
          test_data(2,4,1) = 341.6286828d0
          test_data(2,3,1) = 5.979210546d0
          test_data(2,2,1) = -33.0248749d0
          test_data(2,1,1) = 5.796465574d0
                                      
          test_data(4,5,1) = -35.2692175d0
          test_data(4,4,1) = 341.6286828d0
          test_data(4,3,1) = 5.979210546d0
          test_data(4,2,1) = -33.0248749d0
          test_data(4,1,1) = 5.796465574d0
                                      
          test_data(5,5,1) = 44.10934585d0
          test_data(5,4,1) = -15.2165587d0
          test_data(5,3,1) = 11.25514976d0
          test_data(5,2,1) = 79.09267117d0
          test_data(5,1,1) = -110.744439d0

                                 
          !momentum-x
          test_data(1,5,2) = -623.274248d0
          test_data(1,4,2) = -4.23812689d0 
          test_data(1,3,2) = -88.0603233d0
          test_data(1,2,2) = -254.873883d0
          test_data(1,1,2) = -1056.25085d0

          test_data(2,5,2) = 0.818692465d0
          test_data(2,4,2) = -13.9166020d0
          test_data(2,3,2) = 27.35376184d0
          test_data(2,2,2) = -32.5438124d0
          test_data(2,1,2) = -213.608694d0
                             
          test_data(4,5,2) = -0.81869246d0
          test_data(4,4,2) = 13.91660206d0  
          test_data(4,3,2) = -27.3537618d0
          test_data(4,2,2) = 32.54381247d0  
          test_data(4,1,2) = 213.6086941d0
                             
          test_data(5,5,2) = 623.274248d0
          test_data(5,4,2) = 4.23812689d0
          test_data(5,3,2) = 88.0603233d0
          test_data(5,2,2) = 254.873883d0
          test_data(5,1,2) = 1056.25085d0

                                 
          !momentum-y
          test_data(1,5,3) = 1256.48673d0
          test_data(1,4,3) = 7.01407643d0 
          test_data(1,3,3) = 448.268623d0
          test_data(1,2,3) = 372.943861d0
          test_data(1,1,3) = -61.835447d0
                             
          test_data(2,5,3) = 245.0656909d0
          test_data(2,4,3) = -85.0659158d0
          test_data(2,3,3) = -302.770573d0
          test_data(2,2,3) = 62.04892821d0
          test_data(2,1,3) = 7.667075478d0
                             
          test_data(4,5,3) = 245.0656909d0
          test_data(4,4,3) = -85.0659158d0
          test_data(4,3,3) = -302.770573d0
          test_data(4,2,3) = 62.04892821d0
          test_data(4,1,3) = 7.667075478d0
                             
          test_data(5,5,3) = 1256.48673d0 
          test_data(5,4,3) = 7.01407643d0
          test_data(5,3,3) = 448.268623d0 
          test_data(5,2,3) = 372.943861d0 
          test_data(5,1,3) = -61.835447d0

          !total energy
          test_data(1,5,4) = 23715.67526d0
          test_data(1,4,4) = -58.1465572d0
          test_data(1,3,4) = 11030.34145d0
          test_data(1,2,4) = -419.433190d0
          test_data(1,1,4) = 7872.233274d0
                             
          test_data(2,5,4) = -621.786095d0
          test_data(2,4,4) = 209.4959552d0
          test_data(2,3,4) = -1092.99628d0
          test_data(2,2,4) = -14.1760078d0
          test_data(2,1,4) = 9393.828732d0
                             
          test_data(4,5,4) = -621.786095d0
          test_data(4,4,4) = 209.4959552d0
          test_data(4,3,4) = -1092.99628d0
          test_data(4,2,4) = -14.1760078d0
          test_data(4,1,4) = 9393.828732d0
                             
          test_data(5,5,4) = 23715.67526d0
          test_data(5,4,4) = -58.1465572d0
          test_data(5,3,4) = 11030.34145d0
          test_data(5,2,4) = -419.433190d0
          test_data(5,1,4) = 7872.233274d0

        end subroutine get_test_data_for_lodi_outflow_timedevx


        subroutine get_test_data_for_lodi_inflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data


          !L1
          test_data(1,5,1) = 540.4347568d0
          test_data(1,4,1) = 0.562369759d0
          test_data(1,3,1) = 306.7471758d0
          test_data(1,2,1) = 24.23568908d0
          test_data(1,1,1) = 3.421977833d0
                             
          test_data(2,5,1) = 1.70950990d0
          test_data(2,4,1) = -0.1563713d0
          test_data(2,3,1) = 2.83158308d0
          test_data(2,2,1) = -0.7406068d0
          test_data(2,1,1) = 3.06662209d0
                             
          test_data(4,5,1) = 1.70950990d0
          test_data(4,4,1) = -0.1563713d0
          test_data(4,3,1) = 2.83158308d0
          test_data(4,2,1) = -0.7406068d0
          test_data(4,1,1) = 3.06662209d0
                             
          test_data(5,5,1) = 540.4347568d0
          test_data(5,4,1) = 0.562369759d0
          test_data(5,3,1) = 306.7471758d0
          test_data(5,2,1) = 24.23568908d0
          test_data(5,1,1) = 3.421977833d0

          !L2
          test_data(1,5,2) = 54.23088137d0
          test_data(1,4,2) = -12.8598311d0
          test_data(1,3,2) = 1217.924473d0
          test_data(1,2,2) = 30.54118812d0
          test_data(1,1,2) = -43.3197689d0
                                      
          test_data(2,5,2) =  25.57906023d0
          test_data(2,4,2) = -21.11020646d0
          test_data(2,3,2) =   9.36048316d0
          test_data(2,2,2) = -57.19780745d0
          test_data(2,1,2) =  478.8211866d0 
                                         
          test_data(4,5,2) =  25.57906023d0
          test_data(4,4,2) = -21.11020646d0
          test_data(4,3,2) =   9.36048316d0
          test_data(4,2,2) = -57.19780745d0
          test_data(4,1,2) =  478.8211866d0
                                      
          test_data(5,5,2) = 54.23088137d0
          test_data(5,4,2) = -12.8598311d0
          test_data(5,3,2) = 1217.924473d0
          test_data(5,2,2) = 30.54118812d0
          test_data(5,1,2) = -43.3197689d0

          !L3
          test_data(1,5,3) = -3872.8478d0
          test_data(1,4,3) = 69.6020753d0 
          test_data(1,3,3) = -2763.6543d0
          test_data(1,2,3) = 1973.37048d0
          test_data(1,1,3) = -6295.0801d0
                                         
          test_data(2,5,3) = 1367.715676d0 
          test_data(2,4,3) = -299.003693d0 
          test_data(2,3,3) = 497.1727008d0 
          test_data(2,2,3) = -11.7202103d0
          test_data(2,1,3) = -15463.5656d0   

          test_data(4,5,3) = 12.14325378d0
          test_data(4,4,3) = 5.450606392d0
          test_data(4,3,3) = -1.16793344d0
          test_data(4,2,3) = 7.404333897d0
          test_data(4,1,3) = 33.89325625d0
                                          
          test_data(5,5,3) = -19.5795799d0
          test_data(5,4,3) = 7.642473002d0
          test_data(5,3,3) = 1.365742398d0
          test_data(5,2,3) = -15.7670170d0
          test_data(5,1,3) = -5.92464050d0

          !L4
          test_data(1,5,4) =-27.0795799d0
          test_data(1,4,4) =0.142473002d0
          test_data(1,3,4) =-6.13425760d0
          test_data(1,2,4) =-23.2670170d0
          test_data(1,1,4) =-13.4246405d0
                                         
          test_data(2,5,4) =4.643253785d0
          test_data(2,4,4) =-2.04939360d0
          test_data(2,3,4) =-8.66793344d0
          test_data(2,2,4) =-0.09566610d0
          test_data(2,1,4) =26.39325625d0
                                           
          test_data(4,5,4) =1367.715676d0
          test_data(4,4,4) =-299.003693d0
          test_data(4,3,4) =497.1727008d0
          test_data(4,2,4) =-11.7202103d0
          test_data(4,1,4) =-15463.5656d0 
                                         
          test_data(5,5,4) =-3872.8478d0 
          test_data(5,4,4) =69.6020753d0 
          test_data(5,3,4) =-2763.6543d0 
          test_data(5,2,4) =1973.37048d0 
          test_data(5,1,4) =-6295.0801d0 

        end subroutine get_test_data_for_lodi_inflow_x


        subroutine get_test_data_for_lodi_inflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !mass
          test_data(1,5,1) = 62.9278839d0
          test_data(1,4,1) = -3.2465581d0
          test_data(1,3,1) = 0.56177567d0
          test_data(1,2,1) = -26.811349d0
          test_data(1,1,1) = 586.903034d0
                                 
          test_data(2,5,1) = -10.5411492d0
          test_data(2,4,1) = 171.8291839d0
          test_data(2,3,1) = -6.85357331d0
          test_data(2,2,1) = 108.2349372d0
          test_data(2,1,1) = 7.748555568d0
                   
          test_data(4,5,1) = -10.5966868d0
          test_data(4,4,1) = 168.0749796d0
          test_data(4,3,1) = -6.95491241d0
          test_data(4,2,1) = 101.8031772d0
          test_data(4,1,1) = 7.744542029d0
                                     
          test_data(5,5,1) = 62.80340463d0
          test_data(5,4,1) = -3.79963590d0
          test_data(5,3,1) = 0.549158673d0
          test_data(5,2,1) = -26.9113328d0
          test_data(5,1,1) = 586.2147348d0
               
                    
          !momentum-x        
          test_data(1,5,2) = -737.343573d0
          test_data(1,4,2) = 10.19423164d0
          test_data(1,3,2) = -80.3408007d0
          test_data(1,2,2) = 303.7707912d0
          test_data(1,1,2) = -2888.13781d0
                             
          test_data(2,5,2) = 60.10132975d0
          test_data(2,4,2) = -82.2078604d0
          test_data(2,3,2) = 51.03983127d0
          test_data(2,2,2) = 92.82001573d0
          test_data(2,1,2) = -191.655149d0

          test_data(4,5,2) = -59.524637d0
          test_data(4,4,2) = 87.4096625d0
          test_data(4,3,2) = -50.563287d0
          test_data(4,2,2) = -81.940820d0
          test_data(4,1,2) = 191.809801d0
                            
          test_data(5,5,2) = 737.2612518d0     
          test_data(5,4,2) = -8.21856887d0     
          test_data(5,3,2) = 80.55008953d0     
          test_data(5,2,2) = -303.683382d0
          test_data(5,1,2) = 2887.935375d0


          !momentum-y
          test_data(1,5,3) = 1811.93718d0
          test_data(1,4,3) = -5.5936922d0
          test_data(1,3,3) = -47.445772d0
          test_data(1,2,3) = -188.49864d0
          test_data(1,1,3) = 505.031667d0
                             
          test_data(2,5,3) = -21.4681830d0
          test_data(2,4,3) = -46.1434760d0
          test_data(2,3,3) = -26.6929531d0
          test_data(2,2,3) = -18.9275509d0
          test_data(2,1,3) = 21.56107056d0
                             
          test_data(4,5,3) = -21.5488667d0
          test_data(4,4,3) = -45.1086632d0
          test_data(4,3,3) = -26.8992344d0
          test_data(4,2,3) = -17.4243140d0
          test_data(4,1,3) = 21.54828232d0
                             
          test_data(5,5,3) = 1808.13912d0
          test_data(5,4,3) = -5.7322870d0
          test_data(5,3,3) = -47.964715d0
          test_data(5,2,3) = -189.09312d0
          test_data(5,1,3) = 504.430154d0


          !total energy
          test_data(1,5,4) = 32263.07675d0
          test_data(1,4,4) = -42.2159859d0
          test_data(1,3,4) = -297.028576d0
          test_data(1,2,4) = -3334.73191d0
          test_data(1,1,4) = 10512.30477d0
                                    
          test_data(2,5,4) = -894.371848d0
          test_data(2,4,4) = 187.4252157d0
          test_data(2,3,4) = -470.451434d0
          test_data(2,2,4) = 49.86206488d0
          test_data(2,1,4) = 9834.749701d0
                             
          test_data(4,5,4) = -901.174597d0
          test_data(4,4,4) = 179.9288051d0
          test_data(4,3,4) = -475.531847d0
          test_data(4,2,4) = 36.73542412d0
          test_data(4,1,4) = 9827.999751d0
                             
          test_data(5,5,4) = 32201.3574d0
          test_data(5,4,4) = -49.512007d0
          test_data(5,3,4) = -313.18657d0
          test_data(5,2,4) = -3340.2874d0
          test_data(5,1,4) = 10508.2621d0 

        end subroutine get_test_data_for_lodi_inflow_timedevx

      end program test_yoo_ns2d_edge
