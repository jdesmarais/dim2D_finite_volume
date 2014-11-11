      program test_yoo_sym_x

        use bc_operators_class, only :
     $     bc_operators

        use ns2d_parameters, only :
     $       gamma,
     $       mach_infty

        use parameters_constant, only :
     $       left,right,
     $       vector_x,
     $       x_direction,
     $       always_inflow,
     $       always_outflow

        use parameters_input, only :
     $       sigma_P,
     $       flow_direction,
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none


        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind), dimension(nx)         :: x_map
        real(rkind), dimension(ny)         :: y_map
        real(rkind)                        :: t
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        type(pmodel_eq)                    :: p_model
        real(rkind), dimension(nx+1,ny,ne) :: flux_x
        real(rkind), dimension(nx,ny+1,ne) :: flux_y
        real(rkind), dimension(nx,ny,ne)   :: timedev
        type(bc_operators)                 :: bc_used

        character(*), parameter :: FMT='(5F14.5)'

        real(rkind), dimension(nx,ny,ne) :: test_data

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
     $       (.not.(flow_direction.eq.x_direction)).or.
     $       (.not.is_test_validated(p_model%get_mach_ux_infty(left),-mach_infty,detailled)).or.
     $       (.not.is_test_validated(p_model%get_mach_ux_infty(right),mach_infty,detailled)).or.
     $       (.not.is_test_validated(p_model%get_mach_uy_infty(left),0.0d0,detailled)).or.
     $       (.not.is_test_validated(p_model%get_mach_uy_infty(right),0.0d0,detailled))) then

           print '(''the test requires:'')'
           print '(''gamma=5/3: '',L1)', is_test_validated(gamma,5.0d0/3.0d0,detailled)
           print '(''mach_infty=0.2: '',L1)', is_test_validated(mach_infty,0.2d0,detailled)
           print '(''sigma_P=0.25: '',L1)', is_test_validated(sigma_P,0.25d0,detailled)
           print '(''flow_direction=x-direction: '',L1)', (flow_direction.eq.x_direction)
           print '(''ic_choice=sym_x: '')'
           stop ''

        end if


        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call print_nodes(nodes,x_map,y_map)


        !test apply_bc_on_timedev
        print '(''test apply_bc_on_timedev'')'
        print '(''---------------------------------------'')'

        !inflow/inflow edges
        detailled = .false.

        call bc_used%ini(p_model)
        call bc_used%set_obc_type([
     $       always_inflow,
     $       always_inflow,
     $       always_inflow,
     $       always_inflow])

        call bc_used%apply_bc_on_timedev(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x,flux_y,
     $       timedev)

        call get_edge_inflow_inflow_test_data(test_data)

        test_validated = test_edge_x(timedev,test_data,detailled)

        if(.not.detailled) then
           print '(''test_edge_inflow_inflow: '',L1)', test_validated
        end if
        print '()'

        
        !outflow/outflow edges
        detailled = .false.

        call bc_used%ini(p_model)
        call bc_used%set_obc_type([
     $       always_outflow,
     $       always_outflow,
     $       always_outflow,
     $       always_outflow])

        call bc_used%apply_bc_on_timedev(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x,flux_y,
     $       timedev)

        call get_edge_outflow_outflow_test_data(test_data)

        test_validated = test_edge_x(timedev,test_data,detailled)

        if(.not.detailled) then
           print '(''test_edge_inflow_inflow: '',L1)', test_validated
        end if
        print '()'

        call print_timedev(timedev)


        contains


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
             x_map(i) = (i-3)*dx
          end do

          !initialize the y_map
          do i=1,5
             y_map(i) = (i-3)*dy
          end do
          
       end subroutine initialize_nodes


       !print nodes
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


       !print timedev
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


        subroutine get_edge_inflow_inflow_test_data(test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          test_data(1,3,1) =  5.175782381d0
          test_data(2,3,1) = -7.239125586d0
          test_data(4,3,1) = -7.239125586d0
          test_data(5,3,1) =  5.175782381d0

          test_data(1,3,2) = -63.69487469d0
          test_data(2,3,2) =  52.30399907d0
          test_data(4,3,2) = -52.30399907d0
          test_data(5,3,2) =  63.69487469d0

          test_data(1,3,3) = -198.7475544d0
          test_data(2,3,3) =  -27.0085677d0
          test_data(4,3,3) =  -27.0085677d0
          test_data(5,3,3) = -198.7475544d0

          test_data(1,3,4) = -6053.037008d0
          test_data(2,3,4) = -494.9285596d0
          test_data(4,3,4) = -494.9285596d0
          test_data(5,3,4) = -6053.037008d0

        end subroutine get_edge_inflow_inflow_test_data


        subroutine get_edge_outflow_outflow_test_data(test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          test_data(1,3,1) =   18.3065111d0
          test_data(1,3,2) = -72.86211727d0
          test_data(1,3,3) =  568.1087269d0
          test_data(1,3,4) =  14368.95909d0

          test_data(2,3,1) =  6.233415707d0
          test_data(2,3,2) =  28.58576219d0
          test_data(2,3,3) = -305.2461723d0
          test_data(2,3,4) = -1115.994698d0

          test_data(4,3,1) =  6.233415707d0
          test_data(4,3,2) = -28.58576219d0
          test_data(4,3,3) = -305.2461723d0
          test_data(4,3,4) = -1115.994698d0

          test_data(5,3,1) =   18.3065111d0
          test_data(5,3,2) =  72.86211727d0
          test_data(5,3,3) =  568.1087269d0
          test_data(5,3,4) =  14368.95909d0

        end subroutine get_edge_outflow_outflow_test_data


        function test_edge_x(timedev,test_data,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: timedev
          real(rkind), dimension(nx,ny,ne), intent(in) :: test_data
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          logical        :: test_loc
          integer(ikind) :: i,j
          integer        :: k          
          
          test_validated = .true.

          j=3
          do k=1,ne
             do i=1,2
                test_loc = is_test_validated(
     $               timedev(i,j,k), test_data(i,j,k), detailled)
                test_validated = test_validated.and.test_loc
                if(detailled) then
                   print '(''timedev('',I3,'')'')', i,j,k,test_validated
                end if
             end do

             do i=4,5
                test_loc = is_test_validated(
     $               timedev(i,j,k), test_data(i,j,k), detailled)
                test_validated = test_validated.and.test_loc
                if(detailled) then
                   print '(''timedev('',I3,'')'')', i,j,k,test_validated
                end if
             end do
          end do

        end function test_edge_x

      end program test_yoo_sym_x
