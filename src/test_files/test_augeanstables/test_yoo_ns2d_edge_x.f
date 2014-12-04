      program test_yoo_ns2d_edge_x

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

        use test_yoo_ns2d_edge_module, only :
     $       initialize_nodes,
     $       initialize_lodi_intermediate,
     $       get_test_data_for_lodi_outflow_x,
     $       get_test_data_for_lodi_inflow_x,
     $       get_test_data_for_lodi_outflow_timedevx,
     $       get_test_data_for_lodi_inflow_timedevx

        implicit none

        character(*), parameter :: FMT='(5F14.5)'

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        test_validated = .true.

        !verify inputs
        detailled = .false.
        call verify_inputs()


        !compute the lodi vector from the lodi outflow x
        detailled = .false.
        test_loc = test_lodi_outflow_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_outflow_x: '',L1)', test_loc
        print '()'


        !test the computation of the time derivatives for the lodi
        !outflow x
        detailled = .false.
        test_loc = test_lodi_outflow_timedev_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_outflow_timedevx: '',L1)', test_loc
        print '()'


        !compute the lodi vector from the lodi inflow x
        detailled = .false.
        test_loc = test_lodi_inflow_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_inflow_x: '',L1)', test_loc
        print '()'


        !test the computation of the time derivatives for the lodi
        !inflow x
        detailled = .false.
        test_loc = test_lodi_inflow_timedev_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_inflow_timedev_x: '',L1)', test_loc
        print '()'

        print '(''test_yoo_ns2d_edge_x: '',L1)', test_validated 
        print '(''----------------------------'')'

        contains

        !verify the inputs for the test
        subroutine verify_inputs()

          implicit none

          type(pmodel_eq) :: p_model


          if((nx.ne.5).or.(ny.ne.5).or.(ne.ne.4)) then
             print '(''test designed for:'')'
             print '(''nx=5'')'
             print '(''ny=5'')'
             print '(''pm_model=ns2d'')'
             stop 'change inputs'
          end if
        

          if(
     $         (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $         (.not.is_test_validated(mach_infty,0.2d0,detailled)).or.
     $         (.not.is_test_validated(sigma_P,0.25d0,detailled)).or.
     $         (.not.(flow_direction.eq.x_direction)).or.
     $         (.not.is_test_validated(p_model%get_mach_ux_infty(left),mach_infty,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_ux_infty(right),mach_infty,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_uy_infty(left),0.0d0,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_uy_infty(right),0.0d0,detailled))) then

             print '(''the test requires: '')'
             print '(''gamma=5/3'')'
             print '(''mach_infty=0.2'')'
             print '(''sigma_P=0.25'')'
             print '(''flow_direction=x-direction'')'
             print '(''ic_choice=peak'')'
             stop ''

          end if

        end subroutine verify_inputs

      
        function test_lodi_outflow_x(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: transverse_lodi
          real(rkind), dimension(nx,ny,ne)   :: viscous_lodi
          
          type(pmodel_eq)                    :: p_model
          type(lodi_edge_outflow)            :: outflow_bc

          call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)        
          call print_nodes(nodes,x_map,y_map,detailled)
          call get_test_data_for_lodi_outflow_x(test_data)
          
          call outflow_bc%ini()
          
          test_validated = test_lodi_x(
     $         test_data,
     $         outflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_outflow_x


        function test_lodi_outflow_timedev_x(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: transverse_lodi
          real(rkind), dimension(nx,ny,ne)   :: viscous_lodi
          
          type(pmodel_eq)                    :: p_model
          type(lodi_edge_outflow)            :: outflow_bc

          call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)        
          call get_test_data_for_lodi_outflow_timedevx(test_data)        
          call print_nodes(nodes,x_map,y_map,detailled)

          call outflow_bc%ini()
          
          test_validated = test_lodi_timedev_x(
     $         test_data,
     $         outflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_outflow_timedev_x


        function test_lodi_inflow_x(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: transverse_lodi
          real(rkind), dimension(nx,ny,ne)   :: viscous_lodi
          
          type(pmodel_eq)                    :: p_model
          type(lodi_edge_inflow)             :: inflow_bc

          call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)        
          call print_nodes(nodes,x_map,y_map,detailled)
          call get_test_data_for_lodi_inflow_x(test_data)
          
          call inflow_bc%ini()
          
          test_validated = test_lodi_x(
     $         test_data,
     $         inflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_inflow_x


        function test_lodi_inflow_timedev_x(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: transverse_lodi
          real(rkind), dimension(nx,ny,ne)   :: viscous_lodi

          type(pmodel_eq)                    :: p_model
          type(lodi_edge_inflow)             :: inflow_bc

          call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate(transverse_lodi, viscous_lodi)
          call get_test_data_for_lodi_inflow_timedevx(test_data)
          call print_nodes(nodes,x_map,y_map,detailled)

          call inflow_bc%ini()
          
          test_validated = test_lodi_timedev_x(
     $         test_data,
     $         inflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_inflow_timedev_x


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


        subroutine print_nodes(nodes,x_map,y_map,detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          logical                         , intent(in) :: detailled

          integer(ikind) :: j


          if(detailled) then
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
          end if

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
     $     viscous_lodi)
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
          detailled_loc = [.true.,.true.,.true.,.true.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,5
                do i=1,2
                   loc = is_test_validated(lodi(i,j,k),test_data(i,j,k),.false.)
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k).and.(.not.loc)) then
                      print '(''['',3I2,'']: '',F8.3,'' -> '', F8.3)',
     $                     i,j,k,
     $                     lodi(i,j,k),
     $                     test_data(i,j,k)
                   end if
                end do
             
                do i=4,5
                   loc = is_test_validated(lodi(i,j,k),test_data(i,j,k),.false.)
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k).and.(.not.loc)) then
                      print '(''['',3I2,'']: '',F8.3,'' -> '', F8.3)',
     $                     i,j,k,
     $                     lodi(i,j,k),
     $                     test_data(i,j,k)
                   end if
                end do
             end do
             if(.not.detailled_loc(k)) then
                print '(''lodi('',I1,''):'',L3)', k, test_lodi_validated
             end if
          end do

        end function test_lodi_x


        function test_lodi_timedev_x(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     transverse_lodi, viscous_lodi)
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
          detailled_loc = [.true.,.true.,.true.,.true.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,5
                do i=1,2
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),.false.)
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k).and.(.not.loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k, timedev(i,j,k), test_data(i,j,k)
                   end if
                end do
             
                do i=4,5
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),.false.)
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k).and.(.not.loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k, timedev(i,j,k), test_data(i,j,k)
                   end if
                end do
             end do
             if(.not.detailled_loc(k)) then
                print '(''timedev('',I1,''):'',L3)', k, test_lodi_validated
             end if
          end do

        end function test_lodi_timedev_x



      end program test_yoo_ns2d_edge_x
