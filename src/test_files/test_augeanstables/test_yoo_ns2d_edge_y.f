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
     $       y_direction,
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
        test_loc = test_lodi_outflow_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_outflow_y: '',L1)', test_loc
        print '()'


        !test the computation of the time derivatives for the lodi
        !outflow x
        detailled = .false.
        test_loc = test_lodi_outflow_timedev_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_outflow_timedev_y: '',L1)', test_loc
        print '()'


        !compute the lodi vector from the lodi inflow x
        detailled = .false.
        test_loc = test_lodi_inflow_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_inflow_y: '',L1)', test_loc
        print '()'


        !test the computation of the time derivatives for the lodi
        !inflow x
        detailled = .false.
        test_loc = test_lodi_inflow_timedev_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_lodi_inflow_timedev_y: '',L1)', test_loc
        print '()'

        print '(''test_yoo_ns2d_edge_y: '',L1)', test_validated
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
     $         (.not.(flow_direction.eq.y_direction)).or.
     $         (.not.is_test_validated(p_model%get_mach_uy_infty(left),mach_infty,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_uy_infty(right),mach_infty,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_ux_infty(left),0.0d0,detailled)).or.
     $         (.not.is_test_validated(p_model%get_mach_ux_infty(right),0.0d0,detailled))) then

             print '(''the test requires: '')'
             print '(''gamma=5/3'')'
             print '(''mach_infty=0.2'')'
             print '(''sigma_P=0.25'')'
             print '(''flow_direction=y_direction'')'
             print '(''ic_choice=peak'')'
             stop ''

          end if

        end subroutine verify_inputs


        subroutine initialize_nodes_for_test(p_model,nodes,x_map,y_map,dx,dy)

          implicit none

          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy


          integer                       :: k
          real(rkind), dimension(nx,ny) :: tmp_field
          real(rkind), dimension(nx)    :: tmp_map
          real(rkind)                   :: tmp_ds


          call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)

          do k=1,ne
             nodes(:,:,k) = transpose(nodes(:,:,k))
          end do

          tmp_field    = nodes(:,:,2)
          nodes(:,:,2) = nodes(:,:,3)
          nodes(:,:,3) = tmp_field

          tmp_map = x_map
          x_map   = y_map
          y_map   = tmp_map

          tmp_ds = dx
          dx     = dy
          dy     = tmp_ds

        end subroutine initialize_nodes_for_test

      
        function test_lodi_outflow_y(detailled)
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

          call initialize_nodes_for_test(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate_for_test(transverse_lodi, viscous_lodi)        
          call print_nodes(nodes,x_map,y_map,detailled)
          call get_test_data_for_lodi_outflow_y(test_data)
          
          call outflow_bc%ini()
          
          test_validated = test_lodi_y(
     $         test_data,
     $         outflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_outflow_y


        function test_lodi_outflow_timedev_y(detailled)
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

          call initialize_nodes_for_test(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate_for_test(transverse_lodi, viscous_lodi)
          call get_test_data_for_lodi_outflow_timedevy(test_data)
          call print_nodes(nodes,x_map,y_map,detailled)

          call outflow_bc%ini()
          
          test_validated = test_lodi_timedev_y(
     $         test_data,
     $         outflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_outflow_timedev_y


        function test_lodi_inflow_y(detailled)
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

          call initialize_nodes_for_test(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate_for_test(transverse_lodi, viscous_lodi)        
          call print_nodes(nodes,x_map,y_map,detailled)
          call get_test_data_for_lodi_inflow_y(test_data)
          
          call inflow_bc%ini()
          
          test_validated = test_lodi_y(
     $         test_data,
     $         inflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_inflow_y


        function test_lodi_inflow_timedev_y(detailled)
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

          call initialize_nodes_for_test(p_model,nodes,x_map,y_map,dx,dy)
          call initialize_lodi_intermediate_for_test(transverse_lodi, viscous_lodi)
          call get_test_data_for_lodi_inflow_timedevy(test_data)
          call print_nodes(nodes,x_map,y_map,detailled)

          call inflow_bc%ini()
          
          test_validated = test_lodi_timedev_y(
     $         test_data,
     $         inflow_bc,
     $         p_model,
     $         t, nodes, x_map, y_map,
     $         transverse_lodi, viscous_lodi)
          
        end function test_lodi_inflow_timedev_y


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


        subroutine initialize_lodi_intermediate_for_test(transverse_lodi, viscous_lodi)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: transverse_lodi
          real(rkind), dimension(nx,ny,ne), intent(out) :: viscous_lodi

          call initialize_lodi_intermediate(transverse_lodi,viscous_lodi)
          call modify_data_for_edge_y(transverse_lodi)
          call modify_data_for_edge_y(viscous_lodi)

       end subroutine initialize_lodi_intermediate_for_test


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


        function test_lodi_y(
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

          
          j=1
          do i=1,5
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            left,
     $            gradient_y_y_oneside_L0)
          end do

          
          j=2
          do i=1,5
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            left,
     $            gradient_y_y_oneside_L1)
          end do


          j=4
          do i=1,5
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            right,
     $            gradient_y_y_oneside_R1)
          end do


          j=5
          do i=1,5
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi(i,j,:), viscous_lodi(i,j,:),
     $            right,
     $            gradient_y_y_oneside_R0)
          end do


          test_validated = .true.
          detailled_loc = [.true.,.true.,.true.,.true.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,2
                do i=1,5
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

             do j=4,5
                do i=1,5
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

        end function test_lodi_y


        function test_lodi_timedev_y(
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

          
          j=1
          do i=1,5
             
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            left,
     $            gradient_y_y_oneside_L0)
             
          end do


          j=2
          do i=1,5
             
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            left,
     $            gradient_y_y_oneside_L1)
             
          end do


          j=4
          do i=1,5
             
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            right,
     $            gradient_y_y_oneside_R1)
             
          end do


          j=5
          do i=1,5
             
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            transverse_lodi, viscous_lodi,
     $            right,
     $            gradient_y_y_oneside_R0)
             
          end do


          test_validated = .true.
          detailled_loc = [.true.,.true.,.true.,.true.]

          do k=1,4
             test_lodi_validated = .true.
             do j=1,2
                do i=1,5
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),.false.)
                   test_validated = test_validated.and.loc
                   test_lodi_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k).and.(.not.loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k, timedev(i,j,k), test_data(i,j,k)
                   end if
                end do
             end do

             do j=4,5
                do i=1,5
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

        end function test_lodi_timedev_y


        !lodi_outflow_y
        subroutine get_test_data_for_lodi_outflow_y(test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          call get_test_data_for_lodi_outflow_x(test_data)
          call modify_data_for_edge_y(test_data,var=.false.)

        end subroutine get_test_data_for_lodi_outflow_y


        !lodi_outflow_timedevy
        subroutine get_test_data_for_lodi_outflow_timedevy(test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          call get_test_data_for_lodi_outflow_timedevx(test_data)
          call modify_data_for_edge_y(test_data)

        end subroutine get_test_data_for_lodi_outflow_timedevy


        !lodi_inflow_y
        subroutine get_test_data_for_lodi_inflow_y(test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          call get_test_data_for_lodi_inflow_x(test_data)
          call modify_data_for_edge_y(test_data,var=.false.)

        end subroutine get_test_data_for_lodi_inflow_y


        !lodi_inflow_timedevy
        subroutine get_test_data_for_lodi_inflow_timedevy(test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          call get_test_data_for_lodi_inflow_timedevx(test_data)
          call modify_data_for_edge_y(test_data)

        end subroutine get_test_data_for_lodi_inflow_timedevy

       
        !transpose the test data and exchange momentum-x and momentum-y
        subroutine modify_data_for_edge_y(data_modified,var)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: data_modified
          logical    , optional           , intent(in)    :: var

          integer                          :: k
          real(rkind), dimension(nx,ny)    :: tmp_data
          logical :: var_modified

          if(present(var)) then
             var_modified = var
          else
             var_modified = .true.
          end if


          do k=1, ne
             data_modified(:,:,k) = transpose(data_modified(:,:,k))
          end do
          
          if(var_modified) then
             tmp_data             = data_modified(:,:,2)
             data_modified(:,:,2) = data_modified(:,:,3)
             data_modified(:,:,3) = tmp_data
          end if

        end subroutine modify_data_for_edge_y

      end program test_yoo_ns2d_edge_x
