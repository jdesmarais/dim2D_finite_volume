      program test_lodi_operators

        use bc_operators_class, only :
     $       bc_operators

        use lodi_ns2d_class, only : 
     $       lodi_ns2d

        use lodi_inflow_class, only :
     $       lodi_inflow

        use lodi_outflow_class, only :
     $       lodi_outflow

        use parameters_constant, only :
     $       vector_x,
     $       left,
     $       right

        use parameters_input, only :
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
        real(rkind), dimension(nx,ny)      :: nodes_temp
        real(rkind), dimension(nx)         :: x_map
        real(rkind), dimension(ny)         :: y_map
        real(rkind)                        :: t
        real(rkind)                        :: dx
        real(rkind)                        :: dy
        !real(rkind), dimension(nx+1,ny,ne) :: flux_x
        !real(rkind), dimension(nx,ny+1,ne) :: flux_y
        !real(rkind), dimension(nx,ny,ne)   :: timedev
        type(pmodel_eq)                    :: p_model
        !type(bc_operators)                 :: bc_used
        type(lodi_outflow)                 :: outflow_bc
        type(lodi_inflow)                  :: inflow_bc

        character(*), parameter :: FMT='(5F14.5)'
        
        real(rkind), dimension(nx,ny,ne) :: test_data
        logical                          :: detailled
        logical                          :: test_validated

        integer :: k


        if((nx.ne.5).or.(ny.ne.5).or.(ne.ne.4)) then
           print '(''test designed for:'')'
           print '(''nx=5'')'
           print '(''ny=5'')'
           print '(''pm_model=ns2d'')'
           stop 'change inputs'
        end if


        !compute the lodi vector from the lodi outflow x
        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call print_nodes(nodes,x_map,y_map)
        call get_test_data_for_lodi_outflow_x(test_data)

        print '(''test lodi for outflow x'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        call outflow_bc%ini()

        if(.not.is_test_validated(
     $       outflow_bc%get_relaxation_P(),
     $       0.1249875d0,.false.)) then
           stop 'the test requires relaxation_P=0.1249875'
        end if

        test_validated = test_lodi_x(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       detailled)

        print '()'


        !test the computation of the time derivatives for the lodi
        !outflow x

        print '(''test time_dev for outflow x'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call get_test_data_for_lodi_outflow_timedevx(test_data)

        detailled = .false.

        call outflow_bc%ini()

        if(.not.is_test_validated(
     $       outflow_bc%get_relaxation_P(),
     $       0.1249875d0,.false.)) then
           stop 'the test requires relaxation_P=0.1249875'
        end if

        test_validated = test_lodi_timedev_x(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       detailled)

        print '()'


        !compute the lodi vector from the lodi outflow y
        print '(''test lodi for outflow y'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        do k=1,ne
           nodes(:,:,k) = TRANSPOSE(nodes(:,:,k))
        end do
        nodes_temp   = nodes(:,:,2)
        nodes(:,:,2) = nodes(:,:,3)
        nodes(:,:,3) = nodes_temp
        

        call get_test_data_for_lodi_outflow_x(test_data)
        do k=1, ne
           test_data(:,:,k) = TRANSPOSE(test_data(:,:,k))
        end do

        detailled = .false.

        call outflow_bc%ini()

        if(.not.is_test_validated(
     $       outflow_bc%get_relaxation_P(),
     $       0.1249875d0,.false.)) then
           stop 'the test requires relaxation_P=0.1249875'
        end if

        test_validated = test_lodi_y(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, y_map, x_map,
     $       detailled)
        print '()'


        !test the computation of the time derivatives for the lodi
        !outflow y

        print '(''test time_dev for outflow y'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        do k=1,ne
           nodes(:,:,k) = TRANSPOSE(nodes(:,:,k))
        end do
        nodes_temp   = nodes(:,:,2)
        nodes(:,:,2) = nodes(:,:,3)
        nodes(:,:,3) = nodes_temp

        call get_test_data_for_lodi_outflow_timedevx(test_data)
        do k=1, ne
           test_data(:,:,k) = TRANSPOSE(test_data(:,:,k))
        end do
        nodes_temp       = test_data(:,:,2)
        test_data(:,:,2) = test_data(:,:,3)
        test_data(:,:,3) = nodes_temp

        detailled = .false.

        call outflow_bc%ini()

        if(.not.is_test_validated(
     $       outflow_bc%get_relaxation_P(),
     $       0.1249875d0,.false.)) then
           stop 'the test requires relaxation_P=0.1249875'
        end if

        test_validated = test_lodi_timedev_y(
     $       test_data,
     $       outflow_bc,
     $       p_model,
     $       t, nodes, y_map, x_map,
     $       detailled)

        print '()'


        !compute the lodi vector from the lodi inflow x
        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        call get_test_data_for_lodi_inflow_x(test_data)

        print '(''test lodi for inflow x'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        call inflow_bc%ini()

        print '(''*************************************************'')'
        print '(''test only valid (duin/dt,dvin/dt,dTin/dt)=(0,0,0)'')'
        print '(''*************************************************'')'

        test_validated = test_lodi_x(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
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

        print '(''*************************************************'')'
        print '(''test only valid (duin/dt,dvin/dt,dTin/dt)=(0,0,0)'')'
        print '(''*************************************************'')'

        test_validated = test_lodi_timedev_x(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       detailled)

        print '()'


        !compute the lodi vector from the lodi inflow x
        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        do k=1,ne
           nodes(:,:,k) = TRANSPOSE(nodes(:,:,k))
        end do
        nodes_temp   = nodes(:,:,2)
        nodes(:,:,2) = nodes(:,:,3)
        nodes(:,:,3) = nodes_temp

        call get_test_data_for_lodi_inflow_x(test_data)
        do k=1, ne
           test_data(:,:,k) = TRANSPOSE(test_data(:,:,k))
        end do
        nodes_temp       = test_data(:,:,2)
        test_data(:,:,2) = test_data(:,:,3)
        test_data(:,:,3) = nodes_temp

        print '(''test lodi for inflow y'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        call inflow_bc%ini()

        print '(''*************************************************'')'
        print '(''test only valid (duin/dt,dvin/dt,dTin/dt)=(0,0,0)'')'
        print '(''*************************************************'')'

        test_validated = test_lodi_y(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, y_map, x_map,
     $       detailled)

        print '()'


        !test the computation of the time derivatives for the lodi
        !inflow x

        print '(''test time_dev for inflow y'')'
        print '(''---------------------------------------'')'

        call initialize_nodes(p_model,nodes,x_map,y_map,dx,dy)
        do k=1,ne
           nodes(:,:,k) = TRANSPOSE(nodes(:,:,k))
        end do
        nodes_temp   = nodes(:,:,2)
        nodes(:,:,2) = nodes(:,:,3)
        nodes(:,:,3) = nodes_temp

        call get_test_data_for_lodi_inflow_timedevx(test_data)
        do k=1, ne
           test_data(:,:,k) = TRANSPOSE(test_data(:,:,k))
        end do
        nodes_temp       = test_data(:,:,2)
        test_data(:,:,2) = test_data(:,:,3)
        test_data(:,:,3) = nodes_temp

        detailled = .false.

        call inflow_bc%ini()

        print '(''*************************************************'')'
        print '(''test only valid (duin/dt,dvin/dt,dTin/dt)=(0,0,0)'')'
        print '(''*************************************************'')'

        test_validated = test_lodi_timedev_y(
     $       test_data,
     $       inflow_bc,
     $       p_model,
     $       t, nodes, y_map, x_map,
     $       detailled)

        print '()'



c$$$        !test the symmetry for the inflow operators
c$$$        call bc_used%ini(p_model)
c$$$        call bc_used%apply_bc_on_timedev(
c$$$     $       p_model,
c$$$     $       t, nodes, x_map, y_map,
c$$$     $       flux_x, flux_y,
c$$$     $       timedev)
c$$$
c$$$        call print_timedev(timedev)

        


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
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_ns2d)                , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
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
     $            left,
     $            gradient_x_x_oneside_L0)
             
             i=2
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_x_x_oneside_L1)
             
             i=4
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_x_x_oneside_R1)

             i=5
             lodi(i,j,:) = bc_used%compute_x_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
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
                   test_validated = test_lodi_validated.and.loc
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


        function test_lodi_y(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_ns2d)                , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
          logical                         , intent(in)  :: detailled
          logical                                       :: test_validated


          real(rkind), dimension(nx,ny,ne) :: lodi
          logical                          :: loc
          logical                          :: test_lodi_validated
          logical, dimension(ne)           :: detailled_loc

          integer(ikind) :: i,j
          integer        :: k

          
          do i=1,5
             j=1
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_y_y_oneside_L0)
             
             j=2
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_y_y_oneside_L1)
             
             j=4
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_y_y_oneside_R1)

             j=5
             lodi(i,j,:) = bc_used%compute_y_lodi(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_y_y_oneside_R0)

          end do


          test_validated = .true.
          detailled_loc = [.false.,.false.,.false.,.false.]

          do k=1,4
             test_lodi_validated = .true.
             do i=1,5
                do j=1,2
                   loc = is_test_validated(lodi(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''lodi('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             
                do j=4,5
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

        end function test_lodi_y


        function test_lodi_timedev_x(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_ns2d)                , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
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
     $            left,
     $            gradient_x_x_oneside_L0)
             
             i=2
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_x_x_oneside_L1)
             
             i=4
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_x_x_oneside_R1)

             i=5
             timedev(i,j,:) = bc_used%compute_x_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
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
                   test_validated = test_lodi_validated.and.loc
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


        function test_lodi_timedev_y(
     $     test_data,
     $     bc_used,
     $     p_model,
     $     t, nodes, x_map, y_map,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in)  :: test_data
          class(lodi_ns2d)                , intent(in)  :: bc_used
          type(pmodel_eq)                 , intent(in)  :: p_model
          real(rkind)                     , intent(in)  :: t
          real(rkind), dimension(nx,ny,ne), intent(in)  :: nodes
          real(rkind), dimension(nx)      , intent(in)  :: x_map
          real(rkind), dimension(ny)      , intent(in)  :: y_map
          logical                         , intent(in)  :: detailled
          logical                                       :: test_validated


          real(rkind), dimension(nx,ny,ne) :: timedev
          logical                          :: loc
          logical                          :: test_lodi_validated
          logical, dimension(ne)           :: detailled_loc

          integer(ikind) :: i,j
          integer        :: k

          
          do i=1,5
             j=1
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_y_y_oneside_L0)
             
             j=2
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            left,
     $            gradient_y_y_oneside_L1)
             
             j=4
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_y_y_oneside_R1)

             j=5
             timedev(i,j,:) = bc_used%compute_y_timedev(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            right,
     $            gradient_y_y_oneside_R0)

          end do


          test_validated = .true.
          detailled_loc = [.false.,.false.,.false.,.false.]

          do k=1,4
             test_lodi_validated = .true.
             do i=1,5
                do j=1,2
                   loc = is_test_validated(timedev(i,j,k),test_data(i,j,k),detailled_loc(k))
                   test_validated = test_validated.and.loc
                   test_validated = test_lodi_validated.and.loc
                   if(detailled_loc(k)) then
                      print '(''timedev('',I2,I2,I2,''):'',L3)', i,j,k,loc
                   end if
                end do
             
                do j=4,5
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

        end function test_lodi_timedev_y


        subroutine get_test_data_for_lodi_outflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !L1
          test_data(1,5,1) =  446.77785
          test_data(1,4,1) = -1.2737843
          test_data(1,3,1) = 63.7416509
          test_data(1,2,1) = 81.1066497
          test_data(1,1,1) = -15.193722
                                       
          test_data(2,5,1) = -82.306674
          test_data(2,4,1) = -1.1667903
          test_data(2,3,1) = 69.9870196
          test_data(2,2,1) = -6.3174795
          test_data(2,1,1) = 10.5902500
                                       
          test_data(4,5,1) = -82.306674
          test_data(4,4,1) = -1.1667903
          test_data(4,3,1) = 69.9870196
          test_data(4,2,1) = -6.3174795
          test_data(4,1,1) = 10.5902500
                                       
          test_data(5,5,1) = 446.777854
          test_data(5,4,1) = -1.2737843
          test_data(5,3,1) = 63.7416509
          test_data(5,2,1) = 81.1066497
          test_data(5,1,1) = -15.193722

          !L2
          test_data(1,5,2) = 612.011517
          test_data(1,4,2) = 60.8978772
          test_data(1,3,2) = -1973.1928
          test_data(1,2,2) = -3957.7409
          test_data(1,1,2) = 3753.61415
                                       
          test_data(2,5,2) = 1648.38519
          test_data(2,4,2) = -187.99854
          test_data(2,3,2) = -501.81569
          test_data(2,2,2) = 29.5245331
          test_data(2,1,2) = 2106.98704
                                       
          test_data(4,5,2) = 1648.38519
          test_data(4,4,2) = -187.99854
          test_data(4,3,2) = -501.81569
          test_data(4,2,2) = 29.5245331
          test_data(4,1,2) = 2106.98704
                                       
          test_data(5,5,2) = 612.011517
          test_data(5,4,2) = 60.8978772
          test_data(5,3,2) = -1973.1928
          test_data(5,2,2) = -3957.7409
          test_data(5,1,2) = 3753.61415

          !L3
          test_data(1,5,3) = -3872.8478
          test_data(1,4,3) = 69.6020753
          test_data(1,3,3) = -2763.6543
          test_data(1,2,3) = 1973.37048
          test_data(1,1,3) = -6295.0801
                       
          test_data(2,5,3) = 1367.71567
          test_data(2,4,3) = -299.00369
          test_data(2,3,3) = 497.172700
          test_data(2,2,3) = -11.720210
          test_data(2,1,3) = -15463.56565d0
                       
          test_data(4,5,3) = 18.1041049
          test_data(4,4,3) = 0.45929891
          test_data(4,3,3) = 12.3627731
          test_data(4,2,3) = 0.25103836
          test_data(4,1,3) = 71.3446734
                       
          test_data(5,5,3) = 0.32684954
          test_data(5,4,3) = 4.19698320
          test_data(5,3,3) = 5.00151066
          test_data(5,2,3) = 3.25023657
          test_data(5,1,3) = 0.81473509

          !L4
          test_data(1,5,4) = 0.32684954
          test_data(1,4,4) = 4.19698320
          test_data(1,3,4) = 5.00151066
          test_data(1,2,4) = 3.25023657
          test_data(1,1,4) = 0.81473509
                       
          test_data(2,5,4) = 18.1041049
          test_data(2,4,4) = 0.45929891
          test_data(2,3,4) = 12.3627731
          test_data(2,2,4) = 0.25103836
          test_data(2,1,4) = 71.3446734
                       
          test_data(4,5,4) = 1367.71567
          test_data(4,4,4) = -299.00369
          test_data(4,3,4) = 497.172700
          test_data(4,2,4) = -11.720210
          test_data(4,1,4) = -15463.56565d0
                       
          test_data(5,5,4) = -3872.8478
          test_data(5,4,4) = 69.6020753
          test_data(5,3,4) = -2763.6543
          test_data(5,2,4) = 1973.37048
          test_data(5,1,4) = -6295.0801

        end subroutine get_test_data_for_lodi_outflow_x



        subroutine get_test_data_for_lodi_outflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data


          !mass
          test_data(1,5,1) =   43.95776926
          test_data(1,4,1) =  -14.42388585
          test_data(1,3,1) =   11.27966366 
          test_data(1,2,1) =   79.17163958 
          test_data(1,1,1) = -111.3175171

          test_data(2,5,1) =  -34.67459837
          test_data(2,4,1) =  337.6488755 
          test_data(2,3,1) =    6.676164817 
          test_data(2,2,1) =  -40.8029952 
          test_data(2,1,1) =    5.981911497 
                                      
          test_data(4,5,1) =  -34.6745983
          test_data(4,4,1) =  337.6488755
          test_data(4,3,1) =    6.676164817
          test_data(4,2,1) =  -40.8029952
          test_data(4,1,1) =    5.981911497
                                      
          test_data(5,5,1) =   43.95776926 
          test_data(5,4,1) =  -14.42388585
          test_data(5,3,1) =   11.27966366 
          test_data(5,2,1) =   79.17163958 
          test_data(5,1,1) = -111.3175171
                                 
          !momentum-x            
          test_data(1,5,2) = -623.1740058 
          test_data(1,4,2) =   -1.40660041 
          test_data(1,3,2) =  -87.6536906 
          test_data(1,2,2) = -254.8048462 
          test_data(1,1,2) =-1056.082304  

          test_data(2,5,2) =    6.993121655
          test_data(2,4,2) =  -19.4309979
          test_data(2,3,2) =   30.63116387
          test_data(2,2,2) =  -45.7003492
          test_data(2,1,2) = -206.462977
                                      
          test_data(4,5,2) =   -6.993121655   
          test_data(4,4,2) =   19.43099799    
          test_data(4,3,2) =  -30.63116387
          test_data(4,2,2) =   45.7003492     
          test_data(4,1,2) =  206.462977  
                                      
          test_data(5,5,2) =  623.1740058 
          test_data(5,4,2) =    1.40660041
          test_data(5,3,2) =   87.6536906 
          test_data(5,2,2) =  254.8048462 
          test_data(5,1,2) = 1056.082304
                                 
          !momentum-y            
          test_data(1,5,3) = 1251.861906   
          test_data(1,4,3) =    7.212710942 
          test_data(1,3,3) =  449.276891   
          test_data(1,2,3) =  373.413394   
          test_data(1,1,3) =  -62.33626751 
                                 
          test_data(2,5,3) =  245.9295404  
          test_data(2,4,3) =  -83.96891766 
          test_data(2,3,3) = -301.351884   
          test_data(2,2,3) =   63.86683771 
          test_data(2,1,3) =    8.257957095
                                 
          test_data(4,5,3) =  245.9295404  
          test_data(4,4,3) =  -83.96891766 
          test_data(4,3,3) = -301.351884   
          test_data(4,2,3) =   63.86683771 
          test_data(4,1,3) =    8.257957095
                                 
          test_data(5,5,3) = 1251.861906   
          test_data(5,4,3) =    7.212710942
          test_data(5,3,3) =  449.276891   
          test_data(5,2,3) =  373.413394   
          test_data(5,1,3) =  -62.33626751


          !total energy
          test_data(1,5,4) = 23640.52055d0 
          test_data(1,4,4) =   -47.6898741
          test_data(1,3,4) = 11061.73524 
          test_data(1,2,4) =  -415.045323
          test_data(1,1,4) =  7868.867401d0
                                  
          test_data(2,5,4) =  -548.9516464
          test_data(2,4,4) =   201.54906   
          test_data(2,3,4) = -1058.056021
          test_data(2,2,4) =   -30.050448  
          test_data(2,1,4) =  9705.71074d0  
                                  
          test_data(4,5,4) =  -548.9516464
          test_data(4,4,4) =   201.54906   
          test_data(4,3,4) = -1058.056021
          test_data(4,2,4) =   -30.050448  
          test_data(4,1,4) =  9705.71074d0  
                                  
          test_data(5,5,4) = 23640.52055d0
          test_data(5,4,4) =   -47.6898741
          test_data(5,3,4) = 11061.73524 
          test_data(5,2,4) =  -415.045323
          test_data(5,1,4) =  7868.867401d0

        end subroutine get_test_data_for_lodi_outflow_timedevx


        subroutine get_test_data_for_lodi_inflow_x(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !L1
          test_data(1,5,1) = 0.0d0
          test_data(1,4,1) = 0.0d0 
          test_data(1,3,1) = 0.0d0 
          test_data(1,2,1) = 0.0d0 
          test_data(1,1,1) = 0.0d0 
                             
          test_data(2,5,1) = 0.0d0 
          test_data(2,4,1) = 0.0d0 
          test_data(2,3,1) = 0.0d0 
          test_data(2,2,1) = 0.0d0 
          test_data(2,1,1) = 0.0d0 
                             
          test_data(4,5,1) = 0.0d0 
          test_data(4,4,1) = 0.0d0 
          test_data(4,3,1) = 0.0d0 
          test_data(4,2,1) = 0.0d0 
          test_data(4,1,1) = 0.0d0 
                             
          test_data(5,5,1) = 0.0d0 
          test_data(5,4,1) = 0.0d0 
          test_data(5,3,1) = 0.0d0 
          test_data(5,2,1) = 0.0d0 
          test_data(5,1,1) = 0.0d0 

          !L2
          test_data(1,5,2) = -2581.898576   
          test_data(1,4,2) =    46.40138354 
          test_data(1,3,2) = -1842.436265  
          test_data(1,2,2) =  1315.58032   
          test_data(1,1,2) = -4196.720112  
                                      
          test_data(2,5,2) =   911.8104509  
          test_data(2,4,2) =  -199.3357957  
          test_data(2,3,2) =   331.4484672  
          test_data(2,2,2) =    -7.813473556
          test_data(2,1,2) =-10309.04377d0    
                                      
          test_data(4,5,2) =   911.8104509 
          test_data(4,4,2) =  -199.3357957
          test_data(4,3,2) =   331.4484672 
          test_data(4,2,2) =    -7.813473556
          test_data(4,1,2) =-10309.04377d0
                                      
          test_data(5,5,2) = -2581.898576  
          test_data(5,4,2) =    46.40138354
          test_data(5,3,2) = -1842.436265  
          test_data(5,2,2) =  1315.58032   
          test_data(5,1,2) = -4196.720112  

          !L3
          test_data(1,5,3) = -3872.847865  
          test_data(1,4,3) =    69.60207531 
          test_data(1,3,3) = -2763.654398  
          test_data(1,2,3) =  1973.37048   
          test_data(1,1,3) = -6295.080168  
                                            
          test_data(2,5,3) =  1367.715676  
          test_data(2,4,3) =  -299.0036936 
          test_data(2,3,3) =   497.1727008 
          test_data(2,2,3) =   -11.72021033
          test_data(2,1,3) =-15463.56565d0   
                                           
          test_data(4,5,3) =  1367.715676  
          test_data(4,4,3) =  -299.0036936 
          test_data(4,3,3) =   497.1727008 
          test_data(4,2,3) =   -11.72021033
          test_data(4,1,3) =-15463.56565d0   
                                           
          test_data(5,5,3) = -3872.847865  
          test_data(5,4,3) =    69.60207531
          test_data(5,3,3) = -2763.654398  
          test_data(5,2,3) =  1973.37048   
          test_data(5,1,3) = -6295.080168  

          !L4
          test_data(1,5,4) = -3872.847865  
          test_data(1,4,4) =    69.60207531
          test_data(1,3,4) = -2763.654398  
          test_data(1,2,4) =  1973.37048   
          test_data(1,1,4) = -6295.080168  
                                            
          test_data(2,5,4) =  1367.715676  
          test_data(2,4,4) =  -299.0036936 
          test_data(2,3,4) =   497.1727008 
          test_data(2,2,4) =   -11.72021033
          test_data(2,1,4) =-15463.56565d0   
                                           
          test_data(4,5,4) =  1367.715676  
          test_data(4,4,4) =  -299.0036936 
          test_data(4,3,4) =   497.1727008 
          test_data(4,2,4) =   -11.72021033
          test_data(4,1,4) =-15463.56565d0   
                                           
          test_data(5,5,4) = -3872.847865  
          test_data(5,4,4) =    69.60207531
          test_data(5,3,4) = -2763.654398  
          test_data(5,2,4) =  1973.37048   
          test_data(5,1,4) = -6295.080168   

        end subroutine get_test_data_for_lodi_inflow_x


        subroutine get_test_data_for_lodi_inflow_timedevx(
     $     test_data)
        
          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !mass
          test_data(1,5,1) = -214.2620132 
          test_data(1,4,1) =   17.1090493 
          test_data(1,3,1) =  -15.4973465 
          test_data(1,2,1) =   87.69076208
          test_data(1,1,1) =-1925.732738  
                                 
          test_data(2,5,1) =   33.75979961
          test_data(2,4,1) = -498.8982123 
          test_data(2,3,1) =   22.3924598 
          test_data(2,2,1) =  -33.50292477
          test_data(2,1,1) =  -27.58383676
                                     
          test_data(4,5,1) =   33.75979961
          test_data(4,4,1) = -498.8982123 
          test_data(4,3,1) =   22.3924598 
          test_data(4,2,1) =  -33.50292477
          test_data(4,1,1) =  -27.58383676
                                     
          test_data(5,5,1) = -214.2620132 
          test_data(5,4,1) =   17.1090493 
          test_data(5,3,1) =  -15.4973465 
          test_data(5,2,1) =   87.69076208
          test_data(5,1,1) =-1925.732738  
                    
          !momentum-x            
          test_data(1,5,2) = 1317.711381  
          test_data(1,4,2) =   16.56558538
          test_data(1,3,2) =   10.10696511
          test_data(1,2,2) = -460.3765009 
          test_data(1,1,2) = 5061.328     
                                 
          test_data(2,5,2) =   73.14623249
          test_data(2,4,2) = -192.6514635 
          test_data(2,3,2) =  -30.91652283
          test_data(2,2,2) =  -31.08759763
          test_data(2,1,2) = -219.7241899 
                            
          test_data(4,5,2) =   73.14623249
          test_data(4,4,2) = -192.6514635 
          test_data(4,3,2) =  -30.91652283
          test_data(4,2,2) =  -31.08759763
          test_data(4,1,2) = -219.7241899 
                            
          test_data(5,5,2) = 1317.711381  
          test_data(5,4,2) =   16.56558538
          test_data(5,3,2) =   10.10696511
          test_data(5,2,2) = -460.3765009 
          test_data(5,1,2) = 5061.328  

          !momentum-y            
          test_data(1,5,3) = -6537.455417 
          test_data(1,4,3) =     4.2873264
          test_data(1,3,3) =  -637.4125993
          test_data(1,2,3) =   521.3946562
          test_data(1,1,3) = -1682.922958 
                                  
          test_data(2,5,3) =    49.0454866
          test_data(2,4,3) =   137.5168   
          test_data(2,3,3) =    45.5810959
          test_data(2,2,3) =     7.83033  
          test_data(2,1,3) =   -87.8896759
                                 
          test_data(4,5,3) =   49.0454866
          test_data(4,4,3) =  137.5168    
          test_data(4,3,3) =   45.5810959
          test_data(4,2,3) =    7.83033   
          test_data(4,1,3) =  -87.88967594
                                 
          test_data(5,5,3) = -6537.455417 
          test_data(5,4,3) =     4.2873264
          test_data(5,3,3) =  -637.4125993
          test_data(5,2,3) =   521.3946562
          test_data(5,1,3) = -1682.922958 


          !total energy
          test_data(1,5,4) = -109595.0198d0 
          test_data(1,4,4) =     112.95998
          test_data(1,3,4) =  -17257.30602d0
          test_data(1,2,4) =    5718.60689
          test_data(1,1,4) =  -16829.22958d0
                                    
          test_data(2,5,4) =    2166.44136
          test_data(2,4,4) =    -504.65473
          test_data(2,3,4) =     813.49318
          test_data(2,2,4) =     -32.91857
          test_data(2,1,4) =  -24210.49578d0
                                  
          test_data(4,5,4) =    2166.44136
          test_data(4,4,4) =    -504.65473
          test_data(4,3,4) =     813.49318
          test_data(4,2,4) =     -32.91857
          test_data(4,1,4) =  -24210.49578d0
                                  
          test_data(5,5,4) = -109595.0198d0
          test_data(5,4,4) =     112.95998
          test_data(5,3,4) =  -17257.30602d0
          test_data(5,2,4) =    5718.60689
          test_data(5,1,4) =  -16829.22958d0

          test_data(:,:,1) = -test_data(:,:,1)
          test_data(:,:,3) = -test_data(:,:,3)
          test_data(:,:,4) = -test_data(:,:,4)

        end subroutine get_test_data_for_lodi_inflow_timedevx

      end program test_lodi_operators
