      program test_openbc_operators

        use check_data_module, only :
     $       is_real_validated

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right,
     $       inflow_left,
     $       inflow_right,
     $       compute_fluxes_at_the_edges_2ndorder,
     $       add_body_forces

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       x_min,x_max,y_min,y_max
        
        use parameters_kind, only :
     $       ikind, rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq
        
        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0
        
        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_incoming_left(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_incoming_left: '',L1)', test_loc
        print '()'

        test_loc = test_incoming_right(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_incoming_right: '',L1)', test_loc
        print '()'

        test_loc = test_inflow_left(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_inflow_left: '',L1)', test_loc
        print '()'

        test_loc = test_inflow_right(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_inflow_right: '',L1)', test_loc
        print '()'

        test_loc = test_compute_fluxes_at_the_edges_2ndorder(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_fluxes_at_the_edges_2ndorder: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains

        function test_incoming_left(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(2) :: eigenvalues_test
          logical    , dimension(2) :: incoming_test
          logical                   :: incoming
          integer                   :: k

          test_validated = .true.

          eigenvalues_test = [1.0d0,-1.0d0]
          incoming_test    = [.true.,.false.]

          do k=1, size(eigenvalues_test,1)

             !output
             incoming = incoming_left(eigenvalues_test(k))

             !validation
             test_loc = incoming.eqv.incoming_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,''): '',L1,'' -> '',L1)', k, incoming, incoming_test(k)
             end if

          end do

        end function test_incoming_left


        function test_incoming_right(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(2) :: eigenvalues_test
          logical    , dimension(2) :: incoming_test
          logical                   :: incoming
          integer                   :: k

          test_validated = .true.

          eigenvalues_test = [1.0d0,-1.0d0]
          incoming_test    = [.false.,.true.]

          do k=1, size(eigenvalues_test,1)

             !output
             incoming = incoming_right(eigenvalues_test(k))

             !validation
             test_loc = incoming.eqv.incoming_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,''): '',L1,'' -> '',L1)', k, incoming, incoming_test(k)
             end if

          end do

        end function test_incoming_right


        function test_inflow_left(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(2) :: eigenvalues_test
          logical    , dimension(2) :: inflow_test
          logical                   :: inflow
          integer                   :: k

          test_validated = .true.

          eigenvalues_test = [1.0d0,-1.0d0]
          inflow_test    = [.true.,.false.]

          do k=1, size(eigenvalues_test,1)

             !output
             inflow = inflow_left(eigenvalues_test(k))

             !validation
             test_loc = inflow.eqv.inflow_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,''): '',L1,'' -> '',L1)', k, inflow, inflow_test(k)
             end if

          end do

        end function test_inflow_left


        function test_inflow_right(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(2) :: eigenvalues_test
          logical    , dimension(2) :: inflow_test
          logical                   :: inflow
          integer                   :: k

          test_validated = .true.

          eigenvalues_test = [1.0d0,-1.0d0]
          inflow_test    = [.false.,.true.]

          do k=1, size(eigenvalues_test,1)

             !output
             inflow = inflow_right(eigenvalues_test(k))

             !validation
             test_loc = inflow.eqv.inflow_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,''): '',L1,'' -> '',L1)', k, inflow, inflow_test(k)
             end if

          end do

          end function test_inflow_right


        function test_compute_fluxes_at_the_edges_2ndorder(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          type(sd_operators_x_oneside_L0)    :: s_x_L0
          type(sd_operators_x_oneside_L1)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    :: s_x_R1
          type(sd_operators_x_oneside_R0)    :: s_x_R0
          type(sd_operators_y_oneside_L0)    :: s_y_L0
          type(sd_operators_y_oneside_L1)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    :: s_y_R1
          type(sd_operators_y_oneside_R0)    :: s_y_R0
          type(pmodel_eq)                    :: p_model
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map

          integer(ikind) :: i
          integer(ikind) :: j
          integer        :: k

          test_validated = .true.
          

          !maps + nodes initialization
          dx    = (x_max-x_min)/(nx+1)
          dy    = (y_max-y_min)/(ny+1)
          x_map = (/ (x_min+(i-1)*dx,i=1,nx) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,ny) /)

          call p_model%apply_ic(nodes,x_map,y_map)


          !fluxes initialization
          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,nx+1), j=1,ny), k=1,ne) /),
     $         (/ nx+1,ny,ne /))
          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,nx), j=1,ny+1), k=1,ne) /),
     $         (/ nx,ny+1,ne /))

          !output
          call compute_fluxes_at_the_edges_2ndorder(
     $         nodes, dx, dy,
     $         s_x_L0, s_x_L1, s_x_R1, s_x_R0,
     $         s_y_L0, s_y_L1, s_y_R1, s_y_R0,
     $         p_model,
     $         flux_x, flux_y)

          !validation-flux_x
          do k=1, ne
             do j=1,bc_size
                do i=bc_size+1, nx-bc_size+1
                   test_loc = is_real_validated(flux_x(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.(.not.test_loc)
                   if(detailled.and.test_loc) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=ny-bc_size+1,ny
                do i=bc_size+1, nx-bc_size+1
                   test_loc = is_real_validated(flux_x(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.(.not.test_loc)
                   if(detailled.and.test_loc) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=bc_size+1,ny-bc_size
                do i=1, nx
                   test_loc = is_real_validated(flux_x(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=1,ny
                do i=1,bc_size
                   test_loc = is_real_validated(flux_x(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=1,ny
                do i=nx-bc_size+2,nx
                   test_loc = is_real_validated(flux_x(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          !validation-flux_y
          do k=1, ne
             do j=bc_size+1, ny-bc_size+1
                do i=1,bc_size
                   test_loc = is_real_validated(flux_y(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.(.not.test_loc)
                   if(detailled.and.test_loc) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=bc_size+1, ny-bc_size+1
                do i=nx-bc_size+1,nx
                   test_loc = is_real_validated(flux_y(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.(.not.test_loc)
                   if(detailled.and.test_loc) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=1,ny
                do i=bc_size+1,nx-bc_size
                   test_loc = is_real_validated(flux_y(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=1,bc_size
                do i=1,ny
                   test_loc = is_real_validated(flux_y(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

          do k=1, ne
             do j=ny-bc_size+2,ny
                do i=1,nx
                   test_loc = is_real_validated(flux_y(i,j,k),-99.0d0,.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print *, i,j,k
                   end if
                end do
             end do
          end do

        end function test_compute_fluxes_at_the_edges_2ndorder

      end program test_openbc_operators
