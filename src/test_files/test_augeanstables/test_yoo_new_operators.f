      program test_yoo_new_operators

        use bc_operators_class, only :
     $       bc_operators

        use bc_operators_ns2d_class, only :
     $       bc_operators_ns2d

        use ns2d_parameters, only :
     $       viscous_r,
     $       Re,
     $       Pr,
     $       gamma,
     $       mach_infty

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

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


        type(bc_operators_ns2d) :: bc_working
        type(bc_operators)      :: bc_bugged

        logical :: test_validated
        logical :: test_loc
        logical :: detailled
        
        
        detailled = .false.
        if(
     $       (.not.is_test_validated(viscous_r,-2.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(Re,10.0d0,detailled)).or.
     $       (.not.is_test_validated(Pr,1.0d0,detailled)).or.
     $       (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(mach_infty,1.0d0,detailled))) then

           print '(''the test requires: '')'
           print '(''viscous_r=-2/3'')'
           print '(''Re=10.0'')'
           print '(''Pr=1.0'')'
           print '(''gamma=5/3'')'
           print '(''mach_infty=1.0'')'
           stop ''

        end if

        if(
     $       (.not.(nx.eq.5)).or.
     $       (.not.(ny.eq.5))) then

           print '(''the test requires: '')'
           print '(''nx=5'')'
           print '(''ny=5'')'
           stop ''

        end if


        test_validated = .true.

        detailled      = .false.

        !test of the edge fluxes and the LODI terms
        !computation
        test_loc = test_compute_edge_fluxes(bc_bugged, bc_working, detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_edge_fluxes: '',L1)', test_loc


        !test if the two operators gives the same results for
        !the time derivatives
        test_loc = test_apply_bc_on_timedev(bc_bugged, bc_working, detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev: '',L1)', test_loc


        contains


        !test the function bc_used%compute_edge_fluxes
        function test_compute_edge_fluxes(bc_bugged, bc_working, detailled)
     $       result(test_validated)

          implicit none

          class(bc_operators)     , intent(in) :: bc_bugged
          class(bc_operators_ns2d), intent(in) :: bc_working
          logical                 , intent(in) :: detailled
          logical                              :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          logical                       :: test_loc

          type(pmodel_eq)               :: p_model

          real(rkind), dimension(6,5,4) :: flux_x_b
          real(rkind), dimension(5,6,4) :: flux_y_b
          real(rkind), dimension(6,5,4) :: flux_x_w
          real(rkind), dimension(5,6,4) :: flux_y_w

          !bugged values
          real(rkind), dimension(1,2,4) :: transverse_lodi_N_b
          real(rkind), dimension(1,2,4) :: transverse_lodi_S_b
          real(rkind), dimension(2,1,4) :: transverse_lodi_E_b
          real(rkind), dimension(2,1,4) :: transverse_lodi_W_b

          real(rkind), dimension(1,2,4) :: viscous_lodi_N_b
          real(rkind), dimension(1,2,4) :: viscous_lodi_S_b
          real(rkind), dimension(2,1,4) :: viscous_lodi_E_b
          real(rkind), dimension(2,1,4) :: viscous_lodi_W_b

          !working values
          real(rkind), dimension(1,2,4) :: transverse_lodi_N_w
          real(rkind), dimension(1,2,4) :: transverse_lodi_S_w
          real(rkind), dimension(2,1,4) :: transverse_lodi_E_w
          real(rkind), dimension(2,1,4) :: transverse_lodi_W_w

          real(rkind), dimension(1,2,4) :: viscous_lodi_N_w
          real(rkind), dimension(1,2,4) :: viscous_lodi_S_w
          real(rkind), dimension(2,1,4) :: viscous_lodi_E_w
          real(rkind), dimension(2,1,4) :: viscous_lodi_W_w


          integer(ikind) :: i,j
          integer        :: k

          logical :: test_fluxes
          logical :: test_lodi

          test_validated = .true.

          
          !initialize the nodes
          call initialize_nodes(nodes,dx,dy)

          !compute the fluxes for the two boundary conditions
          call bc_bugged%compute_edge_fluxes(
     $         p_model,
     $         nodes,dx,dy,
     $         transverse_lodi_N_b, transverse_lodi_S_b,
     $         transverse_lodi_E_b, transverse_lodi_W_b,
     $         viscous_lodi_N_b, viscous_lodi_S_b,
     $         viscous_lodi_E_b, viscous_lodi_W_b,
     $         flux_x_b, flux_y_b)

          call bc_working%compute_edge_fluxes(
     $         nodes,dx,dy,
     $         transverse_lodi_N_w, transverse_lodi_S_w,
     $         transverse_lodi_E_w, transverse_lodi_W_w,
     $         viscous_lodi_N_w, viscous_lodi_S_w,
     $         viscous_lodi_E_w, viscous_lodi_W_w,
     $         flux_x_w, flux_y_w)


          !test if the y-fluxes are the same
          do k=1,ne
             do j=3,4

                do i=1,2
                   test_loc = is_test_validated(
     $                  flux_y_b(i,j,k),
     $                  flux_y_w(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(.not.test_loc) then
                      print '(''different y-flux at '', 3I2, 2F8.5)',
     $                     i,j,k,
     $                     flux_y_b(i,j,k),
     $                     flux_y_w(i,j,k)
                   end if
                end do

                do i=4,5
                   test_loc = is_test_validated(
     $                  flux_y_b(i,j,k),
     $                  flux_y_w(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(.not.test_loc) then
                      print '(''different y-flux at '', 3I2, 2F8.5)',
     $                     i,j,k,
     $                     flux_y_b(i,j,k),
     $                     flux_y_w(i,j,k)
                   end if
                end do

             end do
          end do


          !test if the x-fluxes are the same
          do k=1,ne

             do j=1,2
                do i=3,4
                   test_loc = is_test_validated(
     $                  flux_x_b(i,j,k),
     $                  flux_x_w(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(.not.test_loc) then
                      print '(''different x-flux at '', 3I2)',i,j,k
                   end if
                end do
             end do

             do j=4,5
                do i=3,4
                   test_loc = is_test_validated(
     $                  flux_x_b(i,j,k),
     $                  flux_x_w(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(.not.test_loc) then
                      print '(''different x-flux at '', 3I2)',i,j,k
                   end if
                end do
             end do

          end do

          test_fluxes = test_validated

          if(.not.detailled) then
             print '(''fluxes: '',L1)', test_fluxes
          end if
          
          
          test_validated = .true.

          !test if the lodi terms are the same

          !East
          do k=1, ne
             j=1
             do i=1,2
                test_loc = is_test_validated(
     $               transverse_lodi_E_b(i,j,k),
     $               transverse_lodi_E_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                 print '(''different trans_lodi_E at '', 3I2,2F10.5)',
     $                  i,j,k,
     $                  transverse_lodi_E_b(i,j,k),
     $                  transverse_lodi_E_w(i,j,k)
                end if
             end do
          end do

          do k=1,ne
             j=1
             do i=1,2
                test_loc = is_test_validated(
     $               viscous_lodi_E_b(i,j,k),
     $               viscous_lodi_E_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different viscous_lodi_E at '', 3I2,2F10.5)',
     $                  i,j,k,
     $                  viscous_lodi_E_b(i,j,k),
     $                  viscous_lodi_E_w(i,j,k)
                end if
             end do
          end do


          !West
          do k=1, ne
             j=1
             do i=1,2
                test_loc = is_test_validated(
     $               transverse_lodi_W_b(i,j,k),
     $               transverse_lodi_W_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different trans_lodi_W at '', 3I2,2F10.5)', 
     $                  i,j,k,
     $                  transverse_lodi_W_b(i,j,k), 
     $                  transverse_lodi_W_w(i,j,k) 
                end if
             end do
          end do

          do k=1,ne
             j=1
             do i=1,2
                test_loc = is_test_validated(
     $               viscous_lodi_W_b(i,j,k),
     $               viscous_lodi_W_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different viscous_lodi_W at '', 3I2,2F10.5)',
     $                  i,j,k,
     $                  viscous_lodi_W_b(i,j,k),
     $                  viscous_lodi_W_w(i,j,k)
                end if
             end do
          end do


          !North
          do k=1, ne
             i=1
             do j=1,2
                test_loc = is_test_validated(
     $               transverse_lodi_N_b(i,j,k),
     $               transverse_lodi_N_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different trans_lodi_N at '', 3I2, 2F10.5)',
     $                  i,j,k,
     $                  transverse_lodi_N_b(i,j,k),
     $                  transverse_lodi_N_w(i,j,k)
                end if
             end do
          end do

          do k=1,ne
             i=1
             do j=1,2
                test_loc = is_test_validated(
     $               viscous_lodi_N_b(i,j,k),
     $               viscous_lodi_N_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different viscous_lodi_N at '', 3I2,2F10.5)',
     $                  i,j,k,
     $                  viscous_lodi_N_b(i,j,k),
     $                  viscous_lodi_N_w(i,j,k)
                end if
             end do
          end do


          !South
          do k=1, ne
             i=1
             do j=1,2
                test_loc = is_test_validated(
     $               transverse_lodi_S_b(i,j,k),
     $               transverse_lodi_S_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated
                if(.not.test_loc) then
                print '(''different trans_lodi_S at '', 3I2,2F10.5)',
     $                  i,j,k,
     $                  transverse_lodi_S_b(i,j,k),
     $                  transverse_lodi_S_w(i,j,k)
                end if
             end do
          end do

          do k=1,ne
             i=1
             do j=1,2
                test_loc = is_test_validated(
     $               viscous_lodi_S_b(i,j,k),
     $               viscous_lodi_S_w(i,j,k),
     $               detailled)
                test_validated = test_loc.and.test_validated 
                if(.not.test_loc) then 
                print '(''different viscous_lodi_S at '', 3I2,2F10.5)', 
     $                  i,j,k,
     $                  viscous_lodi_S_b(i,j,k), 
     $                  viscous_lodi_S_w(i,j,k)
                end if
             end do
          end do

          test_lodi = test_validated

          test_validated = test_fluxes.and.test_lodi

        end function test_compute_edge_fluxes


        !test the function bc_used%compute_edge_fluxes
        function test_apply_bc_on_timedev(bc_bugged, bc_working, detailled)
     $       result(test_validated)

          implicit none

          class(bc_operators)     , intent(inout) :: bc_bugged
          class(bc_operators_ns2d), intent(inout) :: bc_working
          logical                 , intent(in)    :: detailled
          logical                                 :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind), dimension(5)     :: x_map
          real(rkind), dimension(5)     :: y_map
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          logical                       :: test_loc

          real(rkind), dimension(6,5,4) :: flux_x
          real(rkind), dimension(5,6,4) :: flux_y
          real(rkind), dimension(5,5,4) :: timedev_b
          real(rkind), dimension(5,5,4) :: timedev_w

          real(rkind)                   :: t

          integer(ikind)                :: i,j
          integer                       :: k

          type(pmodel_eq)               :: p_model


          test_validated = .true.

          call initialize_nodes(nodes,dx,dy)

          do i=1, size(x_map,1)
             x_map(i) = (i-1)*dx
          end do

          do i=1, size(y_map,1)
             y_map(i) = (i-1)*dy
          end do

          call p_model%apply_ic(nodes,x_map,y_map)

          call bc_bugged%ini(p_model)
          call bc_working%ini(p_model)

          call bc_bugged%apply_bc_on_timedev(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         timedev_b)

          call bc_bugged%apply_bc_on_timedev(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         timedev_w)

          !compare the time derivatives
          do k=1,ne
             do j=1,2
                do i=1,5
                   test_loc = is_test_validated(
     $                  timedev_b(i,j,k),
     $                  timedev_w(i,j,k),
     $                  detailled)
                   if(.not.test_loc) then
                      print '(''time_dev diff: '',3I2,2F10.5)',
     $                     i,j,k,
     $                     timedev_b(i,j,k),
     $                     timedev_w(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                end do
             end do

             j=3
             do i=1,2
                test_loc = is_test_validated(
     $               timedev_b(i,j,k),
     $               timedev_w(i,j,k),
     $               detailled)
                if(.not.test_loc) then
                   print '(''time_dev diff: '',3I2,2F10.5)',
     $                  i,j,k,
     $                  timedev_b(i,j,k),
     $                  timedev_w(i,j,k)
                end if
                test_validated = test_validated.and.test_loc
             end do
             do i=4,5
                test_loc = is_test_validated(
     $               timedev_b(i,j,k),
     $               timedev_w(i,j,k),
     $               detailled)
                if(.not.test_loc) then
                   print '(''time_dev diff: '',3I2,2F10.5)',
     $                  i,j,k,
     $                  timedev_b(i,j,k),
     $                  timedev_w(i,j,k)
                end if
                test_validated = test_validated.and.test_loc
             end do

             do j=4,5
                do i=1,5
                   test_loc = is_test_validated(
     $                  timedev_b(i,j,k),
     $                  timedev_w(i,j,k),
     $                  detailled)
                   if(.not.test_loc) then
                      print '(''time_dev diff: '',3I2,2F10.5)',
     $                     i,j,k,
     $                     timedev_b(i,j,k),
     $                     timedev_w(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do

        end function test_apply_bc_on_timedev

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
     $          int(var*10000.)-
     $          sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
           
         end function is_test_validated         


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind)                  , intent(out) :: dx
          real(rkind)                  , intent(out) :: dy

          integer(ikind) :: j


          !space steps
          dx = 0.5
          dy = 0.6


          !mass
          nodes(1,1,1) = 0.5
          nodes(2,1,1) = 0.2
          nodes(3,1,1) = 1.2
          nodes(4,1,1) = 5.0
          nodes(5,1,1) = 0.6

          nodes(1,2,1) = 3.0
          nodes(2,2,1) = 4.2
          nodes(3,2,1) = 11.0
          nodes(4,2,1) = 10.6
          nodes(5,2,1) = 5.2

          nodes(1,3,1) = -14.2
          nodes(2,3,1) = 23
          nodes(3,3,1) = 9.8
          nodes(4,3,1) = 3.4
          nodes(5,3,1) = 9.1

          nodes(1,4,1) = 2.45
          nodes(2,4,1) = 0.2
          nodes(3,4,1) = 9.0
          nodes(4,4,1) = 5.4
          nodes(5,4,1) =-2.3

          nodes(1,5,1) = 3.6
          nodes(2,5,1) = 0.1
          nodes(3,5,1) = 6.3
          nodes(4,5,1) = 8.9
          nodes(5,5,1) = -4.23


          !momentum-x
          nodes(1,1,2) = 7.012
          nodes(2,1,2) =-6.323
          nodes(3,1,2) = 3.012
          nodes(4,1,2) = 4.5
          nodes(5,1,2) = 9.6
                    
          nodes(1,2,2) = 4.26
          nodes(2,2,2) = 4.23
          nodes(3,2,2) = 4.5
          nodes(4,2,2) = 7.56
          nodes(5,2,2) = 7.21
                    
          nodes(1,3,2) = 0.23
          nodes(2,3,2) = 7.23
          nodes(3,3,2) = 3.1
          nodes(4,3,2) = 8.9
          nodes(5,3,2) = 9.3
                    
          nodes(1,4,2) = 8.23
          nodes(2,4,2) = -3.1
          nodes(3,4,2) = 6.03
          nodes(4,4,2) = 6.25
          nodes(5,4,2) = 5.12
                    
          nodes(1,5,2) = 3.2
          nodes(2,5,2) = 8.12
          nodes(3,5,2) = 8.9
          nodes(4,5,2) = 4.2
          nodes(5,5,2) = 7.8


          !momentum-y
          nodes(1,1,3) = 7.1
          nodes(2,1,3) = 1.052
          nodes(3,1,3) = 1.23
          nodes(4,1,3) = 7.89
          nodes(5,1,3) = 8.0
                    
          nodes(1,2,3) = 8.362
          nodes(2,2,3) = 4.56
          nodes(3,2,3) = 9.6
          nodes(4,2,3) = 8.96
          nodes(5,2,3) = -3.23
                    
          nodes(1,3,3) = 2.53
          nodes(2,3,3) = -3.23
          nodes(3,3,3) = 7.25
          nodes(4,3,3) = 1.02
          nodes(5,3,3) = 9.26
                    
          nodes(1,4,3) = 8.965
          nodes(2,4,3) = 4.789
          nodes(3,4,3) = 4.56
          nodes(4,4,3) = 3.012
          nodes(5,4,3) = -1.45
                    
          nodes(1,5,3) = 6.26
          nodes(2,5,3) = 5.201
          nodes(3,5,3) = 2.03
          nodes(4,5,3) = 7.89
          nodes(5,5,3) = 9.889


          !total energy
          nodes(1,1,4) = 6.23
          nodes(2,1,4) = 4.12
          nodes(3,1,4) = -3.6
          nodes(4,1,4) = -6.52
          nodes(5,1,4) = 9.57
                    
          nodes(1,2,4) = -0.12
          nodes(2,2,4) = 8.2
          nodes(3,2,4) = 1.2
          nodes(4,2,4) = 7.89
          nodes(5,2,4) = 5.62
                    
          nodes(1,3,4) = -6.23
          nodes(2,3,4) = 6.201
          nodes(3,3,4) = 6.7
          nodes(4,3,4) = 4.12
          nodes(5,3,4) = 1.29
                    
          nodes(1,4,4) = 1.2
          nodes(2,4,4) = 7.958
          nodes(3,4,4) = 1
          nodes(4,4,4) = -5.62
          nodes(5,4,4) = 0.36
                    
          nodes(1,5,4) = 9.6
          nodes(2,5,4) = 6.12
          nodes(3,5,4) = 8.9
          nodes(4,5,4) = 8.95
          nodes(5,5,4) = 6.3

          print '()'
          print '(''mass_density'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,1)
          end do
          print '()'

          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,2)
          end do
          print '()'

          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,3)
          end do
          print '()'

          print '()'
          print '(''total energy'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,4)
          end do
          print '()'

        end subroutine initialize_nodes

      end program test_yoo_new_operators
