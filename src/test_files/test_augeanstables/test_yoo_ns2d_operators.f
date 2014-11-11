      program test_yoo_ns2d_operators

        use bc_operators_class, only :
     $       bc_operators

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

        type(bc_operators) :: bc_used
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
        test_loc = test_compute_edge_fluxes(bc_used, detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_edge_fluxes: '',L1)', test_loc
        

        detailled      = .false.

        !test of the LODI terms computation
        test_loc = test_compute_lodi_terms_edge(bc_used, detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_lodi_terms_edge: '',L1)', test_loc


        contains

        
        !test the function bc_used%compute_edge_fluxes
        function test_compute_edge_fluxes(bc_used, detailled)
     $       result(test_validated)

          implicit none

          class(bc_operators), intent(in) :: bc_used
          logical            , intent(in) :: detailled
          logical                         :: test_validated

          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy
          logical                       :: test_loc

          type(pmodel_eq)               :: p_model

          real(rkind), dimension(6,5,4) :: flux_x_data
          real(rkind), dimension(5,6,4) :: flux_y_data
          real(rkind), dimension(6,5,4) :: flux_x
          real(rkind), dimension(5,6,4) :: flux_y

          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0

          real(rkind), dimension(2,2,4) :: transverse_lodi_N
          real(rkind), dimension(2,2,4) :: transverse_lodi_S
          real(rkind), dimension(2,2,4) :: transverse_lodi_E
          real(rkind), dimension(2,2,4) :: transverse_lodi_W

          real(rkind), dimension(2,2,4) :: viscous_lodi_N
          real(rkind), dimension(2,2,4) :: viscous_lodi_S
          real(rkind), dimension(2,2,4) :: viscous_lodi_E
          real(rkind), dimension(2,2,4) :: viscous_lodi_W

          integer(ikind) :: i,j
          integer        :: k

          test_validated = .true.

          
          !initialize the nodes
          call initialize_nodes(nodes,dx,dy)

          !test the computation of the edge fluxes
          !the computation of the fluxes should give the same results
          !as the computation using the fluxes_oneside

          !create the flux data
          !y-fluxes
          do j=3,4
             i=1
             flux_y_data(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes, dx, dy, i,j, s_x_L0)

             i=2
             flux_y_data(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes, dx, dy, i,j, s_x_L1)

             i=4
             flux_y_data(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes, dx, dy, i,j, s_x_R1)

             i=5
             flux_y_data(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes, dx, dy, i,j, s_x_R0)

          end do

          !x-fluxes
          j=1
          do i=3,4
             flux_x_data(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes, dx, dy, i,j, s_y_L0)
          end do
          
          j=2
          do i=3,4
             flux_x_data(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes, dx, dy, i,j, s_y_L1)
          end do

          j=4
          do i=3,4
             flux_x_data(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes, dx, dy, i,j, s_y_R1)
          end do

          j=5
          do i=3,4
             flux_x_data(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes, dx, dy, i,j, s_y_R0)
          end do

          !compute the edge fluxes
          call bc_used%compute_edge_fluxes(
     $         p_model,
     $         nodes,dx,dy,
     $         transverse_lodi_N, transverse_lodi_S,
     $         transverse_lodi_E, transverse_lodi_W,
     $         viscous_lodi_N, viscous_lodi_S,
     $         viscous_lodi_E, viscous_lodi_W,
     $         flux_x, flux_y)


          !compare the fluxes
          !in the y-direction
          do k=1,ne
             do j=3,4
                do i=1,2
                   test_loc = is_test_validated(
     $                  flux_y(i,j,k),
     $                  flux_y_data(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do

                do i=4,5
                   test_loc = is_test_validated(
     $                  flux_y(i,j,k),
     $                  flux_y_data(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do

          !in the x-direction
          do k=1,ne
             do j=1,2
                do i=3,4
                   test_loc = is_test_validated(
     $                  flux_x(i,j,k),
     $                  flux_x_data(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do

             do j=4,5
                do i=3,4
                   test_loc = is_test_validated(
     $                  flux_x(i,j,k),
     $                  flux_x_data(i,j,k),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                end do
             end do
          end do

        end function test_compute_edge_fluxes


        !test the functions:
        ! - bc_used%compute_lodi_terms_x_edge()
        ! - bc_used%compute_lodi_terms_y_edge()
        function test_compute_lodi_terms_edge(bc_used, detailled)
     $     result(test_validated)

          implicit none
          
          class(bc_operators), intent(in) :: bc_used
          logical            , intent(in) :: detailled
          logical                         :: test_validated


          real(rkind), dimension(5,5,4) :: nodes
          real(rkind)                   :: dx
          real(rkind)                   :: dy

          type(pmodel_eq)               :: p_model

          integer(ikind)                :: i_offset
          integer(ikind)                :: j_offset

          real(rkind), dimension(2,2,4) :: edge_inviscid_flux_x
          real(rkind), dimension(2,2,4) :: edge_viscid_flux_x

          real(rkind), dimension(2,2,4) :: edge_inviscid_flux_y
          real(rkind), dimension(2,2,4) :: edge_viscid_flux_y

          real(rkind), dimension(1,1,4) :: test_data_transverse_lodi
          real(rkind), dimension(1,1,4) :: test_data_viscous_lodi

          real(rkind), dimension(2,2,4) :: transverse_lodi
          real(rkind), dimension(2,2,4) :: viscous_lodi
          logical                       :: test_loc

          integer :: k

          test_validated = .true.


          !initialize the nodes
          call initialize_nodes(nodes,dx,dy)

          print '(''test compute_lodi_terms: x '')'
          print '(''---------------------------'')'

          !offset data
          i_offset = 3
          j_offset = 3
          
          !inviscid flux data
          edge_inviscid_flux_x(1,1,1) = 2.3
          edge_inviscid_flux_x(1,1,2) = 8.9
          edge_inviscid_flux_x(1,1,3) =-1.2
          edge_inviscid_flux_x(1,1,4) = 6.23
          
          edge_inviscid_flux_x(2,1,1) = 9.125
          edge_inviscid_flux_x(2,1,2) =-6.236
          edge_inviscid_flux_x(2,1,3) = 7.123
          edge_inviscid_flux_x(2,1,4) = 9.10

          !viscid flux data
          edge_viscid_flux_x(1,1,1) = 3.126
          edge_viscid_flux_x(1,1,2) =-2.156
          edge_viscid_flux_x(1,1,3) = 1.023
          edge_viscid_flux_x(1,1,4) = 9.145
          
          edge_viscid_flux_x(2,1,1) =-2.364
          edge_viscid_flux_x(2,1,2) = 7.122
          edge_viscid_flux_x(2,1,3) = 9.183
          edge_viscid_flux_x(2,1,4) =-4.156

          !test_data for lodi terms
          test_data_transverse_lodi(1,1,1) =- 0.668141399 
          test_data_transverse_lodi(1,1,2) =- 0.513611273
          test_data_transverse_lodi(1,1,3) =-26.82262694
          test_data_transverse_lodi(1,1,4) = 16.93005478
          
          test_data_viscous_lodi(1,1,1) =  0.249417951
          test_data_viscous_lodi(1,1,2) = 2.767423823
          test_data_viscous_lodi(1,1,3) =-4.599855541
          test_data_viscous_lodi(1,1,4) =-1.8133795592
          
          !test compute_lodi_terms in x-direction
          call bc_used%compute_lodi_terms_y_edge(
     $         p_model,
     $         nodes,dx,
     $         i_offset,j_offset,
     $         edge_inviscid_flux_x,
     $         edge_viscid_flux_x,
     $         transverse_lodi,
     $         viscous_lodi)

          do k=1, ne
             test_loc = is_test_validated(
     $            transverse_lodi(1,1,k),
     $            test_data_transverse_lodi(1,1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do           

          if(.not.detailled) then
             print '(''I:test_validated: '',L1)', test_loc
          end if

          test_validated =.true.

          do k=1, ne
             test_loc = is_test_validated(
     $            viscous_lodi(1,1,k),
     $            test_data_viscous_lodi(1,1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do  

          if(.not.detailled) then
             print '(''V:test_validated: '',L1)', test_loc
          end if

          print '(''---------------------------'')'
          print '()'


          print '(''test compute_lodi_terms: y '')'
          print '(''---------------------------'')'

          !offset data
          i_offset = 3
          j_offset = 3
          
          !inviscid flux data
          edge_inviscid_flux_y(1,1,1) = 2.3
          edge_inviscid_flux_y(1,1,2) = 8.9
          edge_inviscid_flux_y(1,1,3) =-1.2
          edge_inviscid_flux_y(1,1,4) = 6.23
          
          edge_inviscid_flux_y(1,2,1) = 9.125
          edge_inviscid_flux_y(1,2,2) =-6.236
          edge_inviscid_flux_y(1,2,3) = 7.123
          edge_inviscid_flux_y(1,2,4) = 9.10

          !viscid flux data
          edge_viscid_flux_y(1,1,1) = 3.126
          edge_viscid_flux_y(1,1,2) =-2.156
          edge_viscid_flux_y(1,1,3) = 1.023
          edge_viscid_flux_y(1,1,4) = 9.145
          
          edge_viscid_flux_y(1,2,1) =-2.364
          edge_viscid_flux_y(1,2,2) = 7.122
          edge_viscid_flux_y(1,2,3) = 9.183
          edge_viscid_flux_y(1,2,4) =-4.156

          !test_data for lodi terms
          test_data_transverse_lodi(1,1,1) =  2.941314383
          test_data_transverse_lodi(1,1,2) = -0.428009394
          test_data_transverse_lodi(1,1,3) = -0.670951424
          test_data_transverse_lodi(1,1,4) = -7.572858707
          
          test_data_viscous_lodi(1,1,1) =  0.187323685
          test_data_viscous_lodi(1,1,2) =  2.30618652
          test_data_viscous_lodi(1,1,3) = -3.960424592
          test_data_viscous_lodi(1,1,4) = -1.383934685
          
          !test compute_lodi_terms in x-direction
          call bc_used%compute_lodi_terms_x_edge(
     $         p_model,
     $         nodes,dy,
     $         i_offset,j_offset,
     $         edge_inviscid_flux_y,
     $         edge_viscid_flux_y,
     $         transverse_lodi,
     $         viscous_lodi)

          do k=1, ne
             test_loc = is_test_validated(
     $            transverse_lodi(1,1,k),
     $            test_data_transverse_lodi(1,1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do           

          if(.not.detailled) then
             print '(''I:test_validated: '',L1)', test_loc
          end if

          test_validated =.true.

          do k=1, ne
             test_loc = is_test_validated(
     $            viscous_lodi(1,1,k),
     $            test_data_viscous_lodi(1,1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do  

          if(.not.detailled) then
             print '(''V:test_validated: '',L1)', test_loc
          end if

          print '(''---------------------------'')'
          print '()'
          
        end function test_compute_lodi_terms_edge


        function is_vector_validated(var,cst,detailled) result(test_validated)

           implicit none

           real(rkind), dimension(:), intent(in) :: var
           real(rkind), dimension(:), intent(in) :: cst
           logical                  , intent(in) :: detailled
           logical                               :: test_validated

           logical :: test_loc
           integer :: k

           test_validated = .true.
           
           do k=1, size(var,1)
              test_loc = is_test_validated(
     $             var(k),
     $             cst(k),
     $             detailled)
              test_validated = test_validated.and.test_loc
              if(detailled) then
                 print '(''var('',I2,''):'',L2)', k,test_loc
              end if
           end do
           

         end function is_vector_validated


         function is_matrix_validated(var,cst,detailled) result(test_validated)

           implicit none

           real(rkind), dimension(:,:,:), intent(in) :: var
           real(rkind), dimension(:,:,:), intent(in) :: cst
           logical                      , intent(in) :: detailled
           logical                                   :: test_validated

           logical :: test_loc
           integer :: i,j,k

           test_validated = .true.
           
           do k=1, size(var,3)
              do j=1, size(var,2)
                 do i=1, size(var,1)
                    test_loc = is_test_validated(
     $                   var(i,j,k),
     $                   cst(i,j,k),
     $                   detailled)
                    test_validated = test_validated.and.test_loc
                    if(detailled) then
                       print '(''var('',3I2,''):'',L2)', i,j,k,test_loc
                    end if
                 end do
              end do
           end do
           

         end function is_matrix_validated

      
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


      end program test_yoo_ns2d_operators
