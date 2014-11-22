      program test_lodi_transverse

         use lodi_transverse_module, only :
     $        compute_edge_fluxes,
     $        compute_lodi_terms,
     $        compute_dev_from_flux_x,
     $        compute_dev_from_flux_y

         use ns2d_parameters, only :
     $        viscous_r,
     $        Re,
     $        Pr,
     $        gamma,
     $        mach_infty,
     $        epsilon

         use ns2d_prim_module, only :
     $        speed_of_sound,
     $        cons_lodi_matrix_x,
     $        cons_lodi_matrix_y

         use ns2d_fluxes_module, only :
     $        flux_x_mass_density,
     $        flux_x_momentum_x,
     $        flux_x_momentum_y,
     $        flux_x_total_energy,
     $        flux_y_mass_density,
     $        flux_y_momentum_x,
     $        flux_y_momentum_y,
     $        flux_y_total_energy,
     $        
     $        flux_x_inviscid_momentum_x,
     $        flux_x_inviscid_momentum_y,
     $        flux_x_inviscid_total_energy,
     $        flux_y_inviscid_momentum_x,
     $        flux_y_inviscid_momentum_y,
     $        flux_y_inviscid_total_energy,
     $        
     $        flux_x_viscid_momentum_x,
     $        flux_x_viscid_momentum_y,
     $        flux_x_viscid_total_energy,
     $        flux_y_viscid_momentum_x,
     $        flux_y_viscid_momentum_y,
     $        flux_y_viscid_total_energy

         use parameters_kind, only :
     $        ikind,
     $        rkind

         use sd_operators_class, only :
     $        sd_operators

         implicit none

         real(rkind), dimension(5,5,4) :: nodes
         real(rkind)                   :: dx
         real(rkind)                   :: dy
         type(sd_operators)            :: s

         logical                       :: detailled
         logical                       :: test_loc
         logical                       :: test_validated         

         detailled = .false.
         if(
     $        (.not.is_test_validated(viscous_r,-2.0d0/3.0d0,detailled)).or.
     $        (.not.is_test_validated(Re,10.0d0,detailled)).or.
     $        (.not.is_test_validated(Pr,1.0d0,detailled)).or.
     $        (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $        (.not.is_test_validated(mach_infty,1.0d0,detailled))) then

            print '(''the test requires: '')'
            print '(''viscous_r=-2/3'')'
            print '(''Re=10.0'')'
            print '(''Pr=1.0'')'
            print '(''gamma=5/3'')'
            print '(''mach_infty=1.0'')'
            stop ''

         end if

         
         test_validated = .true.
         detailled = .false.
         

         !initialize the nodes
         call initialize_nodes(nodes,dx,dy)

         !test compute_edge_fluxes
         test_loc = test_compute_edge_fluxes(nodes,s,dx,dy,detailled)
         test_validated = test_loc.and.test_validated

         !test compute_dev_from_flux
         test_loc = test_compute_dev_from_flux(dx,dy,detailled)
         test_validated = test_loc.and.test_validated

         !test compute_lodi_terms
         test_loc = test_compute_lodi_terms(nodes,dx,dy,detailled)
         test_validated = test_loc.and.test_validated


         contains


         function test_compute_edge_fluxes(nodes,s,dx,dy,detailled)
     $        result(test_validated)

           implicit none

           real(rkind), dimension(5,5,4), intent(in) :: nodes
           real(rkind)                  , intent(in) :: dx
           real(rkind)                  , intent(in) :: dy
           type(sd_operators)           , intent(in) :: s
           logical                      , intent(in) :: detailled
           logical                                   :: test_validated

           integer(ikind) :: i,j
           integer(ikind) :: i_min,i_max,i_offset
           integer(ikind) :: j_min,j_max,j_offset
           
           real(rkind), dimension(5,5,4) :: flux_data
           real(rkind), dimension(1,1,4) :: edge_inviscid_flux
           real(rkind), dimension(1,1,4) :: edge_viscid_flux
           real(rkind), dimension(5,5,4) :: flux

           logical :: test_loc

           test_validated = .true.
           

           !compute the fluxes at (3,3) using the dedicated
           !subroutines flux_x_mass_density... from
           !ns2d_fluxes_module and compare the results to the
           !subroutine computing the fluxes as the addition
           !of the inviscid and the viscid parts
           print '(''test compute_edge_fluxes: x'')'
           print '(''---------------------------'')'
           
           i = 3
           j = 3
           
           flux_data(i,j,1) = flux_x_mass_density(nodes,s,i,j)
           flux_data(i,j,2) = flux_x_momentum_x(nodes,s,i,j,dx,dy)
           flux_data(i,j,3) = flux_x_momentum_y(nodes,s,i,j,dx,dy)
           flux_data(i,j,4) = flux_x_total_energy(nodes,s,i,j,dx,dy)
           
           i_min = 1
           i_max = 1
           i_offset = 2
           j_min = 1
           j_max = 1
           j_offset = 2
           
           call compute_edge_fluxes(
     $          nodes,
     $          s,
     $          dx, dy,
     $          i_min, i_max, i_offset,
     $          j_min, j_max, j_offset,
     $          epsilon,
     $          flux_x_mass_density,
     $          flux_x_inviscid_momentum_x,
     $          flux_x_inviscid_momentum_y,
     $          flux_x_inviscid_total_energy,
     $          flux_x_viscid_momentum_x,
     $          flux_x_viscid_momentum_y,
     $          flux_x_viscid_total_energy,
     $          edge_inviscid_flux,
     $          edge_viscid_flux,
     $          flux)
           
           
           test_loc = compare_fluxes(flux_data,flux,detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''test_validated: '',L1)', test_loc
           end if
           print '(''---------------------------'')'
           print '()'
           
           
           !compute the fluxes at (3,3) using the dedicated
           !subroutines flux_x_mass_density... from
           !ns2d_fluxes_module and compare the results to the
           !subroutine computing the fluxes as the addition
           !of the inviscid and the viscid parts
           print '(''test compute_edge_fluxes: y'')'
           print '(''---------------------------'')'
           
           i = 3
           j = 3
           
           flux_data(i,j,1) = flux_y_mass_density(nodes,s,i,j)
           flux_data(i,j,2) = flux_y_momentum_x(nodes,s,i,j,dx,dy)
           flux_data(i,j,3) = flux_y_momentum_y(nodes,s,i,j,dx,dy)
           flux_data(i,j,4) = flux_y_total_energy(nodes,s,i,j,dx,dy)
           
           i_min = 1
           i_max = 1
           i_offset = 2
           j_min = 1
           j_max = 1
           j_offset = 2
           
           call compute_edge_fluxes(
     $          nodes,
     $          s,
     $          dx, dy,
     $          i_min, i_max, i_offset,
     $          j_min, j_max, j_offset,
     $          epsilon,
     $          flux_y_mass_density,
     $          flux_y_inviscid_momentum_x,
     $          flux_y_inviscid_momentum_y,
     $          flux_y_inviscid_total_energy,
     $          flux_y_viscid_momentum_x,
     $          flux_y_viscid_momentum_y,
     $          flux_y_viscid_total_energy,
     $          edge_inviscid_flux,
     $          edge_viscid_flux,
     $          flux)
           
           
           test_loc = compare_fluxes(flux_data,flux,detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''test_validated: '',L1)', test_loc
           end if
           print '(''---------------------------'')'
           print '()'

         end function test_compute_edge_fluxes


         function test_compute_dev_from_flux(dx,dy,detailled)
     $     result(test_validated)

           implicit none

           real(rkind)                  , intent(in) :: dx
           real(rkind)                  , intent(in) :: dy
           logical                      , intent(in) :: detailled
           logical                                   :: test_validated

           real(rkind), dimension(2,1,4) :: edge_flux_x
           real(rkind), dimension(1,2,4) :: edge_flux_y
           real(rkind), dimension(4)     :: test_data_dev_x
           real(rkind), dimension(4)     :: test_data_dev_y
           real(rkind), dimension(4)     :: dev_x
           real(rkind), dimension(4)     :: dev_y
           logical                       :: test_loc

           test_validated = .true.


           print '(''test compute_dev_from_flux'')'
           print '(''------------------------------'')'

           !edge flux x input data
           edge_flux_x(1,1,1) = 2.3
           edge_flux_x(1,1,2) = 8.9
           edge_flux_x(1,1,3) = -1.2
           edge_flux_x(1,1,4) = 6.23
                             
           edge_flux_x(2,1,1) = 9.125
           edge_flux_x(2,1,2) =-6.236
           edge_flux_x(2,1,3) = 7.123
           edge_flux_x(2,1,4) = 9.10

           !edge flux y input data
           edge_flux_y(1,1,1) = 3.126
           edge_flux_y(1,1,2) =-2.156
           edge_flux_y(1,1,3) = 1.023
           edge_flux_y(1,1,4) = 9.145
                     
           edge_flux_y(1,2,1) =-2.364
           edge_flux_y(1,2,2) = 7.122
           edge_flux_y(1,2,3) = 9.183
           edge_flux_y(1,2,4) =-4.156

           !test_data for dev_x and dev_y
           test_data_dev_x(1) = 13.65
           test_data_dev_x(2) =-30.272
           test_data_dev_x(3) = 16.646
           test_data_dev_x(4) =  5.74
           
           test_data_dev_y(1) =- 9.15
           test_data_dev_y(2) = 15.46333333
           test_data_dev_y(3) = 13.6
           test_data_dev_y(4) =-22.16833333
           
           !test compute_dev_from_flux_x
           dev_x = compute_dev_from_flux_x(
     $          edge_flux_x,1,1,dx)

           test_loc = is_vector_validated(
     $          dev_x,
     $          test_data_dev_x,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''x:test_validated: '',L1)', test_loc
           end if

           !test compute_dev_from_flux_y
           dev_y = compute_dev_from_flux_x(
     $          edge_flux_y,1,1,dy)

           test_loc = is_vector_validated(
     $          dev_y,
     $          test_data_dev_y,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''y:test_validated: '',L1)', test_loc
           end if

           print '(''---------------------------'')'
           print '()'

         end function test_compute_dev_from_flux


         function test_compute_lodi_terms(nodes,dx,dy,detailled)
     $     result(test_validated)

           implicit none

           real(rkind), dimension(:,:,:), intent(in) :: nodes
           real(rkind)                  , intent(in) :: dx
           real(rkind)                  , intent(in) :: dy
           logical                      , intent(in) :: detailled
           logical                                   :: test_validated

           integer(ikind)                :: i_offset
           integer(ikind)                :: j_offset
           real(rkind), dimension(2,1,4) :: edge_inviscid_flux_x
           real(rkind), dimension(2,1,4) :: edge_viscid_flux_x
           real(rkind), dimension(1,2,4) :: edge_inviscid_flux_y
           real(rkind), dimension(1,2,4) :: edge_viscid_flux_y
           real(rkind), dimension(1,1,4) :: test_data_transverse_lodi
           real(rkind), dimension(1,1,4) :: test_data_viscous_lodi
           real(rkind), dimension(1,1,4) :: transverse_lodi
           real(rkind), dimension(1,1,4) :: viscous_lodi
           logical                       :: test_loc

           test_validated = .true.


           print '(''test compute_lodi_terms: x '')'
           print '(''---------------------------'')'

           !offset data
           i_offset = 2
           j_offset = 2
           
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
           call compute_lodi_terms(
     $          nodes,dx,i_offset,j_offset,
     $          epsilon,
     $          edge_inviscid_flux_x,
     $          edge_viscid_flux_x,
     $          cons_lodi_matrix_x,
     $          compute_dev_from_flux_x,
     $          transverse_lodi,
     $          viscous_lodi)

           test_loc = is_matrix_validated(
     $          transverse_lodi,
     $          test_data_transverse_lodi,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''I:test_validated: '',L1)', test_loc
           end if

           test_loc = is_matrix_validated(
     $          viscous_lodi,
     $          test_data_viscous_lodi,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''V:test_validated: '',L1)', test_loc
           end if

           print '(''---------------------------'')'
           print '()'


           print '(''test compute_lodi_terms: y '')'
           print '(''---------------------------'')'

           !offset data
           i_offset = 2
           j_offset = 2
           
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
           call compute_lodi_terms(
     $          nodes,dy,i_offset,j_offset,
     $          epsilon,
     $          edge_inviscid_flux_y,
     $          edge_viscid_flux_y,
     $          cons_lodi_matrix_y,
     $          compute_dev_from_flux_y,
     $          transverse_lodi,
     $          viscous_lodi)

           test_loc = is_matrix_validated(
     $          transverse_lodi,
     $          test_data_transverse_lodi,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''I:test_validated: '',L1)', test_loc
           end if

           test_loc = is_matrix_validated(
     $          viscous_lodi,
     $          test_data_viscous_lodi,
     $          detailled)
           test_validated = test_validated.and.test_loc
           if(.not.detailled) then
              print '(''V:test_validated: '',L1)', test_loc
           end if

           print '(''---------------------------'')'
           print '()'

         end function test_compute_lodi_terms


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


         function compare_fluxes(flux_data,flux,detailled)
     $     result(test_validated)

           implicit none

           real(rkind), dimension(:,:,:), intent(in) :: flux_data
           real(rkind), dimension(:,:,:), intent(in) :: flux
           logical                      , intent(in) :: detailled
           logical                                   :: test_validated

           integer :: k
           logical :: test_loc

           test_validated = .true.
           do k=1,4
              
              test_loc = is_test_validated(
     $             flux_data(3,3,k),
     $             flux(3,3,k),detailled)

              test_validated = test_validated.and.test_loc

           end do

         end function compare_fluxes


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

      end program test_lodi_transverse
