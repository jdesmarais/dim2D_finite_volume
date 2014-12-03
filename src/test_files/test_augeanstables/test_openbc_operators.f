      !test the subroutine compute_fluxes_at_the_edges_2nd_order
      program test_openbc_operators

        use openbc_operators_module, only :
     $     compute_fluxes_at_the_edges_2ndorder

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind

        use sd_operators_class, only :
     $       sd_operators

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

        real(rkind), dimension(nx,ny,ne)   :: nodes
        real(rkind), dimension(nx+1,ny)    :: test_data_x_flux
        real(rkind), dimension(nx,ny+1)    :: test_data_y_flux
        real(rkind)                        :: dx,dy
                                           
        type(pmodel_eq)                    :: p_model
                                           
        type(sd_operators)                 :: s_interior
                                           
        type(sd_operators_x_oneside_L0)    :: s_x_L0
        type(sd_operators_x_oneside_L1)    :: s_x_L1
        type(sd_operators_x_oneside_R1)    :: s_x_R1
        type(sd_operators_x_oneside_R0)    :: s_x_R0
                                           
        type(sd_operators_y_oneside_L0)    :: s_y_L0
        type(sd_operators_y_oneside_L1)    :: s_y_L1
        type(sd_operators_y_oneside_R1)    :: s_y_R1
        type(sd_operators_y_oneside_R0)    :: s_y_R0

        real(rkind), dimension(nx+1,ny,ne) :: flux_x
        real(rkind), dimension(nx,ny+1,ne) :: flux_y

        logical :: detailled
        logical :: test_validated
        logical :: loc

        integer :: i,j



        print '(''************************************'')'
        print '(''WARNING'')'
        print '(''************************************'')'
        print '(''for this test, the simpletest'')'
        print '(''model used must compute the'')'
        print '(''x-fluxes using simply sd_operator%f'')'
        print '(''such that the fluxes table'')'
        print '(''corresponds at the end to the'')'
        print '(''x-gradient'')'
        print '(''************************************'')'
        
        if((nx.ne.7).or.(ny.ne.5)) then
           print '(''the test requires nx=7 and ny=5'')'
           stop 'change nx and ny'
        end if


        !> nodes initialization
        call initialize_nodes(nodes,dx,dy)
        

        !> test data
        test_data_x_flux(3,1) =  0.7d0
        test_data_x_flux(4,1) =  3.1d0 
        test_data_x_flux(5,1) =  2.8d0
        test_data_x_flux(6,1) = -1.5d0
        test_data_x_flux(3,2) =  7.6d0
        test_data_x_flux(4,2) = 10.8d0 
        test_data_x_flux(5,2) =  7.9d0
        test_data_x_flux(6,2) =  3.2d0
        test_data_x_flux(3,3) = 16.4d0
        test_data_x_flux(4,3) =  6.6d0 
        test_data_x_flux(5,3) =  6.25d0
        test_data_x_flux(6,3) =  7.9d0
        test_data_x_flux(3,4) =  4.6d0
        test_data_x_flux(4,4) =  7.2d0 
        test_data_x_flux(5,4) =  1.55d0
        test_data_x_flux(6,4) = -0.65d0
        test_data_x_flux(3,5) =  3.2d0
        test_data_x_flux(4,5) =  7.6d0 
        test_data_x_flux(5,5) =  2.335d0
        test_data_x_flux(6,5) =  2.335d0

        test_data_y_flux(1,3) = -5.6d0
        test_data_y_flux(2,3) = 13.6d0 
        test_data_y_flux(3,3) = 10.4d0
        test_data_y_flux(4,3) =  7d0 
        test_data_y_flux(5,3) =  7.15d0
        test_data_y_flux(6,3) =  3.95d0
        test_data_y_flux(7,3) =  6.005d0
        test_data_y_flux(1,4) = -5.875d0
        test_data_y_flux(2,4) = 11.6d0 
        test_data_y_flux(3,4) =  9.4d0
        test_data_y_flux(4,4) =  4.4d0 
        test_data_y_flux(5,4) =  3.4d0
        test_data_y_flux(6,4) =  3.85d0
        test_data_y_flux(7,4) = -0.75d0


        !> compute the fluxes
        flux_x = p_model%compute_flux_x(nodes,dx,dy,s_interior)
        flux_y = p_model%compute_flux_y(nodes,dx,dy,s_interior)

        call compute_fluxes_at_the_edges_2ndorder(
     $       nodes, dx, dy,
     $       s_x_L0, s_x_L1, s_x_R1, s_x_R0,
     $       s_y_L0, s_y_L1, s_y_R1, s_y_R0,
     $       p_model,
     $       flux_x, flux_y)

          
        !> derivative computation + comparison with test_data
        detailled = .false.


        if(detailled) then
           do j=1,5
              do i=3,6
                 loc = is_test_validated(flux_x(i,j,1), test_data_x_flux(i,j), detailled)
                 print '(''test flux_x('',I2,'','',I2,''): '', L3)', i,j, loc
              end do
           end do

           do j=3,4
              do i=1,7
                 loc = is_test_validated(flux_y(i,j,1), test_data_y_flux(i,j), detailled)
                 print '(''test flux_y('',I2,'','',I2,''): '', L3)', i,j, loc
              end do
           end do

        else
           test_validated=.true.
           do j=1,5
              do i=3,6
                 loc = is_test_validated(flux_x(i,j,1), test_data_x_flux(i,j), detailled)
                 test_validated=test_validated.and.loc
              end do
           end do

           do j=3,4
              do i=1,7
                 loc = is_test_validated(flux_y(i,j,1), test_data_y_flux(i,j), detailled)
                 test_validated=test_validated.and.loc
              end do
           end do
           
           print '(''test_compute_fluxes_at_the_edge_2ndorder: '',L3)', test_validated
        end if
          

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


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          
          dx=0.5
          dy=0.6

          nodes(1,1,1)=0.5
          nodes(2,1,1)=0.2
          nodes(3,1,1)=1.2
          nodes(4,1,1)=5.0
          nodes(5,1,1)=0.6
          nodes(6,1,1)=-3.6
          nodes(7,1,1)=-6.52

          nodes(1,2,1)=3.0
          nodes(2,2,1)=4.2
          nodes(3,2,1)=11.0
          nodes(4,2,1)=10.6
          nodes(5,2,1)=5.2
          nodes(6,2,1)=1.2
          nodes(7,2,1)=7.89

          nodes(1,3,1)=-14.2
          nodes(2,3,1)=23
          nodes(3,3,1)=9.8
          nodes(4,3,1)=3.4
          nodes(5,3,1)=9.1
          nodes(6,3,1)=6.7
          nodes(7,3,1)=4.12

          nodes(1,4,1)=2.45
          nodes(2,4,1)=0.2
          nodes(3,4,1)=9.0
          nodes(4,4,1)=5.4
          nodes(5,4,1)=-2.3
          nodes(6,4,1)=1.0
          nodes(7,4,1)=-5.62

          nodes(1,5,1)=3.6
          nodes(2,5,1)=0.1
          nodes(3,5,1)=6.3
          nodes(4,5,1)=8.9
          nodes(5,5,1)=-4.23
          nodes(6,5,1)=8.9
          nodes(7,5,1)=8.95

        end subroutine initialize_nodes

      end program test_openbc_operators
