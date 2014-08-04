      !test the subroutine compute_fluxes_at_the_edges_2nd_order
      program test_openbc_operators

        use openbc_operators_module, only :
     $     compute_fluxes_at_the_edges_2ndorder

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind, rkind

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
        real(rkind), dimension(nx)         :: test_data
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
        
        if((nx.ne.6).or.(ny.ne.5)) then
           print '(''the test requires nx=7 and ny>4'')'
           stop 'change nx and ny'
        end if


        !> nodes initialization
        call initialize_nodes(nodes,dx,dy)
        

        !> test data
        test_data(1) =  74.4d0 !<test dfdx(1)
        test_data(2) =  24d0   !<test dfdx(2)
        test_data(3) = -19.6d0 !<test dfdx(3)
        test_data(4) = -0.7d0  !<test dfdx(4)
        test_data(5) =  3.3d0  !<test dfdx(5)
        test_data(6) = -4.8d0  !<test dfdx(6)


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
        detailled = .true.
        call test_space_derivative(flux_x, dx, 3, 1, test_data, detailled)
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


        subroutine test_space_derivative(
     $       fluxes, dx, j, k, test_data, detailled)
        
          implicit none

          real(rkind), dimension(nx+1,ny,ne), intent(in) :: fluxes
          real(rkind)                       , intent(in) :: dx
          integer(ikind)                    , intent(in) :: j
          integer(ikind)                    , intent(in) :: k
          real(rkind), dimension(nx)        , intent(in) :: test_data
          logical                           , intent(in) :: detailled


          real(rkind), dimension(:), allocatable :: dev
          integer :: i
          logical :: loc

          allocate(dev(size(fluxes,1)-1))

          do i=1, size(dev,1)
             dev(i) = 1.0/dx*(-fluxes(i,j,k)+fluxes(i+1,j,k))
          end do

          if(detailled) then
             do i=1, size(dev,1)
                loc = is_test_validated(dev(i), test_data(i), detailled)
                print '(''test dev('',I2,''): '', L3)', i, loc
             end do
          else
             i=1
             test_validated=is_test_validated(dev(i), test_data(i), detailled)
             do i=2, size(dev,1)
                loc = is_test_validated(dev(i), test_data(i), detailled)
                test_validated = test_validated.and.loc
             end do
             print '(''test validated: '',L3)', loc
          end if

          deallocate(dev)

        end subroutine test_space_derivative


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
