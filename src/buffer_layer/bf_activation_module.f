      !> @file
      !> module encapsulating the subroutines testing whether
      !> a grid point undermine the application of the open
      !> boundary conditions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines testing whether
      !> a grid point undermine the application of the open
      !> boundary conditions
      !
      !> @date
      ! 26_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_activation_module

        use parameters_input, only : ne
        use parameters_kind , only : rkind

        implicit none

        private
        public :: are_openbc_undermined

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the open boundary conditions are undermined
        !> at the point identified by the nodes
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes_var
        !> governing variables at the grid point tested
        !
        !>@return undermined
        !> logical stating whether the open boundary conditions
        !> are undermined
        !--------------------------------------------------------------
        function are_openbc_undermined(nodes_var) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_var
          logical                                :: undermined

          real(rkind) :: d_liq, d_vap

          d_liq = 1.1-0.05*(1.1-0.1)
          d_vap = 0.1+0.05*(1.1-0.1)

          if((nodes_var(1).ge.d_vap).and.(nodes_var(1).le.d_liq)) then
             undermined = .true.
          else
             undermined = .false.
          end if

          !undermined = .true.

        end function are_openbc_undermined

      end module bf_activation_module
      
