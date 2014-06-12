      module bf_activation_module

        use parameters_input, only : ne
        use parameters_kind , only : rkind

        implicit none

        private
        public :: is_openbc_undermined

        contains


        !< check if the open boundary conditions are undermined
        !> at the point identified by the nodes
        function is_openbc_undermined(nodes_var) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes_var
          logical                                :: undermined

          real(rkind) :: d_liq, d_vap

          d_liq = 1.1-0.1*(1.1-0.1)
          d_vap = 0.1+0.1*(1.1-0.1)

          if((nodes_var(1).ge.d_vap).and.(nodes_var(1).le.d_liq)) then
             undermined = .true.
          else
             undermined = .false.
          end if

          !undermined = .true.

        end function is_openbc_undermined

      end module bf_activation_module
      
