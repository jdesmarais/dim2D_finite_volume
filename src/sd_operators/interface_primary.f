      !> @file
      !> module implementing the abstract interface for procedures
      !> computing primary variables
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> abstract interface for procedures computing primary variables
      !> from conservative variables
      !
      !> @date
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module interface_primary

        use field_class    , only : field
        use parameters_kind, only : rkind, ikind

        implicit none

        private
        public :: get_primary_var


        abstract interface


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for procedure computing primary variables at [i,j]
          !> (ex: pressure, temperature)
          !
          !> @date
          !> 07_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the data
          !
          !>@param i
          !> index along x-axis where the data is evaluated
          !
          !>@param j
          !> index along y-axis where the data is evaluated
          !
          !>@param var
          !> primary variable evaluated at [i,j]
          !---------------------------------------------------------------
          function get_primary_var(field_used,i,j) result(var)

            import field
            import ikind
            import rkind

            class(field)       , intent(in) :: field_used
            integer(ikind)     , intent(in) :: i
            integer(ikind)     , intent(in) :: j
            real(rkind)                     :: var

          end function get_primary_var

        end interface

      end module interface_primary
