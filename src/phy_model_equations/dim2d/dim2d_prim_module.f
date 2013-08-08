      !> @file
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables from the Diffuse Interface
      !> Model governing equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables from the Diffuse Interface
      !> Model governing equations
      !
      !> @date
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_prim_module

        use field_class    , only : field
        use parameters_kind, only : ikind, rkind

        implicit none

        private
        public :: mass_density


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the mass density \f$\rho\f$
        !
        !> @date
        !> 08_08_2012 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function mass_density(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var

          var=field_used%nodes(i,j,1)

        end function mass_density


      end module dim2d_prim_module
