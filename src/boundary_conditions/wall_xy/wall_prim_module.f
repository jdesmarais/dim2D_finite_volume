      !> @file
      !> module encapsulating subroutines computing special
      !> variables out of conservative variables for the
      !> Diffuse Interface Model governing equations when
      !> applying wall boundary conditions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing special
      !> variables out of conservative variables for the
      !> Diffuse Interface Model governing equations when
      !> applying wall boundary conditions
      !
      !> @date
      !> 24_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_prim_module
      
        use dim2d_parameters, only : cv_r
        use field_class     , only : field
        use parameters_kind , only : ikind, rkind

        implicit none

        private
        public :: wall_pressure

        contains

        
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the pressure at the wall
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
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
        !> \f$ P_{\textrm{wall}} \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function wall_pressure(field_used,i,j) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          real(rkind)                :: var


          if(rkind.eq.8) then
             var=3.0d0/((3.0d0-field_used%nodes(i,j,1))*cv_r)*(
     $            field_used%nodes(i,j,4)
     $            + 3.0d0*field_used%nodes(i,j,1)**2)
     $            - 3.0d0*field_used%nodes(i,j,1)**2
          else
             var=3./((3.-field_used%nodes(i,j,1))*cv_r)*(
     $            field_used%nodes(i,j,4)
     $            + 3*field_used%nodes(i,j,1)**2)
     $            - 3*field_used%nodes(i,j,1)**2
          end if

        end function wall_pressure

      end module wall_prim_module

        
