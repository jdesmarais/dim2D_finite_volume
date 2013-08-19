      !> @file
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_class

        use parameters_input, only : nx,ny,ne
        use parameters_kind , only : ikind, rkind
        use surrogate_class , only : surrogate

        implicit none


        private
        public :: field


        !> @class field
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !>
        !> @param nodes
        !> variables computed during the simulation
        !> (ex: mass=nodes(:,:,1),
        !> momentum_x=nodes(:,:,2),
        !> momentum_y=nodes(:,:,3),
        !> energy=nodes(:,:,4))
        !>
        !> @param x_map
        !> discretisation map along the x-axis
        !>
        !> @param y_map
        !> discretisation map along the y-axis
        !>
        !> @param dx
        !> space step along the x-axis
        !>
        !> @param dy
        !> space step along the y-axis
        !---------------------------------------------------------------
        type, extends(surrogate) :: field

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind) :: dx
          real(rkind) :: dy

        end type field

      end module field_class
