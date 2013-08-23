      !> @file
      !> class encapsulating subroutines to compute 
      !> boundary layers in case of reflection xy boundary
      !> conditions on a parallel distributed system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute 
      !> boundary layers in case of reflection xy boundary
      !> conditions on a parallel distributed system
      !
      !> @date
      ! 22_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_par_class

        implicit none

        private
        public  :: bc_operators_par

        !> @class bc_operators_par
        !> class encapsulating subroutines to compute 
        !> boundary layers in case of reflection xy boundary
        !> conditions on a parallel distributed system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine applying the reflection xy boundary
        !> conditions on the grid points
        !---------------------------------------------------------------
        type, extends(bc_abstract_par) :: bc_operators_par

          contains

          procedure, pass :: apply_bc_on_nodes


        end type bc_operators_par



      end module bc_operators_par_class
