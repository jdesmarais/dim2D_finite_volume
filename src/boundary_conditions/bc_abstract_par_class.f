      !> @file
      !> class encapsulating subroutines to compute
      !> the boundary layers in a distributed memory system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the boundary layers in a distributed memory system
      !
      !> @date
      !> 23_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_abstract_par_class

        use field_par_class   , only : field_par
        use cg_operators_class, only : cg_operators

        implicit none


        !> @class bc_abstract_par
        !> encapsulating subroutines to compute boundary
        !> layers in a distributed memory system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine computing the boundary layers in a distributed
        !> memory system
        !---------------------------------------------------------------
        type, abstract :: bc_abstract_par

          contains

          procedure(bc_par), pass, deferred :: apply_bc_on_nodes

        end type bc_abstract_par


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the boundary layers in a distributed
          !> memory system
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p
          !> physical model
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          subroutine bc_par()
          end subroutine bc_par



        end interface




      end module bc_abstract_par_class
