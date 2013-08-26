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

        use cg_operators_class , only : cg_operators
        use dim2d_eq_class     , only : dim2d_eq
        use field_par_class    , only : field_par
        use mpi_mg_bc_ext_class, only : mpi_mg_bc_ext
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind,rkind

        implicit none


        !> @class bc_abstract_par
        !> encapsulating subroutines to compute boundary
        !> layers in a distributed memory system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine computing the boundary layers in a distributed
        !> memory system
        !---------------------------------------------------------------
        type, abstract, extends(mpi_mg_bc_ext) :: bc_abstract_par

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
          !>@param f_used
          !> object encapsulating the main variables
          !
          !>@param nodes
          !> main variables of the governing equations
          !
          !>@param s_op
          !> space discretization operators
          !
          !>@param p_model
          !> physical model
          !--------------------------------------------------------------
          subroutine bc_par(this, f_used, nodes, s_op, p_model)

            import bc_abstract_par
            import field_par
            import rkind
            import nx,ny,ne
            import cg_operators
            import dim2d_eq

            class(bc_abstract_par)          , intent(in)    :: this
            type(field_par)                 , intent(inout) :: f_used
            real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
            type(cg_operators)              , intent(in)    :: s_op
            type(dim2d_eq)                  , intent(in)    :: p_model

          end subroutine bc_par

        end interface

      end module bc_abstract_par_class
