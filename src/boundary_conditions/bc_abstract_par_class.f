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

          procedure(nodes_par) , pass, deferred :: apply_bc_on_nodes
          procedure(fluxes_par), pass, deferred :: apply_bc_on_fluxes

        end type bc_abstract_par


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the boundary layers of the gridpoints
          !> in a distributed memory system
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
          subroutine nodes_par(this, f_used, nodes, s_op, p_model)

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

          end subroutine nodes_par


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the boundary layers of the fluxes
          !> in a distributed memory system
          !
          !> @date
          !> 25_09_2013 - initial version - J.L. Desmarais
          !
          !>@param f_used
          !> object encapsulating the main variables
          !
          !>@param s_op
          !> space discretization operators
          !
          !>@param p_model
          !> physical model
          !
          !>@param flux_x
          !> fluxes along the x-direction
          !
          !>@param flux_y
          !> fluxes along the y-direction
          !--------------------------------------------------------------
          subroutine fluxes_par(
     $       this, f_used, s_op,
     $       flux_x, flux_y)

            import bc_abstract_par
            import field_par
            import rkind
            import nx,ny,ne
            import cg_operators
            import dim2d_eq

            class(bc_abstract_par)            , intent(in)    :: this
            type(field_par)                   , intent(in)    :: f_used
            type(cg_operators)                , intent(in)    :: s_op
            real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
            real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          end subroutine fluxes_par

        end interface

      end module bc_abstract_par_class
