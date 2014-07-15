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
      module bc_operators_abstract_par_class

        use sd_operators_class , only : sd_operators
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind,rkind
        use pmodel_eq_class    , only : pmodel_eq


        implicit none


        !> @class bc_operators_abstract_par
        !> encapsulating subroutines to compute boundary
        !> layers in a distributed memory system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine computing the boundary layers in a distributed
        !> memory system
        !---------------------------------------------------------------
        type,abstract :: bc_operators_abstract_par

          contains

          procedure(ini_par)   , pass, deferred :: ini
          procedure(nodes_par) , pass, deferred :: apply_bc_on_nodes
          procedure(fluxes_par), pass, deferred :: apply_bc_on_fluxes

        end type bc_operators_abstract_par


        abstract interface


          subroutine ini_par(this, comm_2d, p_model)
            
            import bc_operators_abstract_par
            import pmodel_eq

            class(bc_operators_abstract_par), intent(inout) :: this
            integer                         , intent(in)    :: comm_2d
            type(pmodel_eq)                 , intent(in)    :: p_model

          end subroutine ini_par



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
          !>@param comm_2d
          !> integer identifying the general communicator
          !
          !>@param usr_rank
          !> integer identifying the processor in the global communicator
          !
          !>@param nodes
          !> main variables of the governing equations
          !--------------------------------------------------------------
          subroutine nodes_par(this, comm_2d, usr_rank, nodes)

            import bc_operators_abstract_par
            import rkind
            import nx,ny,ne

            class(bc_operators_abstract_par), intent(in)    :: this
            integer                         , intent(in)    :: comm_2d
            integer                         , intent(in)    :: usr_rank
            real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

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
     $       this, comm_2d, usr_rank,
     $       nodes, dx, dy, s_op,
     $       flux_x, flux_y)

            import bc_operators_abstract_par
            import rkind
            import nx,ny,ne
            import sd_operators

            class(bc_operators_abstract_par)  , intent(in)    :: this
            integer                           , intent(in)    :: comm_2d
            integer                           , intent(in)    :: usr_rank
            real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
            real(rkind)                       , intent(in)    :: dx
            real(rkind)                       , intent(in)    :: dy
            type(sd_operators)                , intent(in)    :: s_op
            real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
            real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          end subroutine fluxes_par

        end interface

      end module bc_operators_abstract_par_class
