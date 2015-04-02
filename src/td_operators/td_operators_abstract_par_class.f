      !> @file
      !> abstract class encapsulating subroutines to compute
      !> the time derivatives of the main variables from the
      !> governing equations using the space discretisation
      !> operators and the physical model in a distributed
      !> memory system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> in the governing equations
      !
      !> @date
      !> 25_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_operators_abstract_par_class

        use bc_operators_par_class, only :
     $       bc_operators_par

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        implicit none


        private
        public :: td_operators_abstract_par


        !> @class td_operators_abstract_par
        !> abstract class encapsulating operators to compute
        !> the time derivatives of the main variables from 
        !> the governing equations in a distributed memory
        !> system
        !>
        !> @param compute_time_dev
        !> compute the time derivatives
        !---------------------------------------------------------------
        type, abstract :: td_operators_abstract_par

          contains
          procedure(time_proc), nopass, deferred :: compute_time_dev

        end type td_operators_abstract_par


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the time derivatives using the
          !> space discretisation operators and the physical model
          !> in a distributed memory system
          !
          !> @date
          !> 25_09_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p_model
          !> physical model
          !
          !>@param bc_par_used
          !> boundary conditions for a distributed memory system
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          function time_proc(
     $       comm_2d, usr_rank,
     $       time, nodes, x_map, y_map,
     $       s,p_model,bc_par_used)
     $       result(time_dev)

            import bc_operators_par
            import sd_operators
            import pmodel_eq
            import rkind
            import nx,ny,ne

            integer                         , intent(in) :: comm_2d
            integer                         , intent(in) :: usr_rank
            real(rkind)                     , intent(in) :: time
            real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
            real(rkind), dimension(nx)      , intent(in) :: x_map
            real(rkind), dimension(ny)      , intent(in) :: y_map
            type(sd_operators)              , intent(in) :: s
            type(pmodel_eq)                 , intent(in) :: p_model
            type(bc_operators_par)          , intent(in) :: bc_par_used
            real(rkind), dimension(nx,ny,ne)             :: time_dev

          end function time_proc

        end interface

      end module td_operators_abstract_par_class
