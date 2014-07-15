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
      ! 26_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_par_class

        use bc_operators_abstract_par_class, only:
     $     bc_operators_abstract_par

        use sd_operators_class, only :
     $       sd_operators

        use mpi_mg_bc_ext_class, only :
     $       mpi_mg_bc_ext
                                 
        use parameters_constant, only :
     $       x_direction,
     $       y_direction,
     $       only_compute_proc,
     $       compute_and_exchange_proc,
     $       only_exchange_proc

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use reflection_xy_module, only :
     $       reflection_x_prefactor,
     $       reflection_y_prefactor

        use reflection_xy_par_module, only : 
     $       only_compute_along_x,
     $       only_compute_along_y,
     $       only_exchange,
     $       compute_and_exchange_along_x,
     $       compute_and_exchange_along_y


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
        !>
        !> @param apply_bc_on_fluxes
        !> subroutine applying the reflection xy boundary
        !> conditions on the fluxes
        !---------------------------------------------------------------
        type, extends(bc_operators_abstract_par) :: bc_operators_par

          type(mpi_mg_bc_ext)    :: mpi_messenger_used      

          integer, dimension(ne) :: prefactor_x
          integer, dimension(ne) :: prefactor_y

          contains

          procedure, pass :: ini
          procedure, pass :: apply_bc_on_nodes
          procedure, pass :: apply_bc_on_fluxes

        end type bc_operators_par


        contains

        subroutine ini(this, comm_2d, p_model)
        
          implicit none
          
          class(bc_operators_par), intent(inout) :: this
          integer                , intent(in)    :: comm_2d
          type(pmodel_eq)        , intent(in)    :: p_model

          call this%mpi_messenger_used%ini(comm_2d)
          
          this%prefactor_x = reflection_x_prefactor(p_model)
          this%prefactor_y = reflection_y_prefactor(p_model)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the boundary layers using the
        !> reflection boundary conditions in a distributed
        !> memory system
        !
        !> @date
        !> 26_08_2013 - initial version - J.L. Desmarais
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
        subroutine apply_bc_on_nodes(
     $     this, comm_2d, usr_rank, nodes)

          implicit none
          
          class(bc_operators_par)         , intent(in)    :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          
          
          !compute the boundary layers along the x-direction
          select case(this%mpi_messenger_used%proc_x_choice)
          
            case(only_compute_proc)
               call only_compute_along_x(
     $              nodes,this%prefactor_x)
          
            case(compute_and_exchange_proc)
               call compute_and_exchange_along_x(
     $              this%mpi_messenger_used,
     $              comm_2d,usr_rank,nodes,this%prefactor_x,
     $              this%mpi_messenger_used%exchange_id(x_direction))
          
            case(only_exchange_proc)
               call only_exchange(
     $              this%mpi_messenger_used,
     $              comm_2d,usr_rank,nodes,x_direction)
               
          end select
          
          
          !compute the boundary layers along the y-direction
          select case(this%mpi_messenger_used%proc_y_choice)
          
            case(only_compute_proc)
               call only_compute_along_y(
     $              nodes,this%prefactor_y)
          
            case(compute_and_exchange_proc)
               call compute_and_exchange_along_y(
     $              this%mpi_messenger_used,
     $              comm_2d,usr_rank,nodes,this%prefactor_y,
     $              this%mpi_messenger_used%exchange_id(y_direction))
          
            case(only_exchange_proc)
               call only_exchange(
     $              this%mpi_messenger_used,
     $              comm_2d,usr_rank,nodes,y_direction)
               
          end select

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the reflection boundary layers
        !> of the fluxes in a distributed memory system
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
        subroutine apply_bc_on_fluxes(
     $     this, comm_2d, usr_rank,
     $     nodes, dx, dy, s_op,
     $     flux_x, flux_y)

          implicit none

          class(bc_operators_par)           , intent(in)    :: this
          integer                           , intent(in)    :: comm_2d
          integer                           , intent(in)    :: usr_rank  
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s_op
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y


          !these variables are only used to make sure that there
          !are no warning about unused subroutine arguments
          !as the subroutine is not meant to be used for reflection
          !boundary conditions
          integer :: this_s, comm_2d_s, usr_rank_s, bc_size_s
          real(rkind) :: node_s, dx_s, dy_s, flux
          
          stop 'reflection_xy_par: apply_bc_on_fluxes not implemented'

          this_s=this%mpi_messenger_used%proc_x_choice
          comm_2d_s=comm_2d
          usr_rank_s=usr_rank
          bc_size_s=s_op%get_bc_size()
          node_s=nodes(1,1,1)
          dx_s=dx
          dy_s=dy
          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)

        end subroutine apply_bc_on_fluxes

      end module bc_operators_par_class
