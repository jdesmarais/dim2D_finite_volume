      !> @file
      !> class encapsulating subroutines to compute 
      !> boundary layers in case of periodic xy boundary
      !> conditions on a parallel distributed system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute 
      !> boundary layers in case of periodic xy boundary
      !> conditions on a parallel distributed system
      !
      !> @date
      ! 27_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_par_class

        use bc_abstract_par_class, only : bc_abstract_par
        use cg_operators_class   , only : cg_operators
        use dim2d_eq_class       , only : dim2d_eq
        use field_par_class      , only : field_par
        use mpi_process_class    , only : mpi_process
                                 
        use parameters_constant  , only :
     $     x_direction,
     $     y_direction,
     $     only_compute_proc,
     $     only_exchange_proc
        use parameters_input     , only : nx,ny,ne
        use parameters_kind      , only : ikind,rkind

        use periodic_xy_par_module, only : 
     $     only_compute_along_x, only_compute_along_y,
     $     only_exchange


        implicit none

        private
        public  :: bc_operators_par

        !> @class bc_operators_par
        !> class encapsulating subroutines to compute 
        !> boundary layers in case of periodic xy boundary
        !> conditions on a parallel distributed system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine applying the periodic xy boundary
        !> conditions on the grid points
        !---------------------------------------------------------------
        type, extends(bc_abstract_par) :: bc_operators_par

          contains

          procedure, pass :: apply_bc_on_nodes

        end type bc_operators_par


        contains


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> subroutine to compute the boundary layers using the
          !> periodic boundary conditions in a distributed
          !> memory system
          !
          !> @date
          !> 27_08_2013 - initial version - J.L. Desmarais
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
          subroutine apply_bc_on_nodes(this, f_used, nodes, s_op, p_model)

            implicit none
            
            class(bc_operators_par)         , intent(in)    :: this
            type(field_par)                 , intent(inout) :: f_used
            real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
            type(cg_operators)              , intent(in)    :: s_op
            type(dim2d_eq)                  , intent(in)    :: p_model
            
            
            integer           :: bc_size,neq
            type(mpi_process) :: mpi_op
            
            
            bc_size = s_op%get_bc_size()
            neq     = p_model%get_eq_nb()
            
            
            !compute the boundary layers along the x-direction
            select case(this%proc_x_choice)
            
              case(only_compute_proc)
                 call only_compute_along_x(nodes,bc_size)
            
              case(only_exchange_proc)
                 call only_exchange(this,f_used,nodes,x_direction)

              case default
                 call mpi_op%finalize_mpi()
                 print '(''bc_operators_par_class'')'
                 print '(''periodic_xy_par'')'
                 stop 'proc_x_choice not recognized'
                 
            end select
            
            
            !compute the boundary layers along the y-direction
            select case(this%proc_y_choice)
            
              case(only_compute_proc)
                 call only_compute_along_y(nodes,bc_size)
            
              case(only_exchange_proc)
                 call only_exchange(this,f_used,nodes,y_direction)

              case default
                 call mpi_op%finalize_mpi()
                 print '(''bc_operators_par_class'')'
                 print '(''periodic_xy_par'')'
                 stop 'proc_y_choice not recognized'
                 
            end select

          end subroutine apply_bc_on_nodes

      end module bc_operators_par_class
