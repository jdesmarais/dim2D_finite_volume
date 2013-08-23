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
      module bc_operators_par_class

        use cg_operators_class, only : cg_operators
        use dim2d_eq_class    , only : dim2d_eq
        use field_par_class   , only : field_par
        use parameters_input  , only : nx,ny,ne

        implicit none


        !> @class bc_operators_par
        !> class encapsulating subroutines to compute boundary
        !> layers using reflection boundary conditions in a
        !> distributed memory system
        !>
        !> @param apply_bc_on_nodes
        !> subroutine computing the boundary layers using
        !> reflection boundary conditions in a distributed
        !> memory system
        !---------------------------------------------------------------
        type, extends(mpi_mg_bc_ext) :: bc_operators_par

          contains

          procedure, pass :: apply_bc_on_nodes

        end type bc_operators_par


        contains

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the boundary layers using reflection
          !> boundary conditions in a distributed memory system
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
          !>@param p_model
          !> physical model
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          subroutine apply_bc_on_nodes(this, f_used, nodes, s_op, p_model)

            implicit none

            class(bc_operators_par)         , intent(in)    :: this
            type(field_par)                 , intent(inout) :: f_used
            real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
            type(cg_operators)              , intent(in)    :: s_op
            type(dim2d_eq)                  , intent(in)    :: p_model


            integer :: bc_size
            integer, dimension(ne) :: x_prefactor
            integer, dimension(ne) :: y_prefactor

            
            !< compute the prefactors for the reflection
            x_prefactor = reflection_x_prefactor(p_model)
            y_prefactor = reflection_y_prefactor(p_model)


            !< compute the boundary layers in the x direction
            !> choose the type of procedure for the computation
            select case(this%proc_x_choice)

              case(only_compute_proc)
                 call only_compute_along_x(nodes, p_model)


              case(only_exchange_proc)
                 call only_exchange(this, f_used, nodes, x_direction)
                 



              case(compute_and_exchange_proc)
                 call compute_and_exchange(this)


            end select


            !< compute the boundary layers in the y direction
            !> choose the type of procedure for the computation
            select case()
            end select

          end subroutine apply_bc_on_nodes




      end module bc_operators_par_class
