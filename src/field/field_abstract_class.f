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
      ! 07_08_2013 - initial version               - J.L. Desmarais
      ! 14_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_abstract_class
      
        use bc_operators_class, only : bc_operators
        use parameters_input  , only : nx,ny,ne,bc_size,
     $                                 x_min,x_max,
     $                                 y_min,y_max
        use parameters_kind   , only : ikind, rkind
        use pmodel_eq_class   , only : pmodel_eq
        use sd_operators_class, only : sd_operators
        use surrogate_class   , only : surrogate
        use td_operators_class, only : td_operators

        implicit none


        private
        public :: field_abstract


        !> @class field_abstract
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
        type, extends(surrogate) :: field_abstract

          type(sd_operators) :: sd_operators_used
          type(pmodel_eq)    :: pmodel_eq_used
          type(bc_operators) :: bc_operators_used
          type(td_operators) :: td_operators_used

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind) :: dx
          real(rkind) :: dy

          contains

          procedure, pass :: ini
          procedure, pass :: ini_coordinates
          procedure, pass :: compute_time_dev
          procedure, pass :: apply_bc_on_nodes
          procedure, pass :: compute_integration_step

          procedure, pass :: set_dx    !only for tests
          procedure, pass :: set_dy    !only for tests
          procedure, pass :: get_nodes !only for tests
          procedure, pass :: set_nodes !only for tests

        end type field_abstract


        interface

           subroutine integration_step_proc(
     $          nodes, dt, nodes_tmp, time_dev)
           
             import rkind
             import nx,ny,ne

             real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
             real(rkind)                     , intent(in)    :: dt
             real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
             real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev

           end subroutine integration_step_proc

        end interface


        contains

        !initialize the field_abstract by initializing
        !the boundary conditions bc_operators_used
        subroutine ini(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          call this%bc_operators_used%ini(this%pmodel_eq_used)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to integrate the governing equations using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param x_min
        !> coordinate along the x-axis of the SW border
        !
        !>@param x_max
        !> coordinate along the x-axis of the NE border
        !
        !>@param y_min
        !> coordinate along the y-axis of the SW border
        !
        !>@param y_max
        !> coordinate along the y-axis of the NE border
        !--------------------------------------------------------------
        subroutine ini_coordinates(this)

          implicit none

          class(field_abstract), intent(inout) :: this


          integer(ikind) :: i,j


          !< initialize the space steps along the 
          !> x and y directions
          this%dx = (x_max-x_min)/(nx-1-2*bc_size)
          this%dy = (y_max-y_min)/(ny-1-2*bc_size)


          !< initialize the coordinates along the
          !> x-direction
          do i=1, nx
             this%x_map(i)=x_min + (i-1-bc_size)*this%dx
          end do


          !< initialize the coordinates along the
          !> y-direction
          do j=1, ny
             this%y_map(j)=y_min + (j-1-bc_size)*this%dy
          end do

        end subroutine ini_coordinates


        !compute the time derivative
        function compute_time_dev(this) result(time_dev)

          implicit none

          class(field_abstract), intent(in) :: this
          real(rkind), dimension(nx,ny,ne)  :: time_dev

          !make use of the time discretization operator
          !to compute the time derivative of the field
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%nodes,
     $         this%dx,
     $         this%dy,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used)

        end function compute_time_dev


        !apply the boundary conditions on the grid points
        !of the field
        subroutine apply_bc_on_nodes(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          call this%bc_operators_used%apply_bc_on_nodes(this%nodes)

        end subroutine apply_bc_on_nodes


        !compute the integration step
        subroutine compute_integration_step(
     $     this, dt, nodes_tmp, time_dev, integration_step)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(integration_step_proc)                :: integration_step

          call integration_step(this%nodes, dt, nodes_tmp, time_dev)

        end subroutine compute_integration_step


        !set the attribute dx
        subroutine set_dx(this, dx)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind)          , intent(in)    :: dx

          this%dx = dx

        end subroutine set_dx

      
        !set the attribute dy
        subroutine set_dy(this, dy)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind)          , intent(in)    :: dy

          this%dy = dy

        end subroutine set_dy


        !get the attribute nodes
        subroutine get_nodes(this,nodes)

          implicit none

          class(field_abstract)           , intent(in)  :: this
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes

          nodes = this%nodes

        end subroutine get_nodes


        !set the attribute nodes
        subroutine set_nodes(this,nodes)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes

          this%nodes = nodes

        end subroutine set_nodes

      end module field_abstract_class
