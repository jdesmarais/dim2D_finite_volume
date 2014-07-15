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
        use io_operators_class, only : io_operators
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

          type(sd_operators), private :: sd_operators_used
          type(pmodel_eq)   , private :: pmodel_eq_used
          type(bc_operators), private :: bc_operators_used
          type(td_operators), private :: td_operators_used
          type(io_operators), private :: io_operators_used

          real(rkind), dimension(nx,ny,ne), private :: nodes
          real(rkind), dimension(nx)      , private :: x_map
          real(rkind), dimension(ny)      , private :: y_map
          real(rkind)                     , private :: dx
          real(rkind)                     , private :: dy

          contains

          procedure, pass          :: ini
          procedure, pass, private :: check_inputs
          procedure, pass, private :: ini_coordinates
          procedure, pass, private :: apply_initial_conditions
          procedure, pass          :: compute_time_dev
          procedure, pass          :: apply_bc_on_nodes
          procedure, pass          :: compute_integration_step
          procedure, pass          :: write_data

          procedure, pass          :: set_dx    !only for tests
          procedure, pass          :: set_dy    !only for tests
          procedure, pass          :: get_nodes !only for tests
          procedure, pass          :: set_nodes !only for tests
          procedure, pass          :: get_x_map !only for tests
          procedure, pass          :: get_y_map !only for tests

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

        !initialize the field_abstract by:
        ! - initializing the boundary conditions bc_operators_used
        ! - initializing the coordinates
        ! - applying the initial conditions
        ! - initializing the i/o operator io_operators_used
        subroutine ini(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          call this%check_inputs()

          call this%bc_operators_used%ini(this%pmodel_eq_used)
          call this%ini_coordinates()
          call this%apply_initial_conditions()
          call this%io_operators_used%ini()

        end subroutine ini


        subroutine check_inputs(this)

          implicit none

          class(field_abstract), intent(in) :: this


          if(ne.ne.this%pmodel_eq_used%get_eq_nb()) then
             stop 'ne is not correct considering the physical model'
          end if

          if(bc_size.ne.this%sd_operators_used%get_bc_size()) then
             stop 'bc_size is not correct considering spatial operator'
          end if
          
        end subroutine check_inputs


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


        !apply the initial conditions
        subroutine apply_initial_conditions(this)

          implicit none

          class(field_abstract), intent(inout) :: this
          
          call this%pmodel_eq_used%apply_ic(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map)

        end subroutine apply_initial_conditions


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


        !write the data on an output file
        subroutine write_data(this, time)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind)          , intent(in)    :: time

          call this%io_operators_used%write_data(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%pmodel_eq_used,
     $         time)

        end subroutine write_data


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
        function get_nodes(this) result(nodes)

          implicit none

          class(field_abstract)           , intent(in)  :: this
          real(rkind), dimension(nx,ny,ne)              :: nodes

          nodes = this%nodes

        end function get_nodes


        !set the attribute nodes
        subroutine set_nodes(this,nodes)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes

          this%nodes = nodes

        end subroutine set_nodes


        !get the attribute x_map
        function get_x_map(this) result(x_map)

          implicit none

          class(field_abstract)     , intent(in)  :: this
          real(rkind), dimension(nx)              :: x_map

          x_map = this%x_map

        end function get_x_map


        !get the attribute y_map
        function get_y_map(this) result(y_map)

          implicit none

          class(field_abstract)     , intent(in)  :: this
          real(rkind), dimension(ny)              :: y_map

          y_map = this%y_map

        end function get_y_map      

      end module field_abstract_class
