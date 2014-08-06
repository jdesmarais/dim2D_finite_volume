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
      
        use bc_operators_class        , only : bc_operators
        use interface_integration_step, only : timeInt_step,
     $                                         timeInt_step_nopt
        use io_operators_class        , only : io_operators
        use parameters_constant       , only : periodic_xy_choice,
     $                                         reflection_xy_choice,
     $                                         wall_xy_choice,
     $                                         wall_x_reflection_y_choice,
     $                                         hedstrom_xy_choice,
     $                                         hedstrom_xy_corners_choice,
     $                                         hedstrom_x_reflection_y_choice
        use parameters_input          , only : nx,ny,ne,bc_size,
     $                                         x_min,x_max,
     $                                         y_min,y_max,
     $                                         bc_choice,
     $                                         bcx_type_choice,
     $                                         bcy_type_choice
        use parameters_kind           , only : ikind, rkind
        use pmodel_eq_class           , only : pmodel_eq
        use sd_operators_class        , only : sd_operators
        use surrogate_class           , only : surrogate
        use td_operators_class        , only : td_operators

        implicit none


        private
        public :: field_abstract


        !>@class field_abstract
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !
        !>@param sd_operators_used
        !> space discretization operators
        !
        !>@param pmodel_eq_used
        !> physical model governing equations
        !
        !>@param bc_operators_used
        !> boundary conditions operators
        !
        !>@param td_operators_used
        !> time discretization operators
        !
        !>@param td_integrator_used
        !> time integration operators
        !
        !>@param io_integrator_used
        !> i/o operators
        !
        !>@param nodes
        !> variables computed during the simulation
        !> (ex: mass=nodes(:,:,1),
        !> momentum_x=nodes(:,:,2),
        !> momentum_y=nodes(:,:,3),
        !> energy=nodes(:,:,4))
        !
        !>@param x_map
        !> discretisation map along the x-axis
        !
        !>@param y_map
        !> discretisation map along the y-axis
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param ini
        !> initialize the interior domain and its operators
        !
        !>@param check_inputs
        !> check the simulation inputs
        !
        !>@param ini_coordinates
        !> initialize the space discretization maps for the field
        !
        !>@param apply_initial_conditions
        !> apply the initial conditions of the physical model on
        !> the field
        !
        !>@param compute_time_dev
        !> compute the time derivatives of the field
        !
        !>@param compute_time_dev_ext
        !> compute the time derivatives of the field and its extension
        !> if any
        !
        !>@param apply_bc_on_nodes
        !> apply the boundary conditions on the grid points
        !
        !>@param compute_integration_step
        !> compute the integration step of the interior domain
        !
        !>@param compute_integration_step_ext
        !> compute the integration step of the interior domain and
        !> its extension if any
        !
        !>@param write_data
        !> write the interior domain data on output files depending
        !> on the i/o operators used
        !
        !>@param set_dx
        !> set the grid size along the x-axis
        !
        !>@param set_dy
        !> set the grid size along the y-axis
        !
        !>@param get_nodes
        !> get the array containing the grid point data
        !
        !>@param set_nodes
        !> set the array containing the grid point data
        !
        !>@param get_x_map
        !> get the array containing the space discretization map
        !> along the x-axis
        !
        !>@param get_y_map
        !> get the array containing the space discretization map
        !> along the y-axis
        !---------------------------------------------------------------
        type, extends(surrogate) :: field_abstract

          type(sd_operators), private :: sd_operators_used
          type(pmodel_eq)   , private :: pmodel_eq_used
          type(bc_operators), private :: bc_operators_used
          type(td_operators), private :: td_operators_used
          type(io_operators), private :: io_operators_used

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind)                      :: dx
          real(rkind)                      :: dy

          contains

          procedure, pass          :: ini
          procedure, pass, private :: check_inputs
          procedure, pass, private :: ini_coordinates
          procedure, pass, private :: apply_initial_conditions
          procedure, pass          :: compute_time_dev
          procedure, pass          :: compute_time_dev_ext
          procedure, pass          :: apply_bc_on_nodes
          procedure, pass          :: compute_integration_step
          procedure, pass          :: compute_integration_step_ext
          procedure, pass          :: write_data

          procedure, pass          :: set_dx    !only for tests
          procedure, pass          :: set_dy    !only for tests
          procedure, pass          :: get_nodes !only for tests
          procedure, pass          :: set_nodes !only for tests
          procedure, pass          :: get_x_map !only for tests
          procedure, pass          :: get_y_map !only for tests

        end type field_abstract


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the field_abstract by:
        !> (1) check the simulation inputs
        !> (2) initializing the boundary conditions bc_operators_used
        !> (3) initializing the coordinates
        !> (4) applying the initial conditions
        !> (5) initializing the i/o operators io_operators_used
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          call this%bc_operators_used%ini(this%pmodel_eq_used)
          call this%ini_coordinates()
          call this%apply_initial_conditions()
          call this%io_operators_used%ini()

          call this%check_inputs()

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the inputs of the simulation
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine check_inputs(this)

          implicit none

          class(field_abstract), intent(in) :: this


          if(ne.ne.this%pmodel_eq_used%get_eq_nb()) then
             stop 'ne is not correct considering the physical model'
          end if

          if(bc_size.ne.this%sd_operators_used%get_bc_size()) then
             stop 'bc_size is not correct considering spatial operator'
          end if

          if(bcx_type_choice.ne.this%bc_operators_used%get_bcx_type()) then
             stop 'bcx_type_choice does not match bc_operator'
          end if

          if(bcy_type_choice.ne.this%bc_operators_used%get_bcy_type()) then
             stop 'bcy_type_choice does not match bc_operator'
          end if
          
        end subroutine check_inputs


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the space discretization map for the field
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions of the physical model on
        !> the field
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine apply_initial_conditions(this)

          implicit none

          class(field_abstract), intent(inout) :: this
          
          call this%pmodel_eq_used%apply_ic(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map)

        end subroutine apply_initial_conditions


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the field
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !
        !>@return time_dev
        !> array containing the time derivatives of the grid points
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the field and its extension
        !> if any
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !
        !>@return time_dev
        !> array containing the time derivatives of the interior grid
        !> points
        !--------------------------------------------------------------
        function compute_time_dev_ext(this) result(time_dev)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind), dimension(nx,ny,ne)  :: time_dev

          
          print '(''********************************'')'
          print '(''field_abstract'')'
          print '(''compute_time_dev_ext'')'
          print '(''this procedure is only for field'')'
          print '(''with domain extension'')'
          print '(''********************************'')'
          print '(''subroutine not implemented'')'
          print '(''********************************'')'
          stop

          !make use of the time discretization operator
          !to compute the time derivative of the field
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%nodes,
     $         this%dx,
     $         this%dy,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used)

        end function compute_time_dev_ext


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the boundary conditions on the grid points
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          if((bc_choice.ne.hedstrom_xy_choice).and.
     $       (bc_choice.ne.hedstrom_xy_corners_choice)) then

             call this%bc_operators_used%apply_bc_on_nodes(this%nodes)

          end if


        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step of the interior domain
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array containing the temporary grid points for the
        !> time integration of the interior computational domain
        !
        !>@param time_dev
        !> time derivatives of the interior domain
        !
        !>@param integration_step
        !> procedure for the time integration of the interior domain
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, nodes_tmp, time_dev, integration_step)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step) :: integration_step


          integer(ikind), dimension(2) :: x_borders
          integer(ikind), dimension(2) :: y_borders

          
          select case(bc_choice)

            case(reflection_xy_choice, periodic_xy_choice,
     $           wall_xy_choice, wall_x_reflection_y_choice)
               x_borders=[bc_size+1,nx-bc_size]
               y_borders=[bc_size+1,ny-bc_size]

            case(hedstrom_xy_choice,hedstrom_xy_corners_choice)
               x_borders=[1,nx]
               y_borders=[1,ny]

            case(hedstrom_x_reflection_y_choice)
               x_borders=[1,nx]
               y_borders=[bc_size+1,ny-bc_size]

            case default
               print '(''field_abstract: compute_integration_step'')'
               stop 'bc not implemented'

          end select

          call integration_step(
     $         this%nodes, dt, nodes_tmp, time_dev,
     $         x_borders=x_borders,
     $         y_borders=y_borders)

        end subroutine compute_integration_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step of the interior domain and its
        !> extension if any
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array containing the temporary grid points for the
        !> time integration of the interior computational domain
        !
        !>@param time_dev
        !> time derivatives of the interior domain
        !
        !>@param integration_step
        !> procedure for the time integration of the interior domain
        !
        !>@param integration_step_nopt
        !> procedure for the time integration of the domain extension
        !--------------------------------------------------------------
        subroutine compute_integration_step_ext(
     $     this, dt, nodes_tmp, time_dev,
     $     integration_step, integration_step_nopt)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step)      :: integration_step
          procedure(timeInt_step_nopt) :: integration_step_nopt

          
          integer, dimension(nx,ny) :: grdpts_id

          
          print '(''********************************'')'
          print '(''field_abstract'')'
          print '(''compute_integration_step_ext'')'
          print '(''this procedure is only for field'')'
          print '(''with domain extension'')'
          print '(''********************************'')'
          print '(''subroutine not implemented'')'
          print '(''********************************'')'
          stop

          call integration_step(
     $         this%nodes, dt, nodes_tmp, time_dev)

          call integration_step_nopt(
     $         this%nodes, dt, nodes_tmp, time_dev, grdpts_id)

        end subroutine compute_integration_step_ext


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> write the interior domain data on output files depending
        !> on the i/o operators used
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param time
        !> simulation time
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the grid size along the x-axis
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !--------------------------------------------------------------
        subroutine set_dx(this, dx)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind)          , intent(in)    :: dx

          this%dx = dx

        end subroutine set_dx

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the grid size along the x-axis
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param dy
        !> grid size along the y-axis
        !--------------------------------------------------------------
        subroutine set_dy(this, dy)

          implicit none

          class(field_abstract), intent(inout) :: this
          real(rkind)          , intent(in)    :: dy

          this%dy = dy

        end subroutine set_dy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the array containing the grid point data
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return nodes
        !> grid size along the y-axis
        !--------------------------------------------------------------
        function get_nodes(this) result(nodes)

          implicit none

          class(field_abstract)           , intent(in)  :: this
          real(rkind), dimension(nx,ny,ne)              :: nodes

          nodes = this%nodes

        end function get_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the array containing the grid point data
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param nodes
        !> grid size along the y-axis
        !--------------------------------------------------------------
        subroutine set_nodes(this,nodes)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes

          this%nodes = nodes

        end subroutine set_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the array containing the space discretization map
        !> along the x-axis
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return x_map
        !> discretisation map along the x-axis
        !--------------------------------------------------------------
        function get_x_map(this) result(x_map)

          implicit none

          class(field_abstract)     , intent(in)  :: this
          real(rkind), dimension(nx)              :: x_map

          x_map = this%x_map

        end function get_x_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the array containing the space discretization map
        !> along the y-axis
        !
        !> @date
        !> 18_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return y_map
        !> discretisation map along the y-axis
        !--------------------------------------------------------------
        function get_y_map(this) result(y_map)

          implicit none

          class(field_abstract)     , intent(in)  :: this
          real(rkind), dimension(ny)              :: y_map

          y_map = this%y_map

        end function get_y_map      

      end module field_abstract_class
