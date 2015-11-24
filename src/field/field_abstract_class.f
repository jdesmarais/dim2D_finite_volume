      !> @file
      !> class encapsulating the main tables for the variables
      !> (3D array: nx x ny x ne) and the coordinate maps (1D array:
      !> nx for x-dir, 1D array: ny for y-dir) of the Cartesian mesh
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables
      !> (3D array: nx x ny x ne) and the coordinate maps (1D array:
      !> nx for x-dir, 1D array: ny for y-dir) of the Cartesian mesh
      !
      !> @date
      ! 07_08_2013 - initial version               - J.L. Desmarais
      ! 14_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_abstract_class
      
        use bc_operators_gen_class, only :
     $       bc_operators_gen

        use bc_operators_module, only :
     $       shall_bc_on_nodes_be_applied

        use cmd_operators_class, only :
     $       cmd_operators

        use interface_integration_step, only :
     $       timeInt_step,
     $       timeInt_step_nopt

        use io_operators_class, only :
     $       io_operators

        use parameters_constant, only :
     $       N,S,E,W,
     $       bc_timedev_choice

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       x_min,x_max,
     $       y_min,y_max,
     $       bc_choice,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice,
     $       steady_state_limit

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use surrogate_class, only :
     $       surrogate

        use td_operators_class, only :
     $       td_operators

        implicit none


        private
        public :: field_abstract


        !>@class field_abstract
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !---------------------------------------------------------------
        type, extends(surrogate) :: field_abstract

          type(sd_operators)     :: sd_operators_used !< @brief space discretization operators
          type(pmodel_eq)        :: pmodel_eq_used    !< @brief physical model
          type(bc_operators_gen) :: bc_operators_used !< @brief boundary conditions operators
          type(td_operators)     :: td_operators_used !< @brief time discretization operators
          type(io_operators)     :: io_operators_used !< @brief i/o operators

          real(rkind)                      :: time   !<@brief simulation time \f$ t \f$
          real(rkind), dimension(nx,ny,ne) :: nodes  !<@brief governing variables: e.g. \f$ \left(\rho, q_x, q_y, \rho E \right) \f$
          real(rkind), dimension(nx)       :: x_map  !<@brief space discretization map along x-direction: \f$(x_1, \cdots, x_{\textrm{nx}})\f$
          real(rkind), dimension(ny)       :: y_map  !<@brief space discretization map along y-direction: \f$(y_1, \cdots, y_{\textrm{ny}})\f$
          real(rkind)                      :: dx     !<@brief grid-size along x-direction: \f$ dx \f$
          real(rkind)                      :: dy     !<@brief grid-size along y-direction: \f$ dy \f$

          integer(ikind), dimension(2) :: x_borders  !<@brief time integration borders along the x-direction \f$[i_{\textrm{min}},i_{\textrm{max}}] \f$
          integer(ikind), dimension(2) :: y_borders  !<@brief time integration borders along the y-direction \f$[j_{\textrm{min}},j_{\textrm{max}}] \f$

          contains

          procedure, pass          :: ini                             !<@brief initialize the interior domain and its operators
          procedure, pass          :: check_inputs                    !<@brief check the simulation inputs
          procedure, pass, private :: ini_coordinates                 !<@brief initialize the space discretization maps for the computational domain
          procedure, pass          :: ini_for_timeInt                 !<@brief initialize the integration borders

          procedure, pass          :: apply_initial_conditions        !<@brief apply the initial conditions of the physical model on the computational domain
          procedure, pass          :: compute_time_dev                !<@brief compute the time derivatives of the computational domain
          procedure, pass          :: compute_time_dev_ext            !<@brief compute the time derivatives of the computational domain and its extension if any
          procedure, pass          :: apply_bc_on_nodes               !<@brief apply the boundary conditions on the computational domain
          procedure, pass          :: compute_integration_step        !<@brief compute the integration step of the interior domain
          procedure, pass          :: compute_integration_step_ext    !<@brief compute the integration step of the interior domain and its extension if any (only interface)
          procedure, pass          :: write_data                      !<@brief write the interior domain data on output files

          procedure, pass          :: check_steady_state              !<@brief check whether the steady state is reached

          procedure, pass          :: adapt_domain                    !<@brief check whether the computational domain should be extended

          procedure, pass          :: set_dx                          !<@brief set the grid size along the x-axis (only for tests)
          procedure, pass          :: set_dy                          !<@brief set the grid size along the y-axis (only for tests)
          procedure, pass          :: get_nodes                       !<@brief get the array containing the grid point data (only for tests)
          procedure, pass          :: set_nodes                       !<@brief set the array containing the grid point data (only for tests)
          procedure, pass          :: get_x_map                       !<@brief get the array containing the space discretization map along the x-axis (only for tests)
          procedure, pass          :: get_y_map                       !<@brief get the array containing the space discretization map along the y-axis (only for tests)
          procedure, pass          :: set_x_map                       !<@brief set the array containing the space discretization map along the x-axis (only for tests)
          procedure, pass          :: set_y_map                       !<@brief set the array containing the space discretization map along the y-axis (only for tests)
          procedure, pass          :: get_time                        !<@brief get the time corresponding to the variables

        end type field_abstract


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the field_abstract by:
        !> -# check the simulation inputs
        !> -# initializing the boundary conditions bc_operators_used
        !> -# initializing the coordinates
        !> -# applying the initial conditions
        !> -# initializing the i/o operators io_operators_used
        !        
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(field_abstract)  , intent(inout) :: this


          type(cmd_operators) :: cmd_operators_used
          real(rkind), dimension(:,:,:), allocatable :: nodes_tmp


          !analyze the command line arguments
          !========================================
          call cmd_operators_used%analyse_cmd_line_arg()


          !initialize the computational domain
          !========================================
          !1) initialize the integration borders
          call this%ini_for_timeInt()

          !2) initialize the time+x_map,y_map+nodes+io_operators
          if(cmd_operators_used%is_restart_activated()) then

             call this%io_operators_used%read_data(
     $            trim(cmd_operators_used%get_restart_filename()),
     $            this%nodes,
     $            this%x_map,
     $            this%y_map,
     $            this%pmodel_eq_used,
     $            this%time)
             
          else
             
             this%time = 0.0

             call this%ini_coordinates()
             call this%apply_initial_conditions()
             call this%io_operators_used%ini()

          end if

          allocate(nodes_tmp(nx,ny,ne))
          nodes_tmp = this%nodes
          call this%apply_bc_on_nodes(nodes_tmp)
          deallocate(nodes_tmp)

          !3) verify the inputs
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
          
        end subroutine check_inputs


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the integration borders when the governing
        !> variables are integrated in time (check whether a PDE
        !> is solved for the boundary points or not)
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main governing variables
        !--------------------------------------------------------------
        subroutine ini_for_timeInt(this)

          implicit none

          class(field_abstract), intent(inout) :: this


          this%x_borders = [bc_size+1,nx-bc_size]
          this%y_borders = [bc_size+1,ny-bc_size]

          if(bc_W_type_choice.eq.bc_timedev_choice) then
             this%x_borders(1) = this%x_borders(1)-bc_size
          end if

          if(bc_E_type_choice.eq.bc_timedev_choice) then
             this%x_borders(2) = this%x_borders(2)+bc_size
          end if
          
          if(bc_S_type_choice.eq.bc_timedev_choice) then
             this%y_borders(1) = this%y_borders(1)-bc_size
          end if

          if(bc_N_type_choice.eq.bc_timedev_choice) then
             this%y_borders(2) = this%y_borders(2)+bc_size
          end if

        end subroutine ini_for_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the space discretization maps for the
        !> computational domain
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
        !> the computational domain
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
        !> compute the time derivatives of the computational domain
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
          !to compute the time derivative of the computational
          !domain
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%time,
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used)

        end function compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the computational domain
        !> and its extension if any
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
          !to compute the time derivative of the computational
          !domain
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%time,
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
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
        !
        !>@param nodes_tmp
        !> governing variables of the computational domain
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,nodes_tmp)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp

          if(shall_bc_on_nodes_be_applied()) then

             call this%bc_operators_used%apply_bc_on_nodes(
     $            this%time,this%x_map,this%y_map,
     $            nodes_tmp,this%pmodel_eq_used,
     $            this%nodes)

          end if

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step for the computational
        !> domain corresponding to the interior domain
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
     $     this,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     integration_step)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step)                         :: integration_step


          call integration_step(
     $         this%nodes, dt, nodes_tmp, time_dev,
     $         x_borders=this%x_borders,
     $         y_borders=this%y_borders)
          
        end subroutine compute_integration_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step for the computational
        !> domain of the interior domain and its extension if any
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
     $     this, dt,
     $     nodes_tmp,
     $     time_dev,
     $     integration_step,
     $     integration_step_nopt)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(timeInt_step)                         :: integration_step
          procedure(timeInt_step_nopt)                    :: integration_step_nopt

          
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
        !--------------------------------------------------------------
        subroutine write_data(this)

          implicit none

          class(field_abstract), intent(inout) :: this

          call this%io_operators_used%write_data(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%pmodel_eq_used,
     $         this%time)

        end subroutine write_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the computational domain
        !
        !> @date
        !> 14_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables at t
        !
        !>@return steady_state
        !> check whether the steady state is reached for
        !> the simulation: maximum time dev below a threshold
        !--------------------------------------------------------------
        function check_steady_state(this) result(steady_state)

          implicit none

          class(field_abstract), intent(in) :: this
          logical                           :: steady_state

          real(rkind), dimension(nx,ny,ne)  :: time_dev

          integer :: i,j,k

          !make use of the time discretization operator
          !to compute the time derivative of the computational
          !domain
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%time,
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used)

          steady_state = .true.

          do k=1, ne
             do j=1, ny
                do i=1, nx

                   if(abs(time_dev(i,j,k)).gt.(steady_state_limit)) then
                      steady_state = .false.
                      return
                   end if

                end do
             end do
          end do

        end function check_steady_state


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the computational domain
        !
        !> @date
        !> 14_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables at t
        !
        !>@param nodes_tmp
        !> nodes at the previous time step (t-dt)
        !
        !>@param dt
        !> time step
        !--------------------------------------------------------------
        subroutine adapt_domain(this,dt,nodes0)

          implicit none

          class(field_abstract)           , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes0

          real(rkind) :: node_s
          real(rkind) :: dx_s
          real(rkind) :: dt_s

          print '(''field_abstract_class'')'
          print '(''adapt_domain'')'
          stop 'not implemented'

          dx_s   = this%dx
          dt_s   = dt
          node_s = nodes0(1,1,1)

        end subroutine adapt_domain


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
        !> governing variables of the computational domain
        !> in the interior domain
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
        !> governing variables of the computational domain
        !> in the interior domain
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

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the array containing the space discretization map
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
        subroutine set_x_map(this,x_map)

          implicit none

          class(field_abstract)     , intent(inout):: this
          real(rkind), dimension(nx), intent(in)   :: x_map

          this%x_map = x_map

        end subroutine set_x_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the array containing the space discretization map
        !> along the y-axis
        !
        !> @date
        !> 02_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return y_map
        !> discretisation map along the y-axis
        !--------------------------------------------------------------
        subroutine set_y_map(this,y_map)

          implicit none

          class(field_abstract)     , intent(inout) :: this
          real(rkind), dimension(ny), intent(in)    :: y_map

          this%y_map = y_map

        end subroutine set_y_map


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the time of the simulation
        !
        !> @date
        !> 06_12_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@return time
        !> simulation time
        !--------------------------------------------------------------
        function get_time(this) result(time)

          implicit none

          class(field_abstract), intent(in) :: this
          real(rkind)                       :: time

          time = this%time

        end function get_time        

      end module field_abstract_class
