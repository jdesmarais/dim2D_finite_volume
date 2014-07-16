      !> @file
      !> class extending the 'field' class to integrate mpi
      !> attributes that will allows the tile to communicate
      !> with its neighbours
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class extending the 'field' class to integrate mpi
      !> attributes that will allows the tile to communicate
      !> with its neighbours
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_abstract_par_class

        use bc_operators_par_class, only : bc_operators_par      
        use io_operators_par_class, only : io_operators_par
        use mpi
        use mpi_process_class     , only : mpi_process        
        use parameters_input      , only : ntx,nty,npx,npy,bc_size,
     $                                     bc_choice,nx,ny,ne,
     $                                     x_min,x_max,y_min,y_max
        use parameters_kind       , only : ikind, rkind
        use pmodel_eq_class       , only : pmodel_eq
        use sd_operators_class    , only : sd_operators
        use surrogate_class       , only : surrogate
        use td_operators_par_class, only : td_operators_par

        implicit none

        private
        public :: field_abstract_par


        !> @class field_abstract_par
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !>
        !> @param comm_2d
        !> attribute identifying the mpi main communicator between the
        !> tiles
        !>
        !> @param usr_rank
        !> attribute identifying the processor computing the tile
        !---------------------------------------------------------------
        type, extends(surrogate) :: field_abstract_par

          type(sd_operators)    , private :: sd_operators_used
          type(pmodel_eq)       , private :: pmodel_eq_used
          type(bc_operators_par), private :: bc_operators_used
          type(td_operators_par), private :: td_operators_used
          type(io_operators_par), private :: io_operators_used

          integer                         , private :: comm_2d
          integer                         , private :: usr_rank

          real(rkind), dimension(nx,ny,ne), private :: nodes
          real(rkind), dimension(nx)      , private :: x_map
          real(rkind), dimension(ny)      , private :: y_map
          real(rkind)                     , private :: dx
          real(rkind)                     , private :: dy

          contains

          procedure, pass          :: ini
          procedure, pass, private :: check_inputs
          procedure, pass, private :: ini_cartesian_communicator
          procedure, pass, private :: ini_coordinates

          procedure, pass, private :: apply_initial_conditions
          procedure, pass          :: compute_time_dev
          procedure, pass          :: apply_bc_on_nodes
          procedure, pass          :: compute_integration_step
          procedure, pass          :: write_data

          procedure, pass          :: get_comm_2d  !only for tests
          procedure, pass          :: get_usr_rank !only for tests
          procedure, pass          :: get_nodes    !only for tests
          procedure, pass          :: get_x_map    !only for tests
          procedure, pass          :: get_y_map    !only for tests      

        end type field_abstract_par


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


        !initialize the field_abstract_par object:
        ! - check the inputs
        ! - initialize the mpi process
        ! - initialize the cartesian communicator between
        !   the tiles
        ! - initialize the coordinates
        ! - initialize the boundary conditions
        ! - initialize the i/o operators
        subroutine ini(this)

          implicit none

          class(field_abstract_par), intent(inout) :: this

          type(mpi_process) :: mpi_op

          call this%check_inputs()
          call mpi_op%ini_mpi()
          call this%ini_cartesian_communicator()
          call this%ini_coordinates()
          call this%bc_operators_used%ini(
     $         this%comm_2d, this%pmodel_eq_used)
          call this%io_operators_used%ini(this%comm_2d, this%usr_rank)

        end subroutine ini


        !check the inputs for initialization of the field:
        ! - check the physical model: number of equations
        ! - check the size of the boundary layer for the
        !   space discretization operator
        subroutine check_inputs(this)

          implicit none

          class(field_abstract_par), intent(in) :: this


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
        !> subroutine initializing the comm_2d and usr_rank
        !> attributes of the tile
        !
        !> @date
        !> 21_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> 'field_abstract_par' object initialized
        !
        !>@param nb_tiles
        !> number of tiles in the x and y directions
        !
        !>@param bc_chosen
        !> boundary conditions chosen
        !--------------------------------------------------------------
        subroutine ini_cartesian_communicator(this)

          implicit none

          !< dummy arguments
          class(field_abstract_par), intent(inout) :: this

          type(mpi_process) :: mpi_op

          call mpi_op%ini_cartesian_communicator(
     $         this%comm_2d,
     $         this%usr_rank)

        end subroutine ini_cartesian_communicator


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

          class(field_abstract_par), intent(inout) :: this


          type(mpi_process)     :: mpi_op
          integer               :: dims_nb
          integer, dimension(2) :: cart_coord
          integer(ikind)        :: i,j
          integer               :: ierror
          real(rkind)           :: x_min_tile
          real(rkind)           :: y_min_tile


          !< initialize the space steps along the 
          !> x and y directions
          if(npx.gt.1) then
             this%dx = (x_max-x_min)/(npx*(nx-2*bc_size)-1)
          else
             this%dx = (x_max-x_min)/(ntx-1-2*bc_size)
          end if

          if(npy.gt.1) then
             this%dy = (y_max-y_min)/(npy*(ny-2*bc_size)-1)
          else
             this%dy = (y_max-y_min)/(nty-1-2*bc_size)
          end if


          !< get the cartesian coordinates of the tile
          dims_nb=2
          call MPI_CART_COORDS(
     $         this%comm_2d, this%usr_rank,
     $         dims_nb, cart_coord,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'ini_coordinates: MPI_CART_RANK failed'
          end if


          !< get the x_min corresponding to the tile
          if(cart_coord(1).eq.0) then
             x_min_tile = x_min
          else
             x_min_tile =
     $            x_min + cart_coord(1)*(nx-2*bc_size)*this%dx
          end if


          !< get the y_min corresponding to the tile
          if(cart_coord(2).eq.0) then
             y_min_tile = y_min
          else
             y_min_tile =
     $            y_min + cart_coord(2)*(ny-2*bc_size)*this%dy
          end if


          !< initialize the coordinates along the
          !> x-direction
          do i=1, nx
             this%x_map(i)=x_min_tile + (i-1-bc_size)*this%dx
          end do


          !< initialize the coordinates along the
          !> y-direction
          do j=1, ny
             this%y_map(j)=y_min_tile + (j-1-bc_size)*this%dy
          end do

        end subroutine ini_coordinates


        !initialize the nodes using the initial conditions
        !of the physical model
        subroutine apply_initial_conditions(this)

          implicit none

          class(field_abstract_par), intent(inout) :: this

          call this%pmodel_eq_used%apply_ic(
     $         this%nodes,
     $         this%x_map,
     $         this%y_map)

        end subroutine apply_initial_conditions


        !compute the time derivatives
        function compute_time_dev(this) result(time_dev)

          implicit none

          class(field_abstract_par), intent(in) :: this
          real(rkind), dimension(nx,ny,ne)      :: time_dev

          !make use of the time discretization operator
          !to compute the time derivative of the field
          time_dev = this%td_operators_used%compute_time_dev(
     $         this%comm_2d,
     $         this%usr_rank,
     $         this%nodes,
     $         this%dx,
     $         this%dy,
     $         this%sd_operators_used,
     $         this%pmodel_eq_used,
     $         this%bc_operators_used)

        end function compute_time_dev


        !apply the boundary conditions on the nodes
        subroutine apply_bc_on_nodes(this)

          implicit none

          class(field_abstract_par), intent(inout) :: this
          
          call this%bc_operators_used%apply_bc_on_nodes(
     $         this%comm_2d,
     $         this%usr_rank,
     $         this%nodes)

        end subroutine apply_bc_on_nodes


        !compute the time derivatives
        subroutine compute_integration_step(
     $     this, dt, nodes_tmp, time_dev, integration_step)

          implicit none

          class(field_abstract_par)       , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          procedure(integration_step_proc) :: integration_step

          call integration_step(this%nodes, dt, nodes_tmp, time_dev)

        end subroutine compute_integration_step


        !write the data on output files
        subroutine write_data(this, time)

          implicit none

          class(field_abstract_par), intent(inout) :: this
          real(rkind)              , intent(in)    :: time

          call this%io_operators_used%write_data(
     $         this%comm_2d,
     $         this%nodes,
     $         this%x_map,
     $         this%y_map,
     $         this%pmodel_eq_used,
     $         time)

        end subroutine write_data


        !get the comm_2d attribute      
        function get_comm_2d(this) result(comm_2d)

          implicit none

          class(field_abstract_par), intent(in) :: this
          integer                               :: comm_2d
          
          comm_2d = this%comm_2d

        end function get_comm_2d


        !get the usr_rank attribute      
        function get_usr_rank(this) result(usr_rank)

          implicit none

          class(field_abstract_par), intent(in) :: this
          integer                               :: usr_rank
          
          usr_rank = this%usr_rank

        end function get_usr_rank

        !get the nodes attribute
        function get_nodes(this) result(nodes)

          implicit none

          class(field_abstract_par), intent(in) :: this
          real(rkind), dimension(nx,ny,ne)      :: nodes
          
          nodes = this%nodes

        end function get_nodes


        !get the x_map attribute
        function get_x_map(this) result(x_map)

          implicit none

          class(field_abstract_par), intent(in) :: this
          real(rkind), dimension(nx)            :: x_map
          
          x_map = this%x_map

        end function get_x_map


        !get the y_map attribute
        function get_y_map(this) result(y_map)

          implicit none

          class(field_abstract_par), intent(in) :: this
          real(rkind), dimension(ny)            :: y_map
          
          y_map = this%y_map

        end function get_y_map

      end module field_abstract_par_class
