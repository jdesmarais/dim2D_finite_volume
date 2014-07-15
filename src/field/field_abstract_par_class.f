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
      
        use surrogate_class    , only : surrogate
        use mpi_process_class  , only : mpi_process
        
        use parameters_input   , only : ntx,nty,nx,ny,npx,npy,bc_size,
     $                                  bc_choice,
     $                                  x_min,x_max,y_min,y_max
        use parameters_kind    , only : ikind, rkind

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

          type(mpi_process) :: mpi_process_used

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

          procedure, pass :: get_comm2d !only for tests

        end type field_abstract_par


        contains


        !initialize the field_abstract_par object:
        ! - initialize the cartesian communicator between
        !   the tiles
        ! - initialize the coordinates
        subroutine ini(this)

          implicit none

          class(field_abstract_par), intent(inout) :: this

          call this%check_inputs()
          call this%mpi_process_used%ini_mpi()
          call this%ini_cartesian_communicator()
          call this%ini_coordinates()

        end subroutine ini


        !check the inputs for initialization of the field:
        ! - check the physical model: number of equations
        ! - check the size of the boundary layer for the
        !   space discretization operator
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

          call this%mpi_process_used%ini_cartesian_communicator(
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


        !get the comm2d attribute
        function get_comm2d(this) result(comm2d)

          implicit none

          class(field_abstract_par), intent(in) :: this
          integer                               :: comm2d

          comm2d = this%comm_2d

        end function get_comm2d

      end module field_abstract_par_class
