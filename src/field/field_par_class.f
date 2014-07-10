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
      module field_par_class
      
        use field_class        , only : field
        use mpi
        use mpi_process_class  , only : mpi_process
        use parameters_constant, only : periodic_xy_choice,
     $                                  reflection_xy_choice,
     $                                  wall_xy_choice,
     $                                  wall_x_reflection_y_choice
        use parameters_input   , only : ntx,nty,nx,ny,npx,npy,bc_size,
     $                                  bc_choice,
     $                                  x_min,x_max,y_min,y_max
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: field_par


        !> @class field
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
        type, extends(field) :: field_par

          integer :: comm_2d
          integer :: usr_rank

          contains

          procedure, pass :: ini_cartesian_communicator
          procedure, pass :: ini_coordinates

        end type field_par


        contains


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
        !> 'field_par' object initialized
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
          class(field_par), intent(inout) :: this

          
          !< local variables
          type(mpi_process)     :: mpi_op
          integer, dimension(2) :: nb_tiles
          integer               :: comm_old
          integer               :: ndims
          logical, dimension(2) :: periods
          logical               :: reorganisation
          integer               :: ierror


          !< set the characteristics of the cartesian grid
          nb_tiles(1)    = npx
          nb_tiles(2)    = npy

          select case(bc_choice)

            case(periodic_xy_choice)
               periods(1)     = .true.
               periods(2)     = .true.

            case(reflection_xy_choice)
               periods(1)     = .false.
               periods(2)     = .false.

            case(wall_xy_choice)
               periods(1)     = .false.
               periods(2)     = .false.

            case(wall_x_reflection_y_choice)
               periods(1)     = .false.
               periods(2)     = .false.

            case default
               print *, 'field_par_class:'
               print *, 'bc_choice not implemented in'
               stop 'splitting the field into tiles'

          end select

          comm_old       = MPI_COMM_WORLD
          ndims          = 2
          reorganisation = .true.


          !< create Cartesian communicator
          call MPI_CART_CREATE(
     $         comm_old,
     $         ndims,
     $         nb_tiles,
     $         periods,
     $         reorganisation,
     $         this%comm_2d,
     $         ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'field_par_class: MPI_CART_CREATE failed'
          end if


          !< get the rank of the processor computing the tile
          call MPI_COMM_RANK(this%comm_2d, this%usr_rank, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             stop 'ini_cartesian_comm: MPI_COMM_RANK failed'
          end if

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

          class(field_par), intent(inout) :: this


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

      end module field_par_class
