      !> @file
      !> class encapsulating subroutines to
      !> write netcdf files
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to
      !> write netcdf files
      !
      !> @date
      !> 28_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module io_operators_abstract_par_class

        use pmodel_eq_class      , only : pmodel_eq
        use mpi
        use mpi_process_class    , only : mpi_process
        use netcdf
        use parameters_input     , only : nx,ny,ne,npx,npy,bc_size
        use parameters_kind      , only : rkind

        implicit none

        private
        public :: io_operators_abstract_par


        !> @class io_operators_abstract_par
        !> class encapsulating subroutines to write netcdf
        !> files in a parallel file management system and 
        !> on a parallel memory distributed system
        !>
        !> @param nb_timesteps_written
        !> attribute saving the total number of timesteps
        !> already written previously to name the output file
        !> correctly
        !>
        !> @param start
        !> attribute saving the indices of the SW border of
        !> the tables written in the netcdf file
        !>
        !> @param count
        !> attribute saving the extends along the x and y
        !> directions of the tables written in the netcdf file
        !>
        !> @param limit
        !> attribute saving the extends of the local nodes
        !> written in the globla computational domain
        !>
        !> @param initialize
        !> initialize the counter to identify the netcdf file
        !> written and the position of the local nodes in the
        !> global domain computed
        !>
        !> @param write_data
        !> write the data encapsulated in field-type objects
        !> on netcdf output files
        !---------------------------------------------------------------
        type,abstract :: io_operators_abstract_par

          integer :: nb_timesteps_written
          
          integer, dimension(2) :: start
          integer, dimension(2) :: count
          integer, dimension(4) :: limit

          contains

          procedure            , pass           :: ini
          procedure            , pass, private  :: ini_writing_borders
          procedure(write_proc), pass, deferred :: write_data
          procedure            , pass           :: increment_counter

        end type io_operators_abstract_par


        interface

           subroutine write_proc(
     $          this, comm_2d,
     $          nodes, x_map, y_map,
     $          p_model, time)

             import io_operators_abstract_par
             import pmodel_eq
             import nx,ny,ne
             import rkind

             class(io_operators_abstract_par), intent(inout) :: this
             integer                         , intent(in)    :: comm_2d
             real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
             real(rkind), dimension(nx)      , intent(in)    :: x_map
             real(rkind), dimension(ny)      , intent(in)    :: y_map
             type(pmodel_eq)                 , intent(in)    :: p_model
             real(rkind)                     , intent(in)    :: time

           end subroutine write_proc

        end interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to initialize the counter identifying
        !> the netcdf file written
        !
        !> @date
        !> 28_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the counter
        !
        !>@param counter
        !> optional argument setting the counter
        !--------------------------------------------------------------
        subroutine ini(this, comm_2d, usr_rank, counter)

          implicit none

          class(io_operators_abstract_par), intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank
          integer      , optional         , intent(in)    :: counter

          !< intialize the 'nb_timesteps_written' attribute
          if(present(counter)) then
             this%nb_timesteps_written=counter
          else
             this%nb_timesteps_written=0
          end if

          !< initialize the 'start', 'count' and 'limit' attributes
          call this%ini_writing_borders(comm_2d, usr_rank)

        end subroutine ini


        !initialize the 'start', 'count' and 'limit' attributes
        subroutine ini_writing_borders(this, comm_2d, usr_rank)

          implicit none

          class(io_operators_abstract_par), intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          integer                         , intent(in)    :: usr_rank


          type(mpi_process)     :: mpi_op
          integer, dimension(2) :: cart_coords
          integer               :: ierror


          !< get the cartesian coordinates of the tile
          call MPI_CART_COORDS(
     $         comm_2d, usr_rank,
     $         2, cart_coords, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             print '(''io_operators_abstract_par: initialize:'')'
             stop ' MPI_CART_COORDS failed'
          end if


          !< initialize count
          this%count(1) = nx-2*bc_size
          this%count(2) = ny-2*bc_size


          !< initialize start
          this%start(1) = 1+bc_size+cart_coords(1)*(nx-2*bc_size)
          this%start(2) = 1+bc_size+cart_coords(2)*(ny-2*bc_size)


          !< initialize limit
          this%limit(1) = 1+bc_size
          this%limit(2) = nx-bc_size
          this%limit(3) = 1+bc_size
          this%limit(4) = ny-bc_size
          

          !< update 'count' depending on the cartesian position
          !> we add the boundary layer if the tile in on the
          !> border of the computational domain
          if(cart_coords(1).eq.0) then
             this%count(1)=this%count(1)+bc_size
          end if

          if(cart_coords(1).eq.(npx-1)) then
             this%count(1)=this%count(1)+bc_size
          end if

          if(cart_coords(2).eq.0) then
             this%count(2)=this%count(2)+bc_size
          end if

          if(cart_coords(2).eq.(npy-1)) then
             this%count(2)=this%count(2)+bc_size
          end if


          !< update 'start' depending on the cartesian position
          !> we begin at the lowest SW border of the tile if we
          !> include the boundary layer in printing the nodes
          if(cart_coords(1).eq.0) then
             this%start(1)=this%start(1)-bc_size
          end if

          if(cart_coords(2).eq.0) then
             this%start(2)=this%start(2)-bc_size
          end if


          !< update limit depending on the cartesian position
          !> depending if the tile in one the edge of the
          !> computational domain, we include the boundary layer
          if(cart_coords(1).eq.0) then
             this%limit(1)=this%limit(1)-bc_size
          end if

          if(cart_coords(1).eq.(npx-1)) then
             this%limit(2)=this%limit(2)+bc_size
          end if
          
          if(cart_coords(2).eq.0) then
             this%limit(3)=this%limit(3)-bc_size
          end if

          if(cart_coords(2).eq.(npy-1)) then
             this%limit(4)=this%limit(4)+bc_size
          end if          

        end subroutine ini_writing_borders


        !increment the counter for tracking the number of timesteps
        !written
        subroutine increment_counter(this)

          implicit none

          class(io_operators_abstract_par), intent(inout) :: this

          this%nb_timesteps_written=this%nb_timesteps_written+1

        end subroutine increment_counter

      end module io_operators_abstract_par_class
