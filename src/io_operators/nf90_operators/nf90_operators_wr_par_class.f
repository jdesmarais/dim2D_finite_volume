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
      module nf90_operators_wr_par_class

        use cg_operators_class   , only : cg_operators
        use dim2d_eq_class       , only : dim2d_eq
        use field_par_class      , only : field_par
        use io_operators_module  , only : get_filename
        use mpi
        use mpi_process_class    , only : mpi_process
        use netcdf
        use nf90_operators_module, only : nf90_handle_err,
     $                                    nf90_write_header,
     $                                    nf90_def_var_model,
     $                                    nf90_put_var_model
        use parameters_input     , only : nx,ny,ne,npx,npy
        use parameters_kind      , only : rkind

        implicit none

        private
        public :: nf90_operators_wr_par


        !> @class nf90_operators_wr_par
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
        type :: nf90_operators_wr_par

          integer :: nb_timesteps_written
          
          integer, dimension(2) :: start
          integer, dimension(2) :: count
          integer, dimension(4) :: limit

          contains

          procedure, pass, non_overridable :: initialize
          procedure, pass, non_overridable :: write_data

        end type nf90_operators_wr_par


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
        subroutine initialize(this, field_used, sd_op, counter)

          implicit none

          class(nf90_operators_wr_par), intent(inout) :: this
          class(field_par)            , intent(in)    :: field_used
          type(cg_operators)          , intent(in)    :: sd_op
          integer           , optional, intent(in)    :: counter


          integer               :: bc_size
          type(mpi_process)     :: mpi_op
          integer, dimension(2) :: cart_coords
          integer               :: ierror


          !< intialize the 'nb_timesteps_written' attribute
          if(present(counter)) then
             this%nb_timesteps_written=counter
          else
             this%nb_timesteps_written=0
          end if


          !< initialize the 'start', 'count' and 'limit' attributes
          
          !< get the cartesian coordinates of the tile
          call MPI_CART_COORDS(
     $         field_used%comm_2d, field_used%usr_rank,
     $         2, cart_coords, ierror)
          if(ierror.ne.MPI_SUCCESS) then
             call mpi_op%finalize_mpi()
             print '(''nf90_operators_wr_par: initialize:'')'
             stop ' MPI_CART_COORDS failed'
          end if


          !< get the size of the boundary layers
          bc_size = sd_op%get_bc_size()


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
          if((cart_coords(1).eq.0).or.(cart_coords(1).eq.(npx-1))) then
             this%count(1)=this%count(1)+bc_size
          end if

          if((cart_coords(2).eq.0).or.(cart_coords(2).eq.(npy-1))) then
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

        end subroutine initialize


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to write data encapsulated in field_used
        !> on netcdf output files
        !
        !> @date
        !> 28_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param p_model
        !> physical model
        !
        !>@param time
        !> reduced simulation time
        !--------------------------------------------------------------
        subroutine write_data(this,f_used,p_model,bc_size,time)

          implicit none

          class(nf90_operators_wr_par), intent(inout) :: this
          class(field_par)            , intent(in)    :: f_used
          type(dim2d_eq)              , intent(in)    :: p_model
          integer                     , intent(in)    :: bc_size
          real(rkind)                 , intent(in)    :: time

          integer                :: ncid
          integer                :: retval
          integer, dimension(3)  :: coord_id
          integer, dimension(ne) :: data_id
          character(len=16)      :: filename


          !<get the name for the netcdf file
          call get_filename(filename, this%nb_timesteps_written)


          !<create the netcdf file with parallel file access
          retval = NF90_CREATE(
     $         trim(filename),
     $         IOR(NF90_NETCDF4, NF90_MPIPOSIX),
     $         ncid,
     $         comm=f_used%comm_2d,
     $         info=MPI_INFO_NULL)
          call nf90_handle_err(retval)


          !<write the header of the file
          call nf90_write_header(ncid,p_model)


          !<define the variables saved in the file
          call nf90_def_var_model(ncid,p_model,bc_size,coord_id,data_id)


          !<put the variables in the file
          call nf90_put_var_model(
     $         ncid,
     $         coord_id, data_id,
     $         time, f_used,
     $         start=this%start,
     $         count=this%count,
     $         limit=this%limit)


          !<close the netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)

          !<increment the internal counter
          this%nb_timesteps_written=this%nb_timesteps_written+1

        end subroutine write_data

      end module nf90_operators_wr_par_class
