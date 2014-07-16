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
      module io_operators_par_class

        use pmodel_eq_class, only :
     $     pmodel_eq

        use io_operators_module, only :
     $       get_filename

        use io_operators_abstract_par_class, only :
     $       io_operators_abstract_par

        use mpi

        use mpi_process_class, only :
     $       mpi_process

        use netcdf

        use nf90_operators_module, only :
     $       nf90_handle_err,
     $       nf90_write_header,
     $       nf90_def_var_model,
     $       nf90_put_var_model,
     $       nf90_close_file

        use parameters_input, only :
     $       nx,ny,ne,npx,npy,bc_size

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: io_operators_par


        !> @class io_operators_par
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
        type, extends(io_operators_abstract_par) :: io_operators_par

          contains

          procedure, pass :: write_data

        end type io_operators_par


        contains


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
        subroutine write_data(
     $       this, comm_2d,
     $       nodes, x_map, y_map,
     $       p_model, time)

          implicit none

          class(io_operators_par)         , intent(inout) :: this
          integer                         , intent(in)    :: comm_2d
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: time

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
     $         comm=comm_2d,
     $         info=MPI_INFO_NULL)
          call nf90_handle_err(retval)


          !<write the header of the file
          call nf90_write_header(ncid,p_model)


          !<define the variables saved in the file
          call nf90_def_var_model(ncid,p_model,coord_id,data_id)


          !<put the variables in the file
          call nf90_put_var_model(
     $         ncid,
     $         coord_id, data_id,
     $         time, nodes, x_map, y_map,
     $         start=this%start,
     $         count=this%count,
     $         limit=this%limit)


          !<close the netcdf file
          call nf90_close_file(ncid)

          !<increment the internal counter
          call this%increment_counter()

        end subroutine write_data

      end module io_operators_par_class
