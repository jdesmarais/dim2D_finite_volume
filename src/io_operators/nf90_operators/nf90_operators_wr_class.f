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
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module nf90_operators_wr_class

        use field_class          , only : field
        use io_operators_module  , only : get_filename
        use netcdf
        use nf90_operators_module, only : nf90_handle_err,
     $                                    nf90_write_header,
     $                                    nf90_def_var_model,
     $                                    nf90_put_var_model
        use parameters_input     , only : ne
        use parameters_kind      , only : rkind
        use phy_model_eq_class   , only : phy_model_eq

        implicit none

        private
        public :: nf90_operators_wr


        !> @class nf90_operators_wr
        !> abstract class encapsulating subroutines to
        !> handle errors when manipulating netcdf files
        !>
        !> @param initialize
        !> initialize the counter to identify the netcdf file
        !> written
        !>
        !> @param write_data
        !> write the data encapsulated in field-type objects
        !> on netcdf output files
        !---------------------------------------------------------------
        type :: nf90_operators_wr

          integer :: nb_timesteps_written

          contains

          procedure,   pass, non_overridable :: initialize
          procedure,   pass, non_overridable :: write_data

        end type nf90_operators_wr


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to initialize the counter identifying
        !> the netcdf file written
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the counter
        !
        !>@param counter
        !> optional argument setting the counter
        !--------------------------------------------------------------
        subroutine initialize(this, counter)

          implicit none

          class(nf90_operators_wr), intent(inout) :: this
          integer       , optional, intent(in)    :: counter

          if(present(counter)) then
             this%nb_timesteps_written=counter
          else
             this%nb_timesteps_written=0
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
        !> 14_08_2013 - initial version - J.L. Desmarais
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
        subroutine write_data(this,f_used,p_model,time)

          implicit none

          class(nf90_operators_wr), intent(inout) :: this
          class(field)            , intent(in)    :: f_used
          class(phy_model_eq)     , intent(in)    :: p_model
          real(rkind)             , intent(in)    :: time

          integer                :: ncid
          integer                :: retval
          integer, dimension(3)  :: coord_id
          integer, dimension(ne) :: data_id
          character(len=16)      :: filename


          !<get the name for the netcdf file
          call get_filename(filename, this%nb_timesteps_written)


          !<create the netcdf file
          retval = NF90_CREATE(trim(filename), NF90_NETCDF4, ncid)
          call nf90_handle_err(retval)


          !<write the header of the file
          call nf90_write_header(ncid,p_model)


          !<define the variables saved in the file
          call nf90_def_var_model(ncid,f_used,p_model,coord_id,data_id)


          !<put the variables in the file
          call nf90_put_var_model(ncid, coord_id, data_id, time, f_used)


          !<close the netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)

          !<increment the internal counter
          this%nb_timesteps_written=this%nb_timesteps_written+1

        end subroutine write_data

      end module nf90_operators_wr_class
