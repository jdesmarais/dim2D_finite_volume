      program extract_1D_variable

        use cmd_operators_extract_var_class, only :
     $     cmd_operators_extract_var

        use netcdf
        
        use nf90_operators_error_module, only :
     $       nf90_handle_err

        use parameters_kind, only :
     $       rkind


        implicit none


        type(cmd_operators_extract_var) :: cmd_operators_used

        character(len=1024) :: filename_nc
        character(len=1024) :: filename_out
        character(len=1024) :: var_name

        real(rkind), dimension(:), allocatable :: var


        ! analyze the command line arguments to extract:
        ! - filename_nc
        ! - filename_out
        ! - var_name
        call cmd_operators_used%analyse_cmd_line_arg()

        if(.not. cmd_operators_used%are_inputs_provided()) then
           call cmd_operators_used%display_help()
           stop ''
        end if

        call cmd_operators_used%get_filename_nc(filename_nc)
        call cmd_operators_used%get_filename_out(filename_out)
        call cmd_operators_used%get_var_name(var_name)


        ! extract and write the file with the variable extracted
        call extract_var(filename_nc,var_name,var)
        call write_var(var,filename_out)


        contains


        subroutine write_var(var,filename_out)

          implicit none

          real(rkind), dimension(:), intent(in) :: var
          character*(*)            , intent(in) :: filename_out

          integer :: ios


          !1) open the binary file to write the data
          open(unit=2,
     $         file=trim(filename_out),
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)


          !2) write data
          if(ios.eq.0) then
             write(unit=2, iostat=ios) var
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          !3) close file
          close(unit=2)

        end subroutine write_var


        subroutine extract_var(filename,var_name,var)

          implicit none

          character*(*)                           , intent(in)  :: filename
          character*(*)                           , intent(in)  :: var_name
          real(rkind)  , dimension(:), allocatable, intent(out) :: var


          integer                            :: retval
          integer                            :: ncid
          integer                            :: var_id
          integer                            :: ndims
          integer, dimension(:), allocatable :: dimids
          integer                            :: var_size


          !1) open the netcdf file for reading
          retval = NF90_OPEN(
     $         trim(filename),
     $         NF90_NOWRITE,
     $         ncid)
          call nf90_handle_err(retval)


          !2) determine the var_id corresponding to the
          !   variable with the name var_name
          retval = NF90_INQ_VARID(
     $         ncid,
     $         var_name,
     $         var_id)
          call nf90_handle_err(retval)
          

          !3) determine the length of the table 
          retval = NF90_INQUIRE_VARIABLE(
     $         ncid,
     $         var_id,
     $         ndims=ndims)
          call nf90_handle_err(retval)

          if(ndims.ne.1) then
             print '(''extract_1D_variable'')'
             print '(''extract_var'')'
             print '(''the variable is not 1D'')'
             stop ''
          end if

          allocate(dimids(ndims))

          retval = NF90_INQUIRE_VARIABLE(
     $         ncid,
     $         var_id,
     $         dimids=dimids)
          call nf90_handle_err(retval)

          retval = NF90_INQUIRE_DIMENSION(
     $         ncid,
     $         dimids(1),
     $         len=var_size)
          call nf90_handle_err(retval)

          deallocate(dimids)


          !4) extract the variable
          allocate(var(var_size))
          retval = NF90_GET_VAR(
     $         ncid,
     $         var_id,
     $         var)
          call nf90_handle_err(retval)


          !5) close the netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)
          

        end subroutine extract_var

      end program extract_1D_variable
