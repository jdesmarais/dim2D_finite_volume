      !---------------------------------------------------------------------------
      !
      ! MODULE: cmd_operators_class
      !
      !> @author 
      !> Julien L. Desmarais
      !
      ! DESCRIPTION:
      !  class extracting information from command line arguments
      !> @brief
      !> class extracting information from command line arguments
      !> to determine the name of netcdf file analyzed as well as
      !> the name of the variable extracted
      !
      ! REVISION HISTORY:
      !>@date
      !> 15_01_2015 - initial version               - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_operators_extract_var_class

        implicit none


        private
        public :: cmd_operators_extract_var


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> class extracting information from command line arguments
        !> @brief
        !> class extracting information from command line arguments
        !> to determine the name of netcdf file analyzed as well as
        !> the name of the variable extracted
        !
        ! REVISION HISTORY:
        ! 15_01_2015 - initial version - J.L. Desmarais
        !
        !>@param : filename_nc               : character(len=1024) giving the path for the
        !                                       netcdf file
        !>@param : filename_out              : character(len=1024) giving the path for the
        !                                       output file
        !>@param : var_name                  : character(len=1024) giving the name of the 
        !                                       variable extracted from the netcdf file
        !
        !>@param : display_help               : display the help for the command line
        !>@param : analyse_cmd_line_arg       : analyse the arguments given in the command line
        !>@param : is_an_option               : determine if the argument is an option
        !>@param : activate_restart           : check if the correct arguments are supplied for the restart option
        !>@param : get_input_filename         : check if the correct arguments are supplied for the input option
        !---------------------------------------------------------------------------
        type :: cmd_operators_extract_var

          character(len=1024) :: filename_nc
          character(len=1024) :: filename_out
          character(len=1024) :: var_name

          logical :: inputs_provided

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: are_inputs_provided
          procedure,   pass :: get_filename_nc
          procedure,   pass :: get_filename_out
          procedure,   pass :: get_var_name

        end type cmd_operators_extract_var

        
        contains


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> display help for the command line
        !> @brief
        !> display the different options available for the command line
        !
        ! REVISION HISTORY:
        ! 15_01_2015 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        subroutine display_help()

          implicit none

          print '(''-h(--help) : display this help'')'
          print '(''-i file.nc : netcdf file analyzed'')'
          print '(''-o file.out: file with variable extracted'')'
          print '(''-v var_name: name of the variable extracted'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators_extract_var)
        !
        ! REVISION HISTORY:
        ! 15_01_2015 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators_extract_var)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators_extract_var), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          logical :: file_nc_provided
          logical :: file_out_provided
          logical :: var_name_provided   
          

          file_nc_provided  = .false.
          file_out_provided = .false.
          var_name_provided = .false.


          !get the total number of arguments given
          !to the executable file
          !---------------------------------------
          arg_nb = command_argument_count()
          
          
          !analyse the arguments given
          !---------------------------
          arg_i=1
          
          do while (arg_i.le.arg_nb)
          

             !get the argument
             call get_command_argument(arg_i, value=cmd_argument)
             

             !check if the option is understood by the program
             select case(trim(cmd_argument))
          
               !input netcdf file
               !=========================================================
               case('-i')
          
                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_nc,
     $                 file_nc_provided)


               !output file
               !=========================================================
               case('-o')

                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_out,
     $                 file_out_provided,
     $                 inquire_file=.false.)


               !name of the variable to be extracted
               !=========================================================
               case('-v')

                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%var_name,
     $                 var_name_provided,
     $                 inquire_file=.false.)


               !help option
               !=========================================================
               case('-h','--help')
          
                  if(arg_i.eq.1) call this%display_help()
          
          
               !option not recognized
               !=========================================================
               case default
                  print '(''****option '', A2, '' ignored****'')',
     $                 trim(cmd_argument)
          
             end select
          
             arg_i=arg_i+1
          
          end do

          this%inputs_provided = 
     $         file_nc_provided.and.
     $         file_out_provided.and.
     $         var_name_provided

        end subroutine analyse_cmd_line_arg


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> check if the argument given is an option
        !> @brief
        !> check if the argument given is an option by testing the first letter
        !> of the argument : if it is a '-', it is an option
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[in] cmd_argument : character*(*)
        !---------------------------------------------------------------------------  
        logical function is_an_option(cmd_argument)

          implicit none

          character*(*), intent(in) :: cmd_argument

          is_an_option = (cmd_argument(1:1)).eq.('-')

        end function is_an_option


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> check if the correct arguments are supplied for the restart option
        !> @brief
        !> check if the correct arguments are supplied for the restart option:
        !> if the file provided for restart, does it exist, is it a netcdf file
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[in]    this  : class(cmd_operators_extract_var)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine get_filename(
     $     arg_i,
     $     arg_nb,
     $     filename,
     $     file_compatible,
     $     inquire_file)

          implicit none

          !i/o variables
          integer                      , intent(inout) :: arg_i
          integer                      , intent(in)    :: arg_nb
          character(len=1024)          , intent(out)   :: filename
          logical                      , intent(out)   :: file_compatible
          logical            , optional, intent(in)    :: inquire_file

          !local variables
          logical :: file_provided
          logical :: netcdf_format
          logical :: file_exists
          logical :: inquire_file_op

          integer             :: lenFilename

          file_provided = .true.
          file_exists   = .true.
          netcdf_format = .true.

          if(present(inquire_file)) then
             inquire_file_op = inquire_file
          else
             inquire_file_op = .true.
          end if


          !increment the argument index to check if
          !it is a potential restart file or not
          arg_i=arg_i+1


          !check if there is a file provided
          if(arg_i.gt.arg_nb) then

             file_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=filename,
     $            LENGTH=lenFilename)
             
             
             !check if the next argument is an option
             if(is_an_option(trim(filename))) then
                
                file_provided = .false.
                arg_i=arg_i-1
                
             else
                
                !check if the file provided exists
                if(inquire_file_op) then
                   INQUIRE(FILE=trim(filename), EXIST=file_exists)

                   if(file_exists) then
                   
                      !check if the file provided uses netcdf format
                      if(filename(lenFilename-2:lenFilename).ne.'.nc') then
                         netcdf_format=.false.
                      end if
                   
                   end if

                else
                   file_exists   = .true.
                   netcdf_format = .true.
                end if                

             end if


             !print if the option is activated or not
             if(.not.file_provided) then
                print '(''file not provided'')'
             end if
             if(.not.netcdf_format) then
                print *, filename, ': no netcdf format'
             end if
             if(.not.file_exists) then
                print *, filename, ':file does not exist'
             end if
             
          end if

          file_compatible =
     $         file_provided.and.
     $         netcdf_format.and.
     $         file_exists

        end subroutine get_filename


        function are_inputs_provided(this)
     $     result(inputs_are_provided)

          implicit none

          class(cmd_operators_extract_var), intent(in) :: this
          logical                                :: inputs_are_provided

          inputs_are_provided = this%inputs_provided

        end function are_inputs_provided


        subroutine get_filename_nc(this,filename)

          implicit none

          class(cmd_operators_extract_var), intent(in) :: this
          character(len=1024)                          :: filename

          filename = this%filename_nc

        end subroutine get_filename_nc


        subroutine get_filename_out(this,filename)

          implicit none

          class(cmd_operators_extract_var), intent(in) :: this
          character(len=1024)                          :: filename

          filename = this%filename_out

        end subroutine get_filename_out


        subroutine get_var_name(this,var_name)

          implicit none

          class(cmd_operators_extract_var), intent(in) :: this
          character(len=1024)                          :: var_name

          var_name = this%var_name

        end subroutine get_var_name

      end module cmd_operators_extract_var_class
      
