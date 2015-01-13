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
      !> to determine the name of the small domain and the large
      !> domain simulation files
      !
      ! REVISION HISTORY:
      !>@date
      !> 12_01_2015 - initial version               - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_operators_error_class

        implicit none


        private
        public :: cmd_operators_error


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> class extracting information from command line arguments
        !> @brief
        !> class extracting information from command line arguments
        !> to determine the name of the small domain and the large
        !> domain simulation files
        !
        ! REVISION HISTORY:
        ! 12_01_2015 - initial version - J.L. Desmarais
        !
        !>@param : filename_sm_domain         : character(len=1024) giving the path for the
        !                                       simulation file on small domain
        !>@param : filename_lg_domain         : character(len=1024) giving the path for the
        !                                       simulation file on large domain
        !>@param : filename_error             : character(len=1024) giving the path for the
        !                                       simulation file on error domain
        !
        !>@param : display_help               : display the help for the command line
        !>@param : analyse_cmd_line_arg       : analyse the arguments given in the command line
        !>@param : is_an_option               : determine if the argument is an option
        !>@param : activate_restart           : check if the correct arguments are supplied for the restart option
        !>@param : get_input_filename         : check if the correct arguments are supplied for the input option
        !---------------------------------------------------------------------------
        type :: cmd_operators_error

          character(len=1024) :: filename_sm_domain
          character(len=1024) :: filename_lg_domain
          character(len=1024) :: filename_error

          logical :: files_provided

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: are_files_provided
          procedure,   pass :: get_filename_sm_domain
          procedure,   pass :: get_filename_lg_domain
          procedure,   pass :: get_filename_error

        end type cmd_operators_error

        
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
        ! 12_01_2015 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        subroutine display_help()

          implicit none

          print '(''-h(--help)              : display this help'')'
          print '(''-s(--sm_domain=) file.nc: file for small domain'')'
          print '(''-l(--lg_domain=) file.nc: file for large domain'')'
          print '(''-o(--output=)    file.nc: file for error output'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators_error)
        !
        ! REVISION HISTORY:
        ! 12_01_2015 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators_error)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators_error), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          logical :: file_sm_domain_provided
          logical :: file_lg_domain_provided
          logical :: file_error_provided   
          

          file_sm_domain_provided = .false.
          file_lg_domain_provided = .false.
          file_error_provided     = .false.


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
          
               !small domain simulation file
               !=========================================================
               case('-s','--small_domain')
          
                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_sm_domain,
     $                 file_sm_domain_provided)


               !large domain simulation file
               !=========================================================
               case('-l','--large_domain')

                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_lg_domain,
     $                 file_lg_domain_provided)


               !output simulation file
               !=========================================================
               case('-o','--output')

                  call get_filename(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_error,
     $                 file_error_provided,
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

          this%files_provided = 
     $         file_sm_domain_provided.and.
     $         file_lg_domain_provided.and.
     $         file_error_provided   

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
        !> @param[in]    this  : class(cmd_operators_error)
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
                else
                   file_exists = .true.
                end if

                if(file_exists) then
                   
                   !check if the file provided uses netcdf format
                   if(filename(lenFilename-2:lenFilename).ne.'.nc') then
                      netcdf_format=.false.
                   end if
                   
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


        function are_files_provided(this)
     $     result(files_are_provided)

          implicit none

          class(cmd_operators_error), intent(in) :: this
          logical                                :: files_are_provided

          files_are_provided = this%files_provided

        end function are_files_provided


        subroutine get_filename_sm_domain(this,filename)

          implicit none

          class(cmd_operators_error), intent(in) :: this
          character(len=1024)                    :: filename

          filename = this%filename_sm_domain

        end subroutine get_filename_sm_domain


        subroutine get_filename_lg_domain(this,filename)

          implicit none

          class(cmd_operators_error), intent(in) :: this
          character(len=1024)                    :: filename

          filename = this%filename_lg_domain

        end subroutine get_filename_lg_domain


        subroutine get_filename_error(this,filename)

          implicit none

          class(cmd_operators_error), intent(in) :: this
          character(len=1024)                    :: filename

          filename = this%filename_error

        end subroutine get_filename_error

      end module cmd_operators_error_class
      
