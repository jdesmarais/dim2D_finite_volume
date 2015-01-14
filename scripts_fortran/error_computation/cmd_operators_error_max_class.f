      !---------------------------------------------------------------------------
      !
      ! MODULE: cmd_operators_error_max_class
      !
      !> @author 
      !> Julien L. Desmarais
      !
      ! DESCRIPTION:
      !  class extracting information from command line arguments
      !> @brief
      !> class extracting information from command line arguments
      !> to determine the directory for the error files and the
      !> name for the output file
      !
      ! REVISION HISTORY:
      !>@date
      !> 14_01_2015 - initial version - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_operators_error_max_class

        implicit none


        private
        public :: cmd_operators_error_max


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> class extracting information from command line arguments
        !> @brief
        !> class extracting information from command line arguments
        !> to determine the directory for the error files and the
        !> name for the output file
        !
        ! REVISION HISTORY:
        ! 14_01_2015 - initial version - J.L. Desmarais
        !
        !>@param : directory_error_files      : character(len=1024) giving the path for the
        !                                       directory with the error files analyzed
        !>@param : filename_error_max         : character(len=1024) giving the path for the
        !                                       simulation file on error domain
        !>@param : nb_files                   : total number of files analyzed
        !
        !>@param : display_help               : display the help for the command line
        !>@param : analyse_cmd_line_arg       : analyse the arguments given in the command line
        !>@param : is_an_option               : determine if the argument is an option
        !>@param : activate_restart           : check if the correct arguments are supplied for the restart option
        !>@param : get_input_filename         : check if the correct arguments are supplied for the input option
        !---------------------------------------------------------------------------
        type :: cmd_operators_error_max

          character(len=1024) :: directory_error_files
          character(len=1024) :: filename_error_max
          integer             :: nb_files
          
          logical :: inputs_provided

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: are_inputs_provided
          procedure,   pass :: get_directory_error_files
          procedure,   pass :: get_filename_error_max
          procedure,   pass :: get_nb_files

        end type cmd_operators_error_max

        
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

          print '(''-h(--help)   : display this help'')'
          print '(''-i <dir>     : directoy with the error files'')'
          print '(''-o <file.nc> : file for error output'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators_error_max)
        !
        ! REVISION HISTORY:
        ! 12_01_2015 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators_error_max)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators_error_max), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          logical :: directory_provided
          logical :: file_error_provided
          logical :: nb_files_provided
          

          directory_provided  = .false.
          file_error_provided = .false.


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
          
               !directory for the error files analyzed
               !=========================================================
               case('-i')
          
                  call get_path(
     $                 arg_i,
     $                 arg_nb,
     $                 this%directory_error_files,
     $                 directory_provided,
     $                 dir_check=.true.)


               !filename for the output file with the error max
               !=========================================================
               case('-o')

                  call get_path(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename_error_max,
     $                 file_error_provided,
     $                 inquire_path=.false.)


               !number of files analyzed for the error max
               !=========================================================
               case('-n')

                  call analyze_nb_files(
     $                 arg_i,
     $                 arg_nb,
     $                 this%nb_files,
     $                 nb_files_provided)


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
     $         directory_provided.and.
     $         file_error_provided.and.
     $         nb_files_provided

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
        !> check if the correct arguments are supplied to the program
        !> @brief
        !> check if the path provided exists if it is a directory for the
        !> output file, check if the outfile file is of netcdf format
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[in]    this  : class(cmd_operators_error_max)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine get_path(
     $     arg_i,
     $     arg_nb,
     $     path,
     $     path_compatible,
     $     dir_check,
     $     inquire_path)

          implicit none

          !i/o variables
          integer                      , intent(inout) :: arg_i
          integer                      , intent(in)    :: arg_nb
          character(len=1024)          , intent(out)   :: path
          logical                      , intent(out)   :: path_compatible
          logical            , optional, intent(in)    :: dir_check
          logical            , optional, intent(in)    :: inquire_path

          !local variables
          logical :: path_provided
          logical :: netcdf_format
          logical :: path_exists
          logical :: inquire_path_op
          logical :: dir_check_op

          integer :: lenPath

          path_provided = .true.
          path_exists   = .true.
          netcdf_format = .true.

          if(present(inquire_path)) then
             inquire_path_op = inquire_path
          else
             inquire_path_op = .true.
          end if

          if(present(dir_check)) then
             dir_check_op = dir_check
          else
             dir_check_op = .false.
          end if


          !increment the argument index to check if
          !it is a potential restart path or not
          arg_i=arg_i+1


          !check if there is a path provided
          if(arg_i.gt.arg_nb) then

             path_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=path,
     $            LENGTH=lenPath)
             
             
             !check if the next argument is an option
             if(is_an_option(trim(path))) then
                
                path_provided = .false.
                arg_i=arg_i-1
                
             else
                
                ! check if the path provided exists
                if(inquire_path_op) then

                   ! check if the directory provided exists
                   if(dir_check_op) then

                      INQUIRE (DIRECTORY=trim(path), EXIST=path_exists)

                   ! check if the path provided exists
                   else

                      INQUIRE(FILE=trim(path), EXIST=path_exists)

                   end if

                else

                   path_exists = .true.

                end if


                if(path_exists) then
                   
                   ! if this is not a directory which is checked
                   ! check if the path provided uses netcdf format
                   if(.not.dir_check_op) then

                      if(path(lenPath-2:lenPath).ne.'.nc') then
                         netcdf_format=.false.
                      end if

                   end if

                end if

             end if


             !print if the option is activated or not
             if(.not.path_provided) then
                print '(''path not provided'')'
             end if
             if(.not.path_exists) then
                print *, path, ':path does not exist'
             end if
             if(.not.netcdf_format) then
                print *, path, ': no netcdf format'
             end if
             
             
          end if

          path_compatible =
     $         path_provided.and.
     $         path_exists.and.
     $         netcdf_format

        end subroutine get_path


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> check if the correct arguments are supplied to the program
        !> @brief
        !> check if the path provided exists if it is a directory for the
        !> output file, check if the outfile file is of netcdf format
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[in]    this  : class(cmd_operators_error_max)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine analyze_nb_files(
     $     arg_i,
     $     arg_nb,
     $     nb_files,
     $     nb_files_provided)

          implicit none

          !i/o variables
          integer, intent(inout) :: arg_i
          integer, intent(in)    :: arg_nb
          integer, intent(out)   :: nb_files
          logical, intent(out)   :: nb_files_provided


          !local variables
          character(len=1024) :: nb_files_str
          integer             :: len_nb_files_str
          integer             :: integer_size
          character(len=4)    :: format_string

          ! increment the argument index to check if
          ! the argument after cooresponds to a
          ! potential number of files
          arg_i=arg_i+1


          ! check if there is a path provided
          if(arg_i.gt.arg_nb) then

             nb_files_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=nb_files_str,
     $            LENGTH=len_nb_files_str)             
             
             ! check if the next argument is an option
             if(is_an_option(trim(nb_files_str))) then
                
                nb_files_provided = .false.
                arg_i=arg_i-1
                
             else
                
                ! size of the integer
                integer_size = len(trim(nb_files_str))

                ! convert the string into integer
                write(format_string, '(A2,I1,A1)')
     $               '(I',integer_size,')'

                ! get the integer
                read(nb_files_str,format_string) nb_files
             
                nb_files_provided = .true.

             end if

          end if


        end subroutine analyze_nb_files


        function are_inputs_provided(this)
     $     result(inputs_are_provided)

          implicit none

          class(cmd_operators_error_max), intent(in) :: this
          logical                                    :: inputs_are_provided

          inputs_are_provided = this%inputs_provided

        end function are_inputs_provided


        subroutine get_directory_error_files(this,directory)

          implicit none

          class(cmd_operators_error_max), intent(in)  :: this
          character(len=1024)           , intent(out) :: directory

          directory = this%directory_error_files

        end subroutine get_directory_error_files


        subroutine get_filename_error_max(this,filename)

          implicit none

          class(cmd_operators_error_max), intent(in)  :: this
          character(len=1024)           , intent(out) :: filename

          filename = this%filename_error_max

        end subroutine get_filename_error_max


        subroutine get_nb_files(this,nb_files)

          implicit none

          class(cmd_operators_error_max), intent(in)  :: this
          integer                       , intent(out) :: nb_files

          nb_files = this%nb_files

        end subroutine get_nb_files

      end module cmd_operators_error_max_class
      
