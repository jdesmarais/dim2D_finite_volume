      !---------------------------------------------------------------------------
      !
      ! MODULE: cmd_verify_symmetry_class
      !
      !> @author 
      !> Julien L. Desmarais
      !
      ! DESCRIPTION:
      !  class extracting information from command line arguments
      !> @brief
      !> class extracting information from command line arguments
      !> to determine the restart option, the input file name and
      !> the output root for the file name
      !
      ! REVISION HISTORY:
      !>@date
      !> 12_04_2013 - initial version               - J.L. Desmarais
      !> 05_12_2014 - adjustments for augeanstables - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_verify_symmetry_class

        use parameters_bf_layer, only :
     $     N,S,E,W

        implicit none


        private
        public :: cmd_verify_symmetry


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> class extracting information from command line arguments
        !> @brief
        !> class extracting information from command line arguments
        !> to determine the restart option, the input file name and
        !> the output root for the file name
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !>@param : restart_activated          : logical determining in the restart option is activated or not
        !>@param : restart_filename           : character(len=1024) giving the path for the file use for restart
        !
        !>@param : display_help               : display the help for the command line
        !>@param : analyze_cmd_line_arg       : analyze the arguments given in the command line
        !>@param : is_an_option               : determine if the argument is an option
        !>@param : activate_restart           : check if the correct arguments are supplied for the restart option
        !>@param : get_input_filename         : check if the correct arguments are supplied for the input option
        !---------------------------------------------------------------------------
        type :: cmd_verify_symmetry

          logical               :: analyze_symmetry_activated
          character(len=1024)   :: symmetry_filename1
          character(len=1024)   :: symmetry_filename2

          logical               :: analyze_two_files
          logical               :: analyze_bf_file

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyze_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: activate_symmetry_analyze

          procedure,   pass :: is_analyze_activated
          procedure,   pass :: get_symmetry_filename1
          procedure,   pass :: get_symmetry_filename2
          procedure,   pass :: get_analyze_two_files
          procedure,   pass :: get_analyze_bf_file

        end type cmd_verify_symmetry

        
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
        ! 12_04_2013 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        subroutine display_help()

          implicit none

          print '(''-h(--help)          : display this help'')'
          print '(''-i(--input=) file.nc: file analyzed'')'
          print '(''-bf(--buffer_layer) : analyze a buffer layer file'')'
          print '(''-s(--second_input=) file.nc: second file analyzed'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyze the arguments given in the command line
        !> @brief
        !> analyze the arguments given in the command line and initialize the
        !> attributes of class(cmd_verify_symmetry)
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_verify_symmetry)
        !---------------------------------------------------------------------------  
        subroutine analyze_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_verify_symmetry), intent(inout) :: this


          !local variables
          integer :: arg_i
          integer :: arg_nb
          character(len=1024) :: cmd_argument          
          

          !set the default options for the input
          !file
          !--------------------------------------
          this%analyze_symmetry_activated = .false.
          this%analyze_two_files          = .false.
          this%analyze_bf_file            = .false.


          !get the total number of arguments given
          !to the executable file
          !---------------------------------------
          arg_nb = command_argument_count()
          
          
          !analyze the arguments given
          !---------------------------
          arg_i=1
          
          do while (arg_i.le.arg_nb)
          
             !get the argument
             call get_command_argument(arg_i, value=cmd_argument)
          
             !check if the option is understood by the program
             select case(trim(cmd_argument))
          
               !input file
               !=========================================================
               case('-i','--input=')
          
                  call this%activate_symmetry_analyze(
     $                 arg_i,
     $                 arg_nb,
     $                 this%symmetry_filename1)

               !buffer layer file
               !=========================================================
               case('-bf','--buffer_layer')

                  this%analyze_bf_file = .true.

               !second input file
               !=========================================================
               case('-s','--second_input=')

                  call this%activate_symmetry_analyze(
     $                 arg_i,
     $                 arg_nb,
     $                 this%symmetry_filename2)

                  this%analyze_two_files = .true.
               
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

        end subroutine analyze_cmd_line_arg


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
        !> @param[in]    this  : class(cmd_verify_symmetry)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine activate_symmetry_analyze(
     $     this,
     $     arg_i, arg_nb,
     $     file_analyzed)


          implicit none


          !i/o variables
          class(cmd_verify_symmetry), intent(inout) :: this
          integer                   , intent(inout) :: arg_i
          integer                   , intent(in)    :: arg_nb
          character*(*)             , intent(out)   :: file_analyzed


          !local variables
          logical :: file_provided
          logical :: netcdf_format
          logical :: file_exists

          integer :: flength


          file_provided = .true.
          file_exists   = .true.
          netcdf_format = .true.


          !increment the argument index to check if
          !it is a potential restart file or not
          arg_i=arg_i+1


          !check if there is a file provided
          if(arg_i.gt.arg_nb) then

             file_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=file_analyzed,
     $            LENGTH=flength)
             
             
             !check if the next argument is an option
             if(is_an_option(trim(file_analyzed))) then
                
                file_provided = .false.
                arg_i=arg_i-1
                
             else
                
                !check if the file provided exists
                INQUIRE(FILE=trim(file_analyzed), EXIST=file_exists)
                if(file_exists) then
                   
                   !check if the file provided uses netcdf format
                   if(file_analyzed(flength-2:flength).ne.'.nc') then
                      netcdf_format=.false.
                   end if
                   
                end if

             end if


             !check if restart is activated or not
             this%analyze_symmetry_activated = file_provided.and.netcdf_format.and.file_exists

             
             !print if the option is activated or not
             if(.not.file_provided) then
                print '(''-i: file not provided'')'
             end if
             if(.not.netcdf_format) then
                print '(''-i: no netcdf format'')'
             end if
             if(.not.file_exists) then
                print '(''-i: file does not exist'')'
             end if
             if(this%analyze_symmetry_activated) then
                print '(''-i: analyze activated'')'
             end if
             
          end if

        end subroutine activate_symmetry_analyze


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the argument provided for a specific option
        !> @brief
        !> get the argument provided for a specific option
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[inout] arg_i               : integer indexing the command line arguments tested
        !> @param[inout] arg_nb              : integer giving the total number of arguments
        !> @param[inout] option_arg_provided : logical indicating if an argument is passed to the option
        !> @param[inout] option_arg          : character*(*) saving the argument passed to the option
        !---------------------------------------------------------------------------  
        subroutine get_option_arg(
     $     arg_i, arg_nb,
     $     option_arg_provided, option_arg)

          implicit none


          !i/o variables
          integer            , intent(inout):: arg_i
          integer            , intent(in)   :: arg_nb
          logical            , intent(out)  :: option_arg_provided
          character(len=1024), intent(out)  :: option_arg

          
          option_arg_provided    = .true.


          !increment the argument index to check if
          !it is a potential restart file or not
          arg_i=arg_i+1


          !check if there is a file provided
          if(arg_i.gt.arg_nb) then

             option_arg_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=option_arg)
             
             
             !check if the next argument is an option
             if(is_an_option(trim(option_arg))) then
                
                option_arg_provided = .false.
                arg_i=arg_i-1

             end if

          end if

        end subroutine get_option_arg


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> check whether the restart function is activated
        !> @brief
        !> check whether the restart function is activated
        !
        ! REVISION HISTORY:
        ! 05_12_2013 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        function is_analyze_activated(this)
     $     result(restart_activated)

          implicit none

          class(cmd_verify_symmetry), intent(in) :: this
          logical                          :: restart_activated

          restart_activated = this%analyze_symmetry_activated

        end function is_analyze_activated


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the filename for the restart
        !> @brief
        !> get the filename for the restart
        !
        ! REVISION HISTORY:
        ! 05_12_2013 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        function get_symmetry_filename1(this) result(restart_filename)

          implicit none

          class(cmd_verify_symmetry), intent(in) :: this
          character(len=1024)              :: restart_filename

          restart_filename = this%symmetry_filename1

        end function get_symmetry_filename1


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the filename for the restart
        !> @brief
        !> get the filename for the restart
        !
        ! REVISION HISTORY:
        ! 05_12_2013 - initial version - J.L. Desmarais
        !---------------------------------------------------------------------------
        function get_symmetry_filename2(this) result(restart_filename)

          implicit none

          class(cmd_verify_symmetry), intent(in) :: this
          character(len=1024)              :: restart_filename

          restart_filename = this%symmetry_filename2

        end function get_symmetry_filename2


        function get_analyze_two_files(this) result(analyze_two_files)

          implicit none

          class(cmd_verify_symmetry), intent(in) :: this
          logical                                :: analyze_two_files

          analyze_two_files = this%analyze_two_files

        end function get_analyze_two_files


        function get_analyze_bf_file(this) result(analyze_bf_file)

          implicit none

          class(cmd_verify_symmetry), intent(in) :: this
          logical                                :: analyze_bf_file

          analyze_bf_file = this%analyze_bf_file

        end function get_analyze_bf_file

      end module cmd_verify_symmetry_class
      
