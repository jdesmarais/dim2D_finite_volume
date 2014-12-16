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
      !> to determine the restart option, the input file name and
      !> the output root for the file name
      !
      ! REVISION HISTORY:
      !>@date
      !> 12_04_2013 - initial version               - J.L. Desmarais
      !> 05_12_2014 - adjustments for augeanstables - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_operators_class

        use parameters_bf_layer, only :
     $     N,S,E,W

        implicit none


        private
        public :: cmd_operators


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
        !>@param : analyse_cmd_line_arg       : analyse the arguments given in the command line
        !>@param : is_an_option               : determine if the argument is an option
        !>@param : activate_restart           : check if the correct arguments are supplied for the restart option
        !>@param : get_input_filename         : check if the correct arguments are supplied for the input option
        !---------------------------------------------------------------------------
        type :: cmd_operators

          logical               :: restart_activated
          character(len=1024)   :: restart_filename

          integer, dimension(4) :: nb_bf_layers

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: activate_restart

          procedure,   pass :: is_restart_activated
          procedure,   pass :: get_restart_filename
          procedure,   pass :: get_nb_bf_layers
          procedure,   pass :: get_nb_layers

        end type cmd_operators

        
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

          print '(''-h(--help)            : display this help'')'
          print '(''-r(--restart=) file.nc: use file to restart'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators)
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators), intent(inout) :: this


          !local variables
          integer :: arg_i
          integer :: arg_nb
          character(len=1024) :: cmd_argument          
          

          !set the default options for the input
          !file
          !--------------------------------------
          this%restart_activated=.false.


          !set the default options for the number
          !of buffer layers
          !--------------------------------------
          this%nb_bf_layers = [0,0,0,0]


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
          
               !restart option
               !=========================================================
               case('-r','--restart')
          
                  call this%activate_restart(
     $                 arg_i,
     $                 arg_nb)

               !number North buffer layers
               !=========================================================
               case('--nb_bf_N')

                  this%nb_bf_layers(N) = this%get_nb_layers(
     $                 arg_i,
     $                 arg_nb)

               !number South buffer layers
               !=========================================================
               case('--nb_bf_S')

                  this%nb_bf_layers(S) = this%get_nb_layers(
     $                 arg_i,
     $                 arg_nb)
               
               !number East buffer layers
               !=========================================================
               case('--nb_bf_E')

                  this%nb_bf_layers(E) = this%get_nb_layers(
     $                 arg_i,
     $                 arg_nb)

               !number West buffer layers
               !=========================================================
               case('--nb_bf_W')

                  this%nb_bf_layers(W) = this%get_nb_layers(
     $                 arg_i,
     $                 arg_nb)
          
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
        !> @param[in]    this  : class(cmd_operators)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine activate_restart(
     $     this,
     $     arg_i, arg_nb)


          implicit none


          !i/o variables
          class(cmd_operators), intent(inout) :: this
          integer             , intent(inout) :: arg_i
          integer             , intent(in)    :: arg_nb


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
     $            VALUE=this%restart_filename,
     $            LENGTH=flength)
             
             
             !check if the next argument is an option
             if(is_an_option(trim(this%restart_filename))) then
                
                file_provided = .false.
                arg_i=arg_i-1
                
             else
                
                !check if the file provided exists
                INQUIRE(FILE=trim(this%restart_filename), EXIST=file_exists)
                if(file_exists) then
                   
                   !check if the file provided uses netcdf format
                   if(this%restart_filename(flength-2:flength).ne.'.nc') then
                      netcdf_format=.false.
                   end if
                   
                end if

             end if


             !check if restart is activated or not
             this%restart_activated = file_provided.and.netcdf_format.and.file_exists

             
             !print if the option is activated or not
             if(.not.file_provided) then
                print '(''-r: file not provided'')'
             end if
             if(.not.netcdf_format) then
                print '(''-r: no netcdf format'')'
             end if
             if(.not.file_exists) then
                print '(''-r: file does not exist'')'
             end if
             if(this%restart_activated) then
                print '(''-r: restart activated'')'
             end if
             
          end if

        end subroutine activate_restart


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
        function is_restart_activated(this)
     $     result(restart_activated)

          implicit none

          class(cmd_operators), intent(in) :: this
          logical                          :: restart_activated

          restart_activated = this%restart_activated

        end function is_restart_activated


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
        function get_restart_filename(this) result(restart_filename)

          implicit none

          class(cmd_operators), intent(in) :: this
          character(len=1024)              :: restart_filename

          restart_filename = this%restart_filename

        end function get_restart_filename


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the number of buffer layers given by the option
        !> @brief
        !> get the number of buffer layers given by the option
        !
        ! REVISION HISTORY:
        ! 16_12_2013 - initial version - J.L. Desmarais
        !
        !> @param[in]    this  : class(cmd_operators)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        function get_nb_layers(
     $     this,
     $     arg_i, arg_nb)
     $     result(nb_bf_layers)


          implicit none


          !i/o variables
          class(cmd_operators), intent(inout) :: this
          integer             , intent(inout) :: arg_i
          integer             , intent(in)    :: arg_nb
          integer                             :: nb_bf_layers

          logical             :: option_arg_provided
          character(len=1024) :: option_arg
          
          call get_option_arg(
     $         arg_i, arg_nb,
     $         option_arg_provided, option_arg)
          
          if(option_arg_provided) then
             read( option_arg, '(i10)' ) nb_bf_layers
          else
             nb_bf_layers = 0
          end if

        end function get_nb_layers


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the number of buffer layers given by the user
        !> @brief
        !> get the number of buffer layers given by the user
        !
        ! REVISION HISTORY:
        ! 16_12_2014 - initial version - J.L. Desmarais
        !
        !> @param[in]  this        : class(cmd_operators)
        !> @param[out] nb_bf_layers: number of buffer layers
        !---------------------------------------------------------------------------  
        function get_nb_bf_layers(
     $     this)
     $     result(nb_bf_layers)


          implicit none


          !i/o variables
          class(cmd_operators), intent(in) :: this
          integer, dimension(4)            :: nb_bf_layers

          nb_bf_layers = this%nb_bf_layers

        end function get_nb_bf_layers

      end module cmd_operators_class
      
