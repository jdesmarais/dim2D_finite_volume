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
      !> to determine the filename where the governing variables
      !> (mass, momentum, energy) are saved
      !
      ! REVISION HISTORY:
      !>@date
      !> 11_08_2015 - initial version - J.L. Desmarais
      !---------------------------------------------------------------------------
      module cmd_operators_extract_class

        use parameters_kind, only :
     $     rkind

        implicit none


        private
        public :: cmd_operators_extract


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
        ! 08_04_2015 - initial version - J.L. Desmarais
        !
        !>@param : nb_tiles_x   : number of tiles along x-direction
        !>@param : nb_tiles_y   : number of tiles along y-direction
        !>@param : ne           : number of governing variables
        !>@param : bc_size      : size of the boundary edge
        !>@param : nb_timesteps : number of timesteps
        !---------------------------------------------------------------------------
        type :: cmd_operators_extract

          character(len=100) :: filename

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: get_filename

        end type cmd_operators_extract

        
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

          print '(''-h(--help) : display this help'')'
          print '(''-i         : input filename'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators_truncate)
        !
        ! REVISION HISTORY:
        ! 12_01_2015 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators_truncate)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators_extract), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          
          !get the total number of arguments given
          !to the executable file
          !---------------------------------------
          arg_nb = command_argument_count()
          if(arg_nb.eq.0) call display_help()
          
          !analyse the arguments given
          !---------------------------
          arg_i=1

          do while (arg_i.le.arg_nb)
          

             !get the argument
             call get_command_argument(arg_i, value=cmd_argument)
             

             !check if the option is understood by the program
             select case(trim(cmd_argument))
          
               !number of timesteps truncated
               !=========================================================
               case('-i')

                  call get_character_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%filename)

                  
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
        function is_an_option(cmd_argument)

          implicit none

          character*(*), intent(in) :: cmd_argument
          logical                   :: is_an_option

          is_an_option = (cmd_argument(1:1)).eq.('-')

        end function is_an_option


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> determine a character-like argument
        !> @brief
        !> check if the argument is supplied and extract it
        !> from the command line
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[inout] arg_i  : integer indexing the command line arguments tested
        !> @param[inout] arg_nb : integer giving the total number of arguments
        !> @param[out]   arg    : character argument supplied
        !---------------------------------------------------------------------------  
        subroutine get_character_arg(
     $     arg_i,
     $     arg_nb,
     $     arg)

          implicit none

          !i/o variables
          integer           , intent(inout) :: arg_i
          integer           , intent(in)    :: arg_nb
          character(len=100), intent(out)   :: arg

          !local variables
          character(len=100) :: arg_str
          logical            :: arg_provided

          arg_provided = .true.


          !increment the argument index to check if
          !it is a potential restart file or not
          arg_i=arg_i+1


          !check if there is a file provided
          if(arg_i.gt.arg_nb) then

             arg_provided = .false.
             
          else
                   
             call get_command_argument(
     $            arg_i,
     $            VALUE=arg)
             
             !check if the next argument is an option
             if(is_an_option(trim(arg_str))) then
                
                arg_provided = .false.
                arg_i=arg_i-1
                
             end if
             
          end if

          if(.not.arg_provided) then
             print '(''missing option argument'')'
             print '()'
             call display_help()
             print '()'
             stop ''
          end if

        end subroutine get_character_arg


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> get the filename from the cnd operator
        !> @brief
        !> get the filename from the cnd operator
        !
        ! REVISION HISTORY:
        ! 11_08_2015 - initial version - J.L. Desmarais
        !
        !> @param[in]  this     : cmd operator
        !> @param[out] filename : netcdf file path
        !---------------------------------------------------------------------------  
        function get_filename(this) result(filename)

          implicit none

          class(cmd_operators_extract), intent(in) :: this
          character(len=100)                       :: filename

          filename = this%filename

        end function get_filename

      end module cmd_operators_extract_class
