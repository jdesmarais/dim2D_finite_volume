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
      module cmd_operators_merge_class

        implicit none


        private
        public :: cmd_operators_merge


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
        type :: cmd_operators_merge

          integer :: nb_tiles_x
          integer :: nb_tiles_y
          integer :: ne
          integer :: bc_size
          integer :: nb_timesteps

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: get_nb_tiles_x
          procedure,   pass :: get_nb_tiles_y
          procedure,   pass :: get_ne
          procedure,   pass :: get_bc_size
          procedure,   pass :: get_nb_timesteps

        end type cmd_operators_merge

        
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

          print '(''-h(--help)        : display this help'')'
          print '(''-x(--nb_tiles_x=) : nb tiles along x-direction'')'
          print '(''-y(--nb_tiles_y=) : nb tiles along y-direction'')'
          print '(''-e(--ne=)         : nb governing variables'')'
          print '(''-b(--bc_size=)    : size of boundary layer'')'
          print '(''-t(--nt=)         : number of timesteps analyzed'')'

        end subroutine display_help


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> analyse the arguments given in the command line
        !> @brief
        !> analyse the arguments given in the command line and initialize the
        !> attributes of class(cmd_operators_merge)
        !
        ! REVISION HISTORY:
        ! 12_01_2015 - initial version - J.L. Desmarais
        !
        !> @param[inout] this : class(cmd_operators_merge)
        !---------------------------------------------------------------------------  
        subroutine analyse_cmd_line_arg(this)


          implicit none


          !i/o variables
          class(cmd_operators_merge), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          
          this%nb_tiles_x   = 1
          this%nb_tiles_y   = 1
          this%ne           = 4
          this%bc_size      = 2
          this%nb_timesteps = 1


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
          
               !number tiles along x-direction
               !=========================================================
               case('-x','--nb_tiles_x')
          
                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%nb_tiles_x)


               !number tiles along y-direction
               !=========================================================
               case('-y','--nb_tiles_y')

                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%nb_tiles_y)


               !number of equations
               !=========================================================
               case('-e','--ne')

                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%ne)


               !size of boundary layer
               !=========================================================
               case('-b','--bc_size')

                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%bc_size)


               !number of timesteps merged
               !=========================================================
               case('-t','--nt')

                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%nb_timesteps)

                  
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
        !> check if the correct arguments are supplied for the restart option
        !> @brief
        !> check if the correct arguments are supplied for the restart option:
        !> if the file provided for restart, does it exist, is it a netcdf file
        !
        ! REVISION HISTORY:
        ! 12_04_2013 - initial version - J.L. Desmarais
        !
        !> @param[in]    this  : class(cmd_operators_merge)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine get_integer_arg(
     $     arg_i,
     $     arg_nb,
     $     arg)

          implicit none

          !i/o variables
          integer, intent(inout) :: arg_i
          integer, intent(in)    :: arg_nb
          integer, intent(out)   :: arg

          !local variables
          character(len=10) :: arg_str
          logical           :: arg_provided

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
     $            VALUE=arg_str)
             
             !check if the next argument is an option
             if(is_an_option(trim(arg_str))) then
                
                arg_provided = .false.
                arg_i=arg_i-1
                
             else
                
                read(arg_str,'(i10)') arg

             end if
             
          end if

          if(.not.arg_provided) then
             print '(''missing option argument'')'
             print '()'
             call display_help()
             print '()'
             stop ''
          end if

        end subroutine get_integer_arg


        function get_nb_tiles_x(this)

          implicit none

          class(cmd_operators_merge), intent(in) :: this
          integer                                :: get_nb_tiles_x

          get_nb_tiles_x = this%nb_tiles_x

        end function get_nb_tiles_x


        function get_nb_tiles_y(this)

          implicit none

          class(cmd_operators_merge), intent(in) :: this
          integer                                :: get_nb_tiles_y

          get_nb_tiles_y = this%nb_tiles_y

        end function get_nb_tiles_y


        function get_ne(this)

          implicit none

          class(cmd_operators_merge), intent(in) :: this
          integer                                :: get_ne

          get_ne = this%ne

        end function get_ne


        function get_bc_size(this)

          implicit none

          class(cmd_operators_merge), intent(in) :: this
          integer                                :: get_bc_size

          get_bc_size = this%bc_size

        end function get_bc_size



        function get_nb_timesteps(this)

          implicit none

          class(cmd_operators_merge), intent(in) :: this
          integer                                :: get_nb_timesteps

          get_nb_timesteps = this%nb_timesteps

        end function get_nb_timesteps

      end module cmd_operators_merge_class
      
