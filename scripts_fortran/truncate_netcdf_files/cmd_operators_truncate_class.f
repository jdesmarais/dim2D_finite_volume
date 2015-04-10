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
      module cmd_operators_truncate_class

        use parameters_kind, only :
     $     rkind

        implicit none        


        private
        public :: cmd_operators_truncate


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
        type :: cmd_operators_truncate

          real(rkind) :: x_min
          real(rkind) :: x_max
          real(rkind) :: y_min
          real(rkind) :: y_max
          integer     :: ne
          integer     :: nb_timesteps

          contains

          procedure, nopass :: display_help
          procedure,   pass :: analyse_cmd_line_arg
          procedure, nopass :: is_an_option
          procedure,   pass :: get_borders
          procedure,   pass :: get_ne
          procedure,   pass :: get_nb_timesteps

        end type cmd_operators_truncate

        
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
          print '(''-xmin             : x_min for truncation'')'
          print '(''-xmax             : x_max for truncation'')'
          print '(''-ymin             : y_min for truncation'')'
          print '(''-ymax             : y_max for truncation'')'
          print '(''-e(--ne=)         : nb governing variables'')'
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
          class(cmd_operators_truncate), intent(inout) :: this


          !local variables
          integer             :: arg_i
          integer             :: arg_nb
          character(len=1024) :: cmd_argument
          
          
          this%ne           = 4
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
          
               !x_min
               !=========================================================
               case('-xmin')
          
                  call get_real_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%x_min)


               !x_min
               !=========================================================
               case('-xmax')
          
                  call get_real_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%x_max)


               !x_min
               !=========================================================
               case('-ymin')
          
                  call get_real_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%y_min)


               !y_max
               !=========================================================
               case('-ymax')
          
                  call get_real_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%y_max)


               !number of equations
               !=========================================================
               case('-e','--ne')

                  call get_integer_arg(
     $                 arg_i,
     $                 arg_nb,
     $                 this%ne)


               !number of timesteps truncated
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
        !> @param[in]    this  : class(cmd_operators_truncate)
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
        !> @param[in]    this  : class(cmd_operators_truncate)
        !> @param[inout] arg_i : integer indexing the command line arguments tested
        !> @param[inout] arg_nb: integer giving the total number of arguments
        !---------------------------------------------------------------------------  
        subroutine get_real_arg(
     $     arg_i,
     $     arg_nb,
     $     arg)

          implicit none

          !i/o variables
          integer    , intent(inout) :: arg_i
          integer    , intent(in)    :: arg_nb
          real(rkind), intent(out)   :: arg

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
             
             read(arg_str,'(F20.13)') arg
             
          end if

          if(.not.arg_provided) then
             print '(''missing option argument'')'
             print '()'
             call display_help()
             print '()'
             stop ''
          end if

        end subroutine get_real_arg


        function get_ne(this)

          implicit none

          class(cmd_operators_truncate), intent(in) :: this
          integer                                :: get_ne

          get_ne = this%ne

        end function get_ne


        function get_borders(this)

          implicit none

          class(cmd_operators_truncate), intent(in) :: this
          real(rkind), dimension(2,2)               :: get_borders

          get_borders(1,1) = this%x_min
          get_borders(1,2) = this%x_max
          get_borders(2,1) = this%y_min
          get_borders(2,2) = this%y_max          

        end function get_borders



        function get_nb_timesteps(this)

          implicit none

          class(cmd_operators_truncate), intent(in) :: this
          integer                                :: get_nb_timesteps

          get_nb_timesteps = this%nb_timesteps

        end function get_nb_timesteps

      end module cmd_operators_truncate_class
      
