      !> @file
      !> module encapsulating subroutines useful
      !> for i/o interface
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines useful
      !> for i/o interface
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module io_operators_module

        implicit none

        private
        public :: get_filename


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine writing the name of the file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param filename
        !> character(len=16) for the name of the file
        !
        !>@param timestep
        !> integer identifying the number of timesteps
        !> written before writing this file
        !--------------------------------------------------------------
        subroutine get_filename(filename, timestep)

          implicit none

          integer          , intent(in)  :: timestep
          character(len=16), intent(out) :: filename

          character(len=16) :: format_string
          integer           :: size_timestep

          if (timestep.eq.0) then
             size_timestep  = 1
          else
             size_timestep  = floor(log(float(timestep))/log(10.))+1
          end if

         !create the format string
         write (format_string, "(A5,I1,A4)")
     $         '(A4,I', size_timestep, ',A3)'

         !create the name of the file
         write (filename, trim(format_string)) 'data', timestep, '.nc'

        end subroutine get_filename

      end module io_operators_module
