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
        public ::
     $       get_filename,
     $       get_timestep


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
        subroutine get_filename(filename, timestep, rank)

          implicit none

          integer          , intent(in)  :: timestep
          character(len=16), intent(out) :: filename
          integer, optional, intent(in)  :: rank

          integer           :: size_timestep
          integer           :: size_rank
          character(len=16) :: format_string


          ! determine the character size to write the timestep
          if (timestep.eq.0) then
             size_timestep  = 1
          else
             size_timestep  = floor(log(float(timestep))/log(10.))+1
          end if

          if(present(rank)) then

             ! determine the character size to write the timestep
             if (rank.eq.0) then
                size_rank  = 1
             else
                size_rank  = floor(log(float(rank))/log(10.))+1
             end if

             !create the format string
             write (format_string, "(A5,I1,A5,I1,A4)")
     $            '(A4,I', size_timestep,
     $            ',A1,I', size_rank,
     $            ',A3)'

             !create the name of the file
             write (filename, trim(format_string)) 'data', timestep, '_', rank, '.nc'

          else

             !create the format string
             write (format_string, "(A5,I1,A4)")
     $            '(A4,I', size_timestep, ',A3)'

             !create the name of the file
             write (filename, trim(format_string)) 'data', timestep, '.nc'

          end if

        end subroutine get_filename


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the timestep corresponding to the filename
        !
        !> @date
        !> 05_12_2014 - initial version - J.L. Desmarais
        !
        !>@param filename
        !> character(len=16) for the name of the file
        !
        !>@return timestep
        !> integer identifying the number of timesteps
        !> written before writing this file
        !--------------------------------------------------------------
        function get_timestep(filename) result(timestep)

          implicit none

          character*(*), intent(in) :: filename
          integer                   :: timestep

          integer :: length

          !the name of the input file has the structure:
          !'data'+timestep+'.nc'
          length   = len(filename)
          read(filename(5:length-3),'(i10)') timestep

        end function get_timestep

      end module io_operators_module
