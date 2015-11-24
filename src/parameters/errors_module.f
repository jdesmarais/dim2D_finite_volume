      !> @file
      !> module implementing the subroutines run
      !> when an exception is caught in the program
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the subroutines run
      !> when an exception is caught in the program
      !
      !> @date
      !> 09_06_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module errors_module

        use ISO_FORTRAN_ENV, only :
     $     ERROR_UNIT

        implicit none

        private
        public :: 
     $       error_bc_section_type,
     $       error_bc_choice


        integer, parameter :: ERROR_BC_SECTION_ID_CODE = 0 !<@brief error code for bc_section_id exception
        integer, parameter :: ERROR_BC_CHOICE_ID_CODE  = 1 !<@brief error code for bc_choice_id exception


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if the type of bc_section is
        !> not recognized
        !
        !> @date
        !> 09_06_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param bc_section
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_bc_section_type(
     $     file_name, fct_name, bc_section)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: bc_section

          write(ERROR_UNIT, '(I3)') ERROR_BC_SECTION_ID_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'bc_section not recognized'
          write(ERROR_UNIT, '(''bc_section: '',I2)') bc_section

          stop 'error_bc_section'

        end subroutine error_bc_section_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if bc_choice not recognized
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param bc_choice
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_bc_choice(
     $     file_name, fct_name, bc_choice)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: bc_choice

          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'bc_choice not recognized'
          write(ERROR_UNIT, '(''bc_choice: '',I2)') bc_choice

          stop 'error_bc_choice'

        end subroutine error_bc_choice

      end module errors_module
