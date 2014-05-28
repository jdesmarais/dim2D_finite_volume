      module bf_layer_errors_module

        use ISO_FORTRAN_ENV, only : ERROR_UNIT

        implicit none

        private
        public :: error_mainlayer_id,
     $            error_diff_mainlayer_id

        contains


        !< error printed if mainlayer_id.ne.(N,S,E,W)
        subroutine error_mainlayer_id(
     $     file_name, fct_name, mainlayer_id)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_id

          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'mainlayer_id not recognized'
          write(ERROR_UNIT, '(''mainlayer_id: '',I2)') mainlayer_id

        end subroutine error_mainlayer_id


        !< error printed if mainlayer_id do not match
        subroutine error_diff_mainlayer_id(
     $     file_name, fct_name,
     $     mainlayer_id1, mainlayer_id2)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_id1
          integer     , intent(in) :: mainlayer_id2

          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'mainlayer_id do not match'
          write(ERROR_UNIT, '(''mainlayer_id1: '',I2)') mainlayer_id1
          write(ERROR_UNIT, '(''mainlayer_id2: '',I2)') mainlayer_id2

        end subroutine error_diff_mainlayer_id

      end module bf_layer_errors_module
