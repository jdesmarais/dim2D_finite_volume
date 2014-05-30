      module bf_layer_errors_module

        use ISO_FORTRAN_ENV, only : ERROR_UNIT

        implicit none

        private
        public :: error_mainlayer_id,
     $            error_diff_mainlayer_id,
     $            error_incompatible_neighbor,
     $            error_neighbor_index


        integer, parameter :: ERROR_MAINLAYER_ID_CODE          = 0
        integer, parameter :: ERROR_DIFF_MAINLAYER_ID_CODE     = 1
        integer, parameter :: ERROR_INCOMPATIBLE_NEIGHBOR_CODE = 2
        integer, parameter :: ERROR_NEIGHBOR_INDEX_CODE        = 3


        contains


        !< error printed if mainlayer_id.ne.(N,S,E,W)
        subroutine error_mainlayer_id(
     $     file_name, fct_name, mainlayer_id)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_id

          write(ERROR_UNIT, '(I3)') ERROR_MAINLAYER_ID_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'mainlayer_id not recognized'
          write(ERROR_UNIT, '(''mainlayer_id: '',I2)') mainlayer_id

          stop 'error_mainlayer_id'

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

          write(ERROR_UNIT, '(I3)') ERROR_DIFF_MAINLAYER_ID_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'mainlayer_id do not match'
          write(ERROR_UNIT, '(''mainlayer_id1: '',I2)') mainlayer_id1
          write(ERROR_UNIT, '(''mainlayer_id2: '',I2)') mainlayer_id2

          stop 'error_diff_mainlayer_id'

        end subroutine error_diff_mainlayer_id


        !< error printed if the neighbor is incompatible
        !> with the current buffer layer
        subroutine error_incompatible_neighbor(
     $     file_name, fct_name,
     $     mainlayer_id1, mainlayer_id2)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_id1
          integer     , intent(in) :: mainlayer_id2

          write(ERROR_UNIT, '(I3)') ERROR_INCOMPATIBLE_NEIGHBOR_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'incompatible neighbors'
          write(ERROR_UNIT, '(''mainlayer_id1: '',I2)') mainlayer_id1
          write(ERROR_UNIT, '(''mainlayer_id2: '',I2)') mainlayer_id2
          
          stop 'error_incompatible_neighbor'  

        end subroutine error_incompatible_neighbor


        !< error neighbor index.ne.(1,2)
        subroutine error_neighbor_index(
     $     file_name, fct_name, neighbor_index)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: neighbor_index

          write(ERROR_UNIT, '(I3)') ERROR_NEIGHBOR_INDEX_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'neighbor_index not recognized'
          write(ERROR_UNIT, '(''neighbor_index: '',I2)') neighbor_index

          stop 'error_neighbor_index'

        end subroutine error_neighbor_index

      end module bf_layer_errors_module
