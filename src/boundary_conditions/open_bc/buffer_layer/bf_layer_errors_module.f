      !> @file
      !> module implementing the subroutines run when an exception is
      !> caught in files implementing the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the subroutines run when an exception is
      !> caught in files implementing the buffer layers
      !
      !> @date
      ! 26_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_errors_module

        use ISO_FORTRAN_ENV, only : ERROR_UNIT

        implicit none

        private
        public ::
     $       error_mainlayer_id,
     $       error_diff_mainlayer_id,
     $       error_incompatible_neighbor,
     $       error_neighbor_index,
     $       error_overlap_index,
     $       error_overlap_incompatible,
     $       error_bc_section_type,
     $       error_gradient_type,
     $       error_cpt_overlap_index,
     $       error_cpt_type,
     $       error_mainlayer_interface_type,
     $       error_mainlayer_interface_incompatible


        integer, parameter :: ERROR_MAINLAYER_ID_CODE                     = 0
        integer, parameter :: ERROR_DIFF_MAINLAYER_ID_CODE                = 1
        integer, parameter :: ERROR_INCOMPATIBLE_NEIGHBOR_CODE            = 2
        integer, parameter :: ERROR_NEIGHBOR_INDEX_CODE                   = 3
        integer, parameter :: ERROR_OVERLAP_INDEX_CODE                    = 4
        integer, parameter :: ERROR_OVERLAP_INCOMPATIBLE_CODE             = 5
        integer, parameter :: ERROR_BC_SECTION_TYPE_CODE                  = 6
        integer, parameter :: ERROR_GRADIENT_TYPE_CODE                    = 7
        integer, parameter :: ERROR_CPT_OVERLAP_INDEX_CODE                = 8
        integer, parameter :: ERROR_CPT_TYPE_CODE                         = 9
        integer, parameter :: ERROR_MAINLAYER_INTERFACE_TYPE_CODE         = 10
        integer, parameter :: ERROR_MAINLAYER_INTERFACE_INCOMPATIBLE_CODE = 11
        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if mainlayer_id.ne.(N,S,E,W)
        !
        !> @date
        !> 11_04_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_id
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if mainlayer_id do not match
        !
        !> @date
        !> 11_04_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_id1
        !> value of the parameter trigerring the exception
        !
        !>@param mainlayer_id2
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if the neighbor is incompatible
        !> with the current buffer layer
        !
        !> @date
        !> 11_04_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_id1
        !> value of the parameter trigerring the exception
        !
        !>@param mainlayer_id2
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error neighbor index.ne.(1,2)
        !
        !> @date
        !> 11_04_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param neighbor_index
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error overlap index.ne.(0,1,2,3,4,5,6,7,8)
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param overlap_index
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_overlap_index(
     $     file_name, fct_name, overlap_index)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: overlap_index

          write(ERROR_UNIT, '(I3)') ERROR_OVERLAP_INDEX_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'overlap_index not recognized'
          write(ERROR_UNIT, '(''overlap_index: '',I2)') overlap_index

          stop 'error_overlap_index'

        end subroutine error_overlap_index


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error overlap indices incompatible
        !< (N.and.S) (E.and.W)
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param overlap_index1
        !> value of the parameter trigerring the exception
        !
        !>@param overlap_index2
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_overlap_incompatible(
     $     file_name, fct_name, overlap_index1, overlap_index2)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: overlap_index1
          integer     , intent(in) :: overlap_index2

          write(ERROR_UNIT, '(I3)') ERROR_OVERLAP_INCOMPATIBLE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'overlap indices incompatible'
          write(ERROR_UNIT, '(''overlap_index1: '',I2)') overlap_index1
          write(ERROR_UNIT, '(''overlap_index2: '',I2)') overlap_index2

          stop 'error_overlap_incompatible'

        end subroutine error_overlap_incompatible


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error overlap index.ne.(0,1,2,3,4,5,6,7,8)
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param overlap_index
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_bc_section_type(
     $     file_name, fct_name, bc_section_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: bc_section_type

          write(ERROR_UNIT, '(I3)') ERROR_BC_SECTION_TYPE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'bc_section_type not recognized'
          write(ERROR_UNIT, '(''bc_section_type: '',I2)') bc_section_type

          stop 'error_bc_section_type'

        end subroutine error_bc_section_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error gradient_type
        !
        !> @date
        !> 02_03_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param gradient_type
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_gradient_type(
     $     file_name,
     $     fct_name,
     $     gradient_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: gradient_type

          write(ERROR_UNIT, '(I3)') ERROR_GRADIENT_TYPE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'gradient_type not recognized'
          write(ERROR_UNIT, '(''gradient_type: '',I2)') gradient_type

          stop 'error_gradient_type'

        end subroutine error_gradient_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !< error cpt_overlap_type
        !
        !> @date
        !> 02_03_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param cpt_overlap_type
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_cpt_overlap_index(
     $     file_name,
     $     fct_name,
     $     cpt_overlap_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: cpt_overlap_type

          write(ERROR_UNIT, '(I3)') ERROR_CPT_OVERLAP_INDEX_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'cpt_overlap_type not recognized'
          write(ERROR_UNIT, '(''cpt_overlap_type: '',I2)') cpt_overlap_type

          stop 'error_cpt_overlap_index'

        end subroutine error_cpt_overlap_index


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error cpt_type
        !
        !> @date
        !> 03_03_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param cpt_type
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_cpt_type(
     $     file_name,
     $     fct_name,
     $     cpt_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: cpt_type

          write(ERROR_UNIT, '(I3)') ERROR_CPT_TYPE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'cpt_type not recognized'
          write(ERROR_UNIT, '(''cpt_type: '',I2)') cpt_type

          stop 'error_cpt_type'

        end subroutine error_cpt_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error mainlayer_interface_type
        !
        !> @date
        !> 03_03_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param cpt_type
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_mainlayer_interface_type(
     $     file_name,
     $     fct_name,
     $     mainlayer_interface_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_interface_type

          write(ERROR_UNIT, '(I3)') ERROR_MAINLAYER_INTERFACE_TYPE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'mainlayer_interface_type not recognized'
          write(ERROR_UNIT, '(''mainlayer_interface_type: '',I2)') mainlayer_interface_type

          stop 'error_mainlayer_interface_type'

        end subroutine error_mainlayer_interface_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error incompatible mainlayer_interface
        !
        !> @date
        !> 09_03_2014 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param mainlayer_interface_type
        !> identification of the mainlayer_interface
        !
        !>@param mainlayer_localization
        !> cardinal coordinate of the mainlayer
        !--------------------------------------------------------------
        subroutine error_mainlayer_interface_incompatible(
     $     file_name, fct_name,
     $     mainlayer_interface_type, 
     $     mainlayer_localization)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: mainlayer_interface_type
          integer     , intent(in) :: mainlayer_localization

          write(ERROR_UNIT, '(I3)') ERROR_MAINLAYER_INTERFACE_INCOMPATIBLE_CODE
          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'incompatible mainlayer_interface and localization'
          write(ERROR_UNIT, '(''mainlayer_interface_type: '',I2)') mainlayer_interface_type
          write(ERROR_UNIT, '(''mainlayer_localization: '',I2)') mainlayer_localization
          
          stop 'error_mainlayer_interface_incompatible'  

        end subroutine error_mainlayer_interface_incompatible

      end module bf_layer_errors_module
