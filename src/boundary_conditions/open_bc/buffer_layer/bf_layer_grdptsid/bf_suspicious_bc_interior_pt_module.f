      !> @file
      !> checking the neighbors of the suspicious bc_interior_pt and turn
      !> it into an interior_pt
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing subroutines useful for checking the neighbors
      !> of a suspicious bc_interior_pt that may be turned into an
      !> interior_pt
      !
      !> @date
      ! 27_11_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_suspicious_bc_interior_pt_module

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: verify_if_all_grdpts_exist
        

        contains


        !verify if all the grid points around a central grid points
        !exist in order to be considered an interior_pt
        function verify_if_all_grdpts_exist(i_c,j_c,grdpts_id_tmp)
     $       result(all_grdpts_exist)

          implicit none

          integer(ikind)         , intent(in) :: i_c
          integer(ikind)         , intent(in) :: j_c
          integer, dimension(:,:), intent(in) :: grdpts_id_tmp
          logical                             :: all_grdpts_exist


          integer(ikind) :: i
          integer(ikind) :: j


          all_grdpts_exist = .true.

          do j=j_c-bc_size,j_c+bc_size
             do i=i_c-bc_size,i_c+bc_size

                if(grdpts_id_tmp(i,j).eq.no_pt) then
                   all_grdpts_exist = .false.
                end if

             end do
          end do

        end function verify_if_all_grdpts_exist        

      end module bf_suspicious_bc_interior_pt_module
