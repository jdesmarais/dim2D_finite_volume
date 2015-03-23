      !> @file
      !> sorting functions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> sorting functions
      !
      !> @date
      !> 23_03_2015 - intial version  - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_sorting_module

        use parameters_kind, only :
     $     ikind

        implicit none

        private
        public ::
     $       bubble_sort_grdpts


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> use bubble sorting to order the coordinates of the grdpts_id
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param a
        !> bc_section array sorted
        !--------------------------------------------------------------
        subroutine bubble_sort_grdpts(a)

          implicit none

          integer(ikind), dimension(:,:), intent(inout) :: a

          integer :: i,j
          integer :: n
          integer :: max_j

          n = size(a,2)

          do i=1,n

             max_j = 1

             do j=2, n-i+1
                if(order_grdpts(a(:,j),a(:,max_j))) then
                   max_j = j
                end if
             end do

             call exchange_grdpts(a(:,max_j),a(:,n-i+1))

          end do

        end subroutine bubble_sort_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> operator ordering two grid-points
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param p
        !> general coordinates of a grid-point [i,j]
        !
        !> @param q
        !> general coordinates of a grid-point [i,j]
        !
        !> @return p_larger_than_q
        !> logical indicating whether p>q
        !--------------------------------------------------------------
        function order_grdpts(p,q)
     $     result(p_larger_than_q)
        
          implicit none

          integer(ikind), dimension(2), intent(in) :: p
          integer(ikind), dimension(2), intent(in) :: q
          logical                                  :: p_larger_than_q

          p_larger_than_q = .false.

          !if(p(j)>q(j))
          if (p(2).gt.q(2)) then
             p_larger_than_q = .true.
          else
             if(p(2).eq.q(2)) then

                !if(p(i)>q(i))
                if(p(1).gt.q(1)) then
                   p_larger_than_q = .true.
                end if

             end if
          end if

        end function order_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> permutation of two bc_sections
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !> @param p
        !> general coordinates of a grid-point [i,j]
        !
        !> @param q
        !> general coordinates of a grid-point [i,j]
        !--------------------------------------------------------------
        subroutine exchange_grdpts(p,q)

          implicit none

          integer(ikind), dimension(2), intent(inout) :: p
          integer(ikind), dimension(2), intent(inout) :: q
          integer(ikind), dimension(2)                :: temp

          temp = p
          p    = q
          q    = temp

        end subroutine exchange_grdpts

      end module bf_sorting_module
