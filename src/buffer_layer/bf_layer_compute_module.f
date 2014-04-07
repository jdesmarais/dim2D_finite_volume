      !> @file
      !> module encapsulating the subroutines related to
      !> computation of new grid points for the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to computation of new grid points
      !> for the buffer layer
      !
      !> @date
      ! 04_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_compute_module

        use bf_layer_abstract_class , only : bf_layer_abstract
        use parameters_input        , only : ne
        use parameters_kind         , only : ikind

        implicit none

        private
        public :: compute_new_grdpts


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the new grid points of
        !> the buffer layers using the boundary conditions
        !> applied at the edges
        !
        !> @date
        !> 04_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param list_new_grdpts
        !> table of integers identifying in the buffer layer
        !> the new grid points that should be computed
        !--------------------------------------------------------------
        subroutine compute_new_grdpts(this, list_new_grdpts)

          implicit none

          class(bf_layer_abstract)      , intent(inout) :: this
          integer(ikind), dimension(:,:), intent(in)    :: list_new_grdpts
          
          integer(ikind) :: i,j,k,l

          do l=1, size(list_new_grdpts,1)
             
             i = list_new_grdpts(l,1)
             j = list_new_grdpts(l,2)
          
             do k=1, ne
                this%nodes(i,j,k) = 3
             end do

          end do

        end subroutine compute_new_grdpts

      end module bf_layer_compute_module
