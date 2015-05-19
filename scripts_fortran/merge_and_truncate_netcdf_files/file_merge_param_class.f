      module file_merge_param_class

        use parameters_kind, only :
     $     rkind

        implicit none


        private
        public :: file_merge_param


        type :: file_merge_param

          integer               :: rank
          integer, dimension(2) :: cartesian_comm_ranks

          integer, dimension(2,2) :: extracted_from_file
          integer, dimension(2,2) :: saved_in_lg_domain

          contains

          procedure, pass :: get_rank
          procedure, pass :: get_indices_for_extracting
          procedure, pass :: get_indices_for_saving

          procedure, pass :: set_ranks
          procedure, pass :: determine_indices_for_extracting_and_saving

        end type file_merge_param


        contains


        function get_rank(this)

          implicit none

          class(file_merge_param), intent(in) :: this
          integer                             :: get_rank

          get_rank = this%rank

        end function get_rank


        function get_indices_for_extracting(this)

          implicit none

          class(file_merge_param), intent(in) :: this
          integer, dimension(2,2)             :: get_indices_for_extracting

          get_indices_for_extracting = this%extracted_from_file

        end function get_indices_for_extracting


        function get_indices_for_saving(this)

          implicit none

          class(file_merge_param), intent(in) :: this
          integer, dimension(2,2)             :: get_indices_for_saving

          get_indices_for_saving = this%saved_in_lg_domain

        end function get_indices_for_saving


        subroutine set_ranks(this,rank_x,rank_y,nb_tiles_y)

          implicit none

          class(file_merge_param), intent(inout) :: this
          integer                , intent(in)    :: rank_x
          integer                , intent(in)    :: rank_y
          integer                , intent(in)    :: nb_tiles_y

          this%cartesian_comm_ranks(1) = rank_x
          this%cartesian_comm_ranks(2) = rank_y
          this%rank                    = nb_tiles_y*rank_x + rank_y

        end subroutine set_ranks


        subroutine determine_indices_for_extracting_and_saving(
     $     this,
     $     borders_lg_domain,
     $     borders_extraction,
     $     grid_spacings,
     $     nb_pts_per_tile,
     $     bc_size)

          implicit none

          class(file_merge_param)    , intent(inout) :: this
          real(rkind), dimension(2,2), intent(in)    :: borders_lg_domain
          real(rkind), dimension(2,2), intent(in)    :: borders_extraction
          real(rkind), dimension(2)  , intent(in)    :: grid_spacings
          integer    , dimension(2)  , intent(in)    :: nb_pts_per_tile
          integer                    , intent(in)    :: bc_size

          real(rkind), dimension(2,2) :: borders_tile
          integer                     :: j
          integer                     :: k
          real(rkind)                 :: point


          ! determine the borders of the current tile
          do j=1,2

             borders_tile(j,1) = borders_lg_domain(j,1) +
     $            this%cartesian_comm_ranks(j)*(nb_pts_per_tile(j)-2*bc_size)*
     $            grid_spacings(j)

             borders_tile(j,2) = borders_tile(j,1) +
     $            (nb_pts_per_tile(j)-1)*grid_spacings(j)

          end do

          print '(''borders_tile-x: '',2F7.4)', borders_tile(1,:)
          print '(''borders_tile-y: '',2F7.4)', borders_tile(2,:)


          ! determine the extraction indices
          do j=1,2
             do k=1,2

                this%extracted_from_file(j,k) =
     $               nint((borders_extraction(j,k)-borders_tile(j,1))/
     $               grid_spacings(j))
     $               +1
                
             end do

             this%extracted_from_file(j,1) = max(1,this%extracted_from_file(j,1))
             this%extracted_from_file(j,2) = min(nb_pts_per_tile(j),this%extracted_from_file(j,2))

          end do


          ! determine the indices for saving in the file
          do k=1,2
             do j=1,2

                point = borders_tile(j,1)+
     $                  (this%extracted_from_file(j,k)-1)*grid_spacings(j)

                this%saved_in_lg_domain(j,k) =
     $               nint((point-borders_extraction(j,1))/grid_spacings(j))+1

             end do
          end do

        end subroutine determine_indices_for_extracting_and_saving

      end module file_merge_param_class
