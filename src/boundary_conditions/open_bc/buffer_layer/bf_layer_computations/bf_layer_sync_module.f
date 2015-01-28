      !> @file
      !> procedures identifying the indices synchronizing the content
      !> of the interior data with a temporary array or the content of
      !> a buffer layer with a temporary array
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module with the subroutines identifying the indices
      !> synchronizing the content of the interior data with a
      !> temporary array or the content of a buffer layer with a
      !> temporary array
      !
      !> @date
      ! 28_01_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_layer_sync_module

        use bf_layer_errors_module, only :
     $     error_mainlayer_id

        use parameters_constant, only :
     $     N,S,E,W

        use parameters_input, only : 
     $     nx,ny,bc_size

        use parameters_kind, only :
     $     ikind

        implicit none


        private
        public :: 
     $       get_sync_indices_to_extract_interior_data,
     $       get_sync_indices_to_extract_bf_layer_data,
     $       get_bf_layer_match_table


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the synchronization indices to extract the content
        !> of the interior domain for the tmp array
        !    ___________________
        !   |                  _|_________
        !   |    interior     |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                      .
        !                     /|\
        !                      |
        !                   overlapping 
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param gen_coords
        !> coordinates of the four borders of the tmp array expressed as
        !> SW_corner = [gen_coords(1,1),gen_coords(2,1)]
        !> NE_corner = [gen_coords(1,2),gen_coords(2,2)]
        !> the coordinates are expressed in the general coordinate reference
        !> frame
        !
        !> @param size_x
        !> extent along the x-direction of the overlapping
        !
        !> @param size_y
        !> extent along the y-direction of the overlapping
        !
        !> @param i_recv
        !> x-index of the first gridpoint copied in the tmp array
        !
        !> @param j_recv
        !> y-index of the first gridpoint copied in the tmp array
        !
        !> @param i_send
        !> x-index of the first gridpoint copied from the interior array
        !
        !> @param j_send
        !> y-index of the first gridpoint copied from the interior array
        !--------------------------------------------------------------
        subroutine get_sync_indices_to_extract_interior_data(
     $       gen_coords,
     $       size_x, size_y,
     $       i_recv, j_recv,
     $       i_send, j_send)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: gen_coords
          integer(ikind)                , intent(out) :: size_x
          integer(ikind)                , intent(out) :: size_y
          integer(ikind)                , intent(out) :: i_recv
          integer(ikind)                , intent(out) :: j_recv
          integer(ikind)                , intent(out) :: i_send
          integer(ikind)                , intent(out) :: j_send

          
          integer(ikind) :: i_min,i_max
          integer(ikind) :: j_min,j_max


          i_min = max(1 ,gen_coords(1,1))
          i_max = min(nx,gen_coords(1,2))
          j_min = max(1 ,gen_coords(2,1))
          j_max = min(ny,gen_coords(2,2))

          size_x = i_max-i_min+1
          size_y = j_max-j_min+1 

          i_recv = i_min-gen_coords(1,1)+1
          i_send = i_min

          j_recv = j_min-gen_coords(2,1)+1
          j_send = j_min

        end subroutine get_sync_indices_to_extract_interior_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the synchronization indices to extract the content
        !> of the buffer layer for the tmp array
        !    ___________________
        !   |                  _|_________
        !   |   buffer layer  |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                      .
        !                     /|\
        !                      |
        !                   overlapping 
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_align
        !> relative position of the buffer layer to the interior domain
        !
        !> @param gen_coords
        !> coordinates of the four borders of the tmp array expressed as
        !> SW_corner = [gen_coords(1,1),gen_coords(2,1)]
        !> NE_corner = [gen_coords(1,2),gen_coords(2,2)]
        !> the coordinates are expressed in the general coordinate reference
        !> frame
        !
        !> @param size_x
        !> extent along the x-direction of the overlapping
        !
        !> @param size_y
        !> extent along the y-direction of the overlapping
        !
        !> @param i_recv
        !> x-index of the first gridpoint copied in the tmp array
        !
        !> @param j_recv
        !> y-index of the first gridpoint copied in the tmp array
        !
        !> @param i_send
        !> x-index of the first gridpoint copied from the interior array
        !
        !> @param j_send
        !> y-index of the first gridpoint copied from the interior array
        !--------------------------------------------------------------
        subroutine get_sync_indices_to_extract_bf_layer_data(
     $     bf_align,
     $     gen_coords,
     $     size_x, size_y,
     $     i_recv, j_recv,
     $     i_send, j_send)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_align
          integer(ikind), dimension(2,2), intent(in)  :: gen_coords
          integer(ikind)                , intent(out) :: size_x
          integer(ikind)                , intent(out) :: size_y
          integer(ikind)                , intent(out) :: i_recv
          integer(ikind)                , intent(out) :: j_recv
          integer(ikind)                , intent(out) :: i_send
          integer(ikind)                , intent(out) :: j_send

          integer(ikind) :: i_min, i_max, j_min, j_max


          i_min = max(bf_align(1,1)-bc_size,gen_coords(1,1))
          i_max = min(bf_align(1,2)+bc_size,gen_coords(1,2))
          j_min = max(bf_align(2,1)-bc_size,gen_coords(2,1))
          j_max = min(bf_align(2,2)+bc_size,gen_coords(2,2))

          size_x = i_max-i_min+1
          size_y = j_max-j_min+1 

          i_recv = i_min-gen_coords(1,1)+1
          i_send = i_min-(bf_align(1,1)-bc_size)+1

          j_recv = j_min-gen_coords(2,1)+1
          j_send = j_min-(bf_align(2,1)-bc_size)+1


        end subroutine get_sync_indices_to_extract_bf_layer_data


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the synchronization indices to extract the content
        !> of the buffer layer for the tmp array
        !    ___________________
        !   |                  _|_________
        !   |   buffer layer  |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                      .
        !                     /|\
        !                      |
        !                   overlapping 
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_localization
        !> cardinal coordinate identifying the position of the buffer layer
        !
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !
        !> @return match_table
        !> match_table allowing the conversion of local indices of the buffer
        !> layer into general indices matching the interior domain
        !> 
        !--------------------------------------------------------------
        function get_bf_layer_match_table(
     $     bf_alignment)
     $     result(match_table)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind), dimension(2)               :: match_table

          match_table(1) = bf_alignment(1,1) - bc_size - 1
          match_table(2) = bf_alignment(2,1) - bc_size - 1

        end function get_bf_layer_match_table

      end module bf_layer_sync_module
