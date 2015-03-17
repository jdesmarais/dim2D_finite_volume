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
      module bf_layer_extract_module

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use parameters_bf_layer, only :
     $       no_pt,bc_pt,bc_interior_pt,interior_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only : 
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: 
     $       get_indices_to_extract_interior_data,
     $       get_indices_to_extract_bf_layer_data,
     $       get_bf_layer_match_table,
     $       get_grdpts_id_from_interior,
     $       get_grdpts_id_from_bf_layer,
     $       get_nodes_from_interior,
     $       get_nodes_from_bf_layer,
     $       get_map_from_interior,
     $       set_grdpts_id_in_bf_layer


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
        subroutine get_indices_to_extract_interior_data(
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

        end subroutine get_indices_to_extract_interior_data


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
        subroutine get_indices_to_extract_bf_layer_data(
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


        end subroutine get_indices_to_extract_bf_layer_data


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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id corresponding to the general
        !> coordinates gen_coords for the interior
        !    ___________________
        !   |                  _|_________
        !   |    interior     |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> tmp arrays computed
        !--------------------------------------------------------------
        subroutine get_grdpts_id_from_interior(
     $     tmp_grdpts_id,
     $     gen_coords)
        
          implicit none
          
          integer       , dimension(:,:), intent(out) :: tmp_grdpts_id
          integer(ikind), dimension(2,2), intent(in)  :: gen_coords

          integer                        :: i,j
          integer(ikind), dimension(2,2) :: grdpts_id_coords

          do j=1, size(tmp_grdpts_id,2)
             do i=1, size(tmp_grdpts_id,1)
                tmp_grdpts_id(i,j) = no_pt
             end do
          end do

          !overlap with the bc_pt if any
          grdpts_id_coords(1,1) = 1
          grdpts_id_coords(1,2) = nx
          grdpts_id_coords(2,1) = 1
          grdpts_id_coords(2,2) = ny

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         bc_pt)

          !overlap with the bc_interior_pt if any
          grdpts_id_coords(1,1) = 2
          grdpts_id_coords(1,2) = nx-1
          grdpts_id_coords(2,1) = 2
          grdpts_id_coords(2,2) = ny-1

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         bc_interior_pt)

          !overlap with the interior_pt if any
          grdpts_id_coords(1,1) = bc_size+1
          grdpts_id_coords(1,2) = nx-bc_size
          grdpts_id_coords(2,1) = bc_size+1
          grdpts_id_coords(2,2) = ny-bc_size

          call set_grdptid(
     $         tmp_grdpts_id,
     $         grdpts_id_coords,
     $         gen_coords,
     $         interior_pt)

        end subroutine get_grdpts_id_from_interior

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the content of grdpts_id delimited by the domain
        !> grdpts_id_coords intersecting gen_coords to the constant
        !> pt_type
        !
        !> @date
        !> 27_11_2014 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data
        !
        !>@param grdpts_id_coords
        !> coordinates of the SW and NE corners of the first domain
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> second domain
        !
        !>@param pt_type
        !> constant set in the grdpts_id array
        !--------------------------------------------------------------
        subroutine set_grdptid(
     $     grdpts_id,
     $     grdpts_id_coords,
     $     gen_coords,
     $     pt_type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2), intent(in)    :: grdpts_id_coords
          integer(ikind), dimension(2,2), intent(in)    :: gen_coords
          integer                       , intent(in)    :: pt_type


          integer(ikind) :: i_min,i_max,j_min,j_max
          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,j_recv
          integer(ikind) :: i,j


          i_min = max(grdpts_id_coords(1,1), gen_coords(1,1))
          i_max = min(grdpts_id_coords(1,2), gen_coords(1,2))
          j_min = max(grdpts_id_coords(2,1), gen_coords(2,1))
          j_max = min(grdpts_id_coords(2,2), gen_coords(2,2))

          size_x = i_max-i_min+1
          size_y = j_max-j_min+1 

          i_recv = i_min-gen_coords(1,1)+1
          j_recv = j_min-gen_coords(2,1)+1

          
          do j=1,size_y
             do i=1,size_x
                grdpts_id(i_recv+i-1,j_recv+j-1) = pt_type
             end do
          end do

        end subroutine set_grdptid


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract grdpts_id from a buffer layer
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_grdpts_id
        !> grdpts_id of the buffer layer
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !--------------------------------------------------------------
        subroutine get_grdpts_id_from_bf_layer(
     $     tmp_grdpts_id,
     $     gen_coords,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     extract_param_in,
     $     extract_param_out)

          implicit none

          integer       , dimension(:,:)          , intent(inout) :: tmp_grdpts_id
          integer(ikind), dimension(2,2)          , intent(in)    :: gen_coords
          integer(ikind), dimension(2,2)          , intent(in)    :: bf_alignment
          integer       , dimension(:,:)          , intent(in)    :: bf_grdpts_id
          integer(ikind), dimension(6)  , optional, intent(in)    :: extract_param_in
          integer(ikind), dimension(6)  , optional, intent(out)   :: extract_param_out

          integer(ikind) :: size_x
          integer(ikind) :: size_y
          integer(ikind) :: i_recv
          integer(ikind) :: i_send
          integer(ikind) :: j_recv
          integer(ikind) :: j_send
          integer(ikind) :: i
          integer(ikind) :: j
          
          
          if(present(extract_param_in)) then
             size_x = extract_param_in(1)
             size_y = extract_param_in(2)
             i_recv = extract_param_in(3)
             j_recv = extract_param_in(4)
             i_send = extract_param_in(5)
             j_send = extract_param_in(6)

          else
             call get_indices_to_extract_bf_layer_data(
     $            bf_alignment,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)

             if(present(extract_param_out)) then
                extract_param_out(1) = size_x
                extract_param_out(2) = size_y
                extract_param_out(3) = i_recv
                extract_param_out(4) = j_recv
                extract_param_out(5) = i_send
                extract_param_out(6) = j_send
             end if

          end if

          do j=1, size_y
             do i=1, size_x
                tmp_grdpts_id(i_recv+i-1,j_recv+j-1) =
     $               bf_grdpts_id(i_send+i-1,j_send+j-1)
             end do
          end do

        end subroutine get_grdpts_id_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract nodes from the interior domain
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_nodes
        !> array with the nodes extracted
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !--------------------------------------------------------------
        subroutine get_nodes_from_interior(
     $     tmp_nodes,
     $     gen_coords,
     $     interior_nodes,
     $     extract_param_in,
     $     extract_param_out)

          implicit none

          real(rkind)   , dimension(:,:,:)             , intent(inout) :: tmp_nodes
          integer(ikind), dimension(2,2)               , intent(in)    :: gen_coords
          real(rkind)   , dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer(ikind), dimension(6)       , optional, intent(in)    :: extract_param_in
          integer(ikind), dimension(6)       , optional, intent(out)   :: extract_param_out

          integer(ikind), dimension(2,2) :: alignment


          alignment = reshape((/bc_size+1,bc_size+1,nx-bc_size,ny-bc_size/),(/2,2/))

          
          if(present(extract_param_in)) then
             call get_nodes_from_bf_layer(
     $            tmp_nodes,
     $            gen_coords,
     $            alignment,
     $            interior_nodes,
     $            extract_param_in=extract_param_in)

          else

             if(present(extract_param_out)) then
                call get_nodes_from_bf_layer(
     $               tmp_nodes,
     $               gen_coords,
     $               alignment,
     $               interior_nodes,
     $               extract_param_out=extract_param_out)

             else
                call get_nodes_from_bf_layer(
     $               tmp_nodes,
     $               gen_coords,
     $               alignment,
     $               interior_nodes)

             end if

          end if

        end subroutine get_nodes_from_interior



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract nodes from a buffer layer
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_nodes
        !> array with the nodes extracted
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_nodes
        !> nodes of the buffer layer
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !--------------------------------------------------------------
        subroutine get_nodes_from_bf_layer(
     $     tmp_nodes,
     $     gen_coords,
     $     bf_alignment,
     $     bf_nodes,
     $     extract_param_in,
     $     extract_param_out)

          implicit none

          real(rkind)   , dimension(:,:,:)          , intent(inout) :: tmp_nodes
          integer(ikind), dimension(2,2)            , intent(in)    :: gen_coords
          integer(ikind), dimension(2,2)            , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)          , intent(in)    :: bf_nodes
          integer(ikind), dimension(6)    , optional, intent(in)    :: extract_param_in
          integer(ikind), dimension(6)    , optional, intent(out)   :: extract_param_out

          integer(ikind) :: size_x
          integer(ikind) :: size_y
          integer(ikind) :: i_recv
          integer(ikind) :: i_send
          integer(ikind) :: j_recv
          integer(ikind) :: j_send
          integer(ikind) :: i
          integer(ikind) :: j
          integer        :: k
          
          
          if(present(extract_param_in)) then
             size_x = extract_param_in(1)
             size_y = extract_param_in(2)
             i_recv = extract_param_in(3)
             j_recv = extract_param_in(4)
             i_send = extract_param_in(5)
             j_send = extract_param_in(6)

          else
             call get_indices_to_extract_bf_layer_data(
     $            bf_alignment,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)

             if(present(extract_param_out)) then
                extract_param_out(1) = size_x
                extract_param_out(2) = size_y
                extract_param_out(3) = i_recv
                extract_param_out(4) = j_recv
                extract_param_out(5) = i_send
                extract_param_out(6) = j_send
             end if

          end if

          do k=1,ne
             do j=1, size_y
                do i=1, size_x
                   tmp_nodes(i_recv+i-1,j_recv+j-1,k) =
     $                  bf_nodes(i_send+i-1,j_send+j-1,k)
                end do
             end do
          end do          

        end subroutine get_nodes_from_bf_layer



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the map corresponding to the gen_coords
        !> if the interior map is provided
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_map
        !> array with the nodes extracted
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain extracted
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_nodes
        !> nodes of the buffer layer
        !--------------------------------------------------------------
        subroutine get_map_from_interior(
     $     tmp_map,
     $     gen_coords,
     $     interior_s_map)

          implicit none

          real(rkind), dimension(:), intent(out) :: tmp_map
          integer    , dimension(2), intent(in)  :: gen_coords
          real(rkind), dimension(:), intent(in)  :: interior_s_map

          real(rkind) :: ds
          integer     :: i
          integer     :: size_tmp
          integer     :: size_map


          size_map = size(interior_s_map,1)


          ! --[-----]---|-------|--------
          if(gen_coords(2).le.0) then

             ds = interior_s_map(2)-interior_s_map(1)

             do i=gen_coords(1),gen_coords(2)
                
                tmp_map(i-gen_coords(1)+1) = interior_s_map(1) + (i-1)*ds

             end do

          else
          ! --------[--|--]-----|--------
             if(gen_coords(1).le.0) then

                ds = interior_s_map(2)-interior_s_map(1)

                do i=gen_coords(1),0
                   tmp_map(i-gen_coords(1)+1) = interior_s_map(1) + (i-1)*ds
                end do

                size_tmp = -gen_coords(1) + 1

                if(gen_coords(2).le.size_map) then
                   do i=1,gen_coords(2)
                      tmp_map(size_tmp+i) = interior_s_map(i)
                   end do

                else
                   do i=1,size_map
                      tmp_map(size_tmp+i) = interior_s_map(i)
                   end do

                   size_tmp = size_map+size_tmp

                   ds = interior_s_map(size_map) -
     $                  interior_s_map(size_map-1)

                   do i=size_map+1,gen_coords(2)
                      tmp_map(size_tmp+i-size_map) = interior_s_map(size_map) +
     $                                               (i-size_map)*ds
                   end do

                end if


          ! -----------|-[-----]|--------
             else
                if(gen_coords(2).le.size_map) then

                   do i=gen_coords(1),gen_coords(2)
                      tmp_map(i-gen_coords(1)+1) = interior_s_map(i)
                   end do

          ! -----------|-----[--|--]-----
                else

                   ds = interior_s_map(size_map) -
     $                  interior_s_map(size_map-1)

                   if(gen_coords(1).le.size_map) then
                      
                      do i=gen_coords(1),size_map
                         tmp_map(i-gen_coords(1)+1) = interior_s_map(i)
                      end do

                      size_tmp = size_map - gen_coords(1) + 1

                      do i=size_map+1,gen_coords(2)
                         tmp_map(size_tmp+i-size_map) = interior_s_map(size_map) +
     $                                         (i-size_map)*ds
                      end do

           ! -----------|--------|-[-----]
                   else

                      do i=gen_coords(1),gen_coords(2)
                         tmp_map(i-gen_coords(1)+1) = interior_s_map(size_map) +
     $                                                (i-size_map)*ds
                      end do

                   end if
                end if
             end if
          end if


        end subroutine get_map_from_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> insert grdpts_id in a buffer layer
        !
        !> @date
        !> 18_03_2015 - initial version - J.L. Desmarais
        !
        !>@param tmp_grdpts_id
        !> array with the grdpts_id data
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> domain inserted
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_grdpts_id
        !> grdpts_id of the buffer layer
        !
        !>@param extract_param_in
        !> optional argument to avoid computing the parameters
        !> needed for the extraction
        !
        !>@param extract_param_out
        !> optional argument to get the parameters needed for the
        !> extraction
        !--------------------------------------------------------------
        subroutine set_grdpts_id_in_bf_layer(
     $     tmp_grdpts_id,
     $     gen_coords,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     insert_param_in,
     $     insert_param_out)

          implicit none

          integer       , dimension(:,:)          , intent(in)    :: tmp_grdpts_id
          integer(ikind), dimension(2,2)          , intent(in)    :: gen_coords
          integer(ikind), dimension(2,2)          , intent(in)    :: bf_alignment
          integer       , dimension(:,:)          , intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(6)  , optional, intent(in)    :: insert_param_in
          integer(ikind), dimension(6)  , optional, intent(out)   :: insert_param_out

          integer(ikind) :: size_x
          integer(ikind) :: size_y
          integer(ikind) :: i_recv
          integer(ikind) :: i_send
          integer(ikind) :: j_recv
          integer(ikind) :: j_send
          integer(ikind) :: i
          integer(ikind) :: j
          
          
          if(present(insert_param_in)) then
             size_x = insert_param_in(1)
             size_y = insert_param_in(2)
             i_recv = insert_param_in(3)
             j_recv = insert_param_in(4)
             i_send = insert_param_in(5)
             j_send = insert_param_in(6)

          else
             call get_indices_to_extract_bf_layer_data(
     $            bf_alignment,
     $            gen_coords,
     $            size_x, size_y,
     $            i_recv, j_recv,
     $            i_send, j_send)

             if(present(insert_param_out)) then
                insert_param_out(1) = size_x
                insert_param_out(2) = size_y
                insert_param_out(3) = i_recv
                insert_param_out(4) = j_recv
                insert_param_out(5) = i_send
                insert_param_out(6) = j_send
             end if

          end if

          do j=1, size_y
             do i=1, size_x
                bf_grdpts_id(i_send+i-1,j_send+j-1) = 
     $               tmp_grdpts_id(i_recv+i-1,j_recv+j-1)
             end do
          end do

        end subroutine set_grdpts_id_in_bf_layer

      end module bf_layer_extract_module
