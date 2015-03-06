      module bf_newgrdpt_extract_module

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_interior_data,
     $       get_grdpts_id_from_interior

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
     $       get_interior_data_for_newgrdpt,
     $       are_intermediate_newgrdpt_data_needed,
     $       get_x_map_for_newgrdpt,
     $       get_y_map_for_newgrdpt

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the grdpts_id, the coordinate maps and the
        !> nodes at t-dt and t corresponding to the general
        !> coordinates gen_coords
        !    ___________________
        !   |                  _|_________
        !   |    interior     |/|         |
        !   |                 |/|  tmp    |
        !   !                 !/!         !
        !                   overlapping which is copied
        !                     from buffer layer to tmp
        !
        !> @date
        !> 18_11_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_nodes0
        !> interior nodes at t=t-dt
        !
        !>@param interior_nodes1
        !> interior nodes at t=t
        !
        !>@param tmp_grdpts_id0
        !> array with the grdpts_id data
        !
        !>@param tmp_nodes0
        !> array with the grid points data at t-dt
        !
        !>@param tmp_nodes1
        !> array with the grid points data at t
        !
        !>@param gen_coords
        !> coordinates of the SW corner and the NE corners of the
        !> tmp arrays computed
        !--------------------------------------------------------------
        subroutine get_interior_data_for_newgrdpt(
     $     interior_nodes0,
     $     interior_nodes1,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_coords)

          implicit none

          real(rkind)    , dimension(nx,ny,ne)                          , intent(in)    :: interior_nodes0
          real(rkind)    , dimension(nx,ny,ne)                          , intent(in)    :: interior_nodes1
          integer        , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1)   , intent(inout) :: tmp_grdpts_id0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes1
          integer(ikind) , dimension(2,2)                               , intent(in)    :: gen_coords


          integer(ikind) :: size_x,size_y
          integer(ikind) :: i_recv,i_send,j_recv,j_send
          integer(ikind) :: i,j
          integer        :: k


          !synchronize the nodes
          call get_indices_to_extract_interior_data(
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)


          do k=1,ne
             do j=1, size_y
                do i=1, size_x

                   tmp_nodes0(i_recv+i-1,j_recv+j-1,k) =
     $                  interior_nodes0(i_send+i-1,j_send+j-1,k)

                   tmp_nodes1(i_recv+i-1,j_recv+j-1,k) =
     $                  interior_nodes1(i_send+i-1,j_send+j-1,k)

                end do
             end do
          end do


          !synchronize the grdpts_id
          !initialize with no_pt
          call get_grdpts_id_from_interior(
     $         tmp_grdpts_id0,
     $         gen_coords)          

        end subroutine get_interior_data_for_newgrdpt        


        function are_intermediate_newgrdpt_data_needed(
     $     mainlayer_id, gen_coords)
     $     result(needed)

          implicit none

          integer                       , intent(in) :: mainlayer_id
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          logical                                    :: needed
          

          select case(mainlayer_id)
            case(N)
               needed = .not.(gen_coords(2,1).gt.ny)

            case(S)
               needed = .not.(gen_coords(2,2).lt.1)
               
            case(E)
               needed = .not.(
     $              (gen_coords(1,1).gt.nx).and.
     $              (gen_coords(2,1).gt.bc_size).and.
     $              (gen_coords(2,2).lt.(ny-bc_size+1)))
               
            case(W)
               needed = .not.(
     $              (gen_coords(1,2).lt.1).and.
     $              (gen_coords(2,1).gt.bc_size).and.
     $              (gen_coords(2,2).lt.(ny-bc_size+1)))
               
            case default
               call error_mainlayer_id(
     $              'bf_layer_newgrdpt_procedure_module',
     $              'are_intermediate_newgrdpt_data_needed',
     $              mainlayer_id)

          end select

        end function are_intermediate_newgrdpt_data_needed


        function get_x_map_for_newgrdpt(
     $     interior_x_map, gen_coords)
     $     result(tmp_x_map)

          implicit none

          real(rkind)   , dimension(nx) , intent(in) :: interior_x_map
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_x_map

          real(rkind) :: dx
          integer     :: i
          integer     :: size

          ! --[-----]---|-------|--------
          if(gen_coords(1,2).le.0) then

             dx = interior_x_map(2)-interior_x_map(1)

             do i=gen_coords(1,1),gen_coords(1,2)
                
                tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(1) + (i-1)*dx

             end do

          else
          ! --------[--|--]-----|--------
             if(gen_coords(1,1).le.0) then

                dx = interior_x_map(2)-interior_x_map(1)

                do i=gen_coords(1,1),0
                   tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(1) + (i-1)*dx
                end do

                size = -gen_coords(1,1) + 1

                do i=1,gen_coords(1,2)
                   tmp_x_map(size+i) = interior_x_map(i)
                end do

          ! -----------|-[-----]|--------
             else
                if(gen_coords(1,2).le.nx) then

                   do i=gen_coords(1,1),gen_coords(1,2)
                      tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(i)
                   end do

          ! -----------|-----[--|--]-----
                else

                   dx = interior_x_map(nx) - interior_x_map(nx-1)

                   if(gen_coords(1,1).le.nx) then
                      
                      do i=gen_coords(1,1),nx
                         tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(i)
                      end do

                      size = -gen_coords(1,1) + 1

                      do i=nx+1,gen_coords(1,2)
                         tmp_x_map(size+i) = interior_x_map(nx) + (i-nx)*dx
                      end do

           ! -----------|--------|-[-----]
                   else

                      do i=gen_coords(1,1),gen_coords(1,2)
                         tmp_x_map(i-gen_coords(1,1)+1) = interior_x_map(nx) + (i-nx)*dx
                      end do

                   end if
                end if
             end if
          end if

        end function get_x_map_for_newgrdpt


        function get_y_map_for_newgrdpt(
     $     interior_y_map, gen_coords)
     $     result(tmp_y_map)

          implicit none

          real(rkind)   , dimension(ny) , intent(in) :: interior_y_map
          integer(ikind), dimension(2,2), intent(in) :: gen_coords
          real(rkind)   , dimension(2*(bc_size+1)+1) :: tmp_y_map

          real(rkind) :: dy
          integer     :: i
          integer     :: size

          ! --[-----]---|-------|--------
          if(gen_coords(2,2).le.0) then

             dy = interior_y_map(2)-interior_y_map(1)

             do i=gen_coords(2,1),gen_coords(2,2)
                
                tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(1) + (i-1)*dy

             end do

          else
          ! --------[--|--]-----|--------
             if(gen_coords(2,1).le.0) then

                dy = interior_y_map(2)-interior_y_map(1)

                do i=gen_coords(2,1),0
                   tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(1) + (i-1)*dy
                end do

                size = -gen_coords(2,1) + 1

                do i=1,gen_coords(2,2)
                   tmp_y_map(size+i) = interior_y_map(i)
                end do

          ! -----------|-[-----]|--------
             else
                if(gen_coords(2,2).le.ny) then

                   do i=gen_coords(2,1),gen_coords(2,2)
                      tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(i)
                   end do

          ! -----------|-----[--|--]-----
                else

                   dy = interior_y_map(ny) - interior_y_map(ny-1)

                   if(gen_coords(2,1).le.ny) then
                      
                      do i=gen_coords(2,1),ny
                         tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(i)
                      end do

                      size = -gen_coords(2,1) + 1

                      do i=ny+1,gen_coords(2,2)
                         tmp_y_map(size+i) = interior_y_map(ny) + (i-ny)*dy
                      end do

           ! -----------|--------|-[-----]
                   else

                      do i=gen_coords(2,1),gen_coords(2,2)
                         tmp_y_map(i-gen_coords(2,1)+1) = interior_y_map(ny) + (i-ny)*dy
                      end do

                   end if
                end if
             end if
          end if

        end function get_y_map_for_newgrdpt

      end module bf_newgrdpt_extract_module
