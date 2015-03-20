      !> @file
      !> module encapsulating the subroutines related to the
      !> exchange of data between the different buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> subroutines related to the exchange of data between
      !> the different buffer layers
      !
      !> @date
      ! 27_05_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_exchange_module

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,y_direction

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind
      
        implicit none

        private
        public ::
     $       do_grdpts_overlap_along_x_dir,
     $       get_match_indices_for_exchange_with_neighbor1,
     $       get_match_indices_for_exchange_with_neighbor2,
     $       copy_from_bf1_to_bf2,
     $       get_sync_indices_with_interior,
     $       sync_nodes,
     $       get_sync_indices_with_neighbor1,
     $       get_sync_indices_with_neighbor2,
     $       update_alignment_for_exchanges

        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy grid points and their role from the bf1 to bf2
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param bf1_i_min
        !> min x-index in copying the grid points from bf1
        !
        !>@param bf1_j_min
        !> min y-index in copying the grid points from bf1
        !
        !>@param bf2_i_min
        !> min x-index in receiving the grid points in bf2
        !
        !>@param bf2_j_min
        !> min y-index in receiving the grid points in bf2
        !
        !>@param bf_copy_size_x
        !> extend of the layer copied in the x-direction
        !
        !>@param bf_copy_size_y
        !> extend of the layer copied in the y-direction
        !
        !>@param bf1_nodes
        !> nodes array for the buffer layer bf1
        !
        !>@param bf1_grdpts_id
        !> grdpts_id array for the buffer layer bf1
        !
        !>@param bf2_nodes
        !> nodes array for the buffer layer bf2
        !
        !>@param bf2_grdpts_id
        !> grdpts_id array for the buffer layer bf2
        !--------------------------------------------------------------
        subroutine copy_from_bf1_to_bf2(
     $       bf1_i_min, bf1_j_min, bf2_i_min, bf2_j_min,
     $       bf_copy_size_x, bf_copy_size_y,
     $       bf1_nodes, bf1_grdpts_id,
     $       bf2_nodes, bf2_grdpts_id)


          implicit none

          integer(ikind)                  , intent(in) :: bf1_i_min
          integer(ikind)                  , intent(in) :: bf1_j_min
          integer(ikind)                  , intent(in) :: bf2_i_min
          integer(ikind)                  , intent(in) :: bf2_j_min
          integer(ikind)                  , intent(in) :: bf_copy_size_x
          integer(ikind)                  , intent(in) :: bf_copy_size_y
          real(rkind)   , dimension(:,:,:), intent(in) :: bf1_nodes
          integer       , dimension(:,:)  , intent(in) :: bf1_grdpts_id
          real(rkind)   , dimension(:,:,:), intent(out):: bf2_nodes
          integer       , dimension(:,:)  , intent(out):: bf2_grdpts_id

          integer(ikind) :: i,j
          integer        :: k

          !nodes copy
          do k=1, ne
             do j=1, bf_copy_size_y
                do i=1, bf_copy_size_x
                   bf2_nodes(bf2_i_min+i-1,
     $                       bf2_j_min+j-1,k) =
     $             bf1_nodes(bf1_i_min+i-1,
     $                       bf1_j_min+j-1,k)
                end do
             end do
          end do

          !grdpts_id copy
          do j=1,bf_copy_size_y
             do i=1, bf_copy_size_x
                bf2_grdpts_id(
     $               bf2_i_min+i-1,
     $               bf2_j_min+j-1) =
     $          bf1_grdpts_id(
     $               bf1_i_min+i-1,
     $               bf1_j_min+j-1)
             end do
          end do

        end subroutine copy_from_bf1_to_bf2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying the layer to be exchanged
        !> between a buffer layer and its neighbor of type 1
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param localization
        !> cardinal coordinate identifying the position of the
        !> first buffer layer
        !
        !>@param bf_alignment
        !> alignment of the first buffer layer
        !
        !>@param bf_size_y
        !> extent in the y-direction of the first buffer layer arrays
        !
        !>@param nbf_alignment
        !> alignment of the neighboring buffer layer of type 1
        !
        !>@param nbf_size_y
        !> extent in the y-direction of the arrays of the neighboring
        !> buffer layer of type 1
        !
        !>@param bf_i_min
        !> min x-index in copying the grid points from the first
        !> buffer layer
        !
        !>@param bf_j_min
        !> min y-index in copying the grid points from the first
        !> buffer layer
        !
        !>@param nbf_i_min
        !> min x-index in receiving the grid points in the neighbor
        !> buffer layer of type 1
        !
        !>@param nbf_j_min
        !> min y-index in receiving the grid points in the neighbor
        !> buffer layer of type 1
        !
        !>@param bf_copy_size_x
        !> extend of the layer copied in the x-direction
        !
        !>@param bf_copy_size_y
        !> extend of the layer copied in the y-direction
        !--------------------------------------------------------------
        subroutine get_match_indices_for_exchange_with_neighbor1(
     $     bf_alignment,
     $     nbf_alignment,
     $     bf_i_min,
     $     bf_j_min,
     $     nbf_i_min,
     $     nbf_j_min,
     $     bf_copy_size_x,
     $     bf_copy_size_y)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer(ikind)                , intent(out) :: bf_i_min
          integer(ikind)                , intent(out) :: bf_j_min
          integer(ikind)                , intent(out) :: nbf_i_min
          integer(ikind)                , intent(out) :: nbf_j_min
          integer(ikind)                , intent(out) :: bf_copy_size_x
          integer(ikind)                , intent(out) :: bf_copy_size_y

          !we need to define the borders identifying the subarrays
          !copied from the tables of neighbors1 to the tables of the
          !current sublayer

          !for the borders along the x-direction, we first use their
          !general coordinates to identify them (that is the coordinates
          !with references to the interior domain)

          !then they are expressed as (bf_i_min, bf_i_max), and 
          !(nbf_i_min, nbf_i_max) for the local coordinates of the tables

          !the borders along the y-direction, we use the assumptions:
          ! - for the N and S buffer layers, the neighbor1 corresponds
          !   to the neighbor on the W main layer
          ! - for the E and W buffer layers, the neighbor1 corresponds
          !   to the neighbor on the S main layer
          
          !using the previous assumptions, it is possible to determine the
          !borders along the y-direction as local coordinates for the tables
          !(bf_j_min, bf_j_max) and (nbf_j_min, nbf_j_max)
          call get_exchange_indices(
     $         bf_alignment, nbf_alignment,
     $         x_direction,
     $         bf_i_min, nbf_i_min,
     $         bf_copy_size_x)

          call get_exchange_indices(
     $         bf_alignment, nbf_alignment,
     $         y_direction,
     $         bf_j_min, nbf_j_min,
     $         bf_copy_size_y)
          
        end subroutine get_match_indices_for_exchange_with_neighbor1    


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying the layer to be exchanged
        !> between a buffer layer and its neighbor of type 2
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param localization
        !> cardinal coordinate identifying the position of the
        !> first buffer layer
        !
        !>@param bf_alignment
        !> alignment of the first buffer layer
        !
        !>@param bf_size_y
        !> extent in the y-direction of the first buffer layer arrays
        !
        !>@param nbf_alignment
        !> alignment of the neighboring buffer layer of type 2
        !
        !>@param nbf_size_y
        !> extent in the y-direction of the arrays of the neighboring
        !> buffer layer of type 2
        !
        !>@param bf_i_min
        !> min x-index in copying the grid points from the first
        !> buffer layer
        !
        !>@param bf_j_min
        !> min y-index in copying the grid points from the first
        !> buffer layer
        !
        !>@param nbf_i_min
        !> min x-index in receiving the grid points in the neighbor
        !> buffer layer of type 2
        !
        !>@param nbf_j_min
        !> min y-index in receiving the grid points in the neighbor
        !> buffer layer of type 2
        !
        !>@param bf_copy_size_x
        !> extend of the layer copied in the x-direction
        !
        !>@param bf_copy_size_y
        !> extend of the layer copied in the y-direction
        !--------------------------------------------------------------
        subroutine get_match_indices_for_exchange_with_neighbor2(
     $     bf_alignment,
     $     nbf_alignment,
     $     bf_i_min,
     $     bf_j_min,
     $     nbf_i_min,
     $     nbf_j_min,
     $     bf_copy_size_x,
     $     bf_copy_size_y)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer(ikind)                , intent(out) :: bf_i_min
          integer(ikind)                , intent(out) :: bf_j_min
          integer(ikind)                , intent(out) :: nbf_i_min
          integer(ikind)                , intent(out) :: nbf_j_min
          integer(ikind)                , intent(out) :: bf_copy_size_x
          integer(ikind)                , intent(out) :: bf_copy_size_y


          !we need to define the borders identifying the subarrays
          !copied from the tables of neighbors1 to the tables of the
          !current sublayer

          !for the borders along the x-direction, we first use their
          !general coordinates to identify thme (that is the coordinates
          !with references to the interior domain)

          !then they are expressed as (bf_i_min, bf_i_max), and 
          !(nbf_i_min, nbf_i_max) for the local coordinates of the tables

          !the borders along the y-direction, we use the assumptions:
          ! - for the N and S buffer layers, the neighbor2 corresponds
          !   to the neighbor on the E main layer
          ! - for the E and W buffer layers, the neighbor2 corresponds
          !   to the neighbor on the N main layer
          
          !using the previous assumptions, it is possible to determine the
          !borders along the y-direction as local coordinates for the tables
          !(bf_j_min, bf_j_max) and (nbf_j_min, nbf_j_max)
          call get_exchange_indices(
     $         bf_alignment, nbf_alignment,
     $         x_direction,
     $         bf_i_min, nbf_i_min,
     $         bf_copy_size_x)
          
          call get_exchange_indices(
     $         bf_alignment, nbf_alignment,
     $         y_direction,
     $         bf_j_min, nbf_j_min,
     $         bf_copy_size_y)

c$$$          !get the local min and max borders along the y-direction as
c$$$          !local coordinates
c$$$          select case(localization)
c$$$            case(N)
c$$$               call get_S_exch_indices(bf_j_min)
c$$$               call get_N_exch_indices(nbf_size_y, nbf_j_min)
c$$$            case(S,E,W)
c$$$               call get_N_exch_indices(bf_size_y, bf_j_min)
c$$$               call get_S_exch_indices(nbf_j_min)
c$$$            case default
c$$$               call error_mainlayer_id(
c$$$     $              'bf_layer_exchange_module',
c$$$     $              'get_match_indices_for_exchange_with_neighbor2',
c$$$     $              localization)
c$$$          end select
c$$$
c$$$          bf_copy_size_y = 2*bc_size

        end subroutine get_match_indices_for_exchange_with_neighbor2


c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> get the min border in the y-direction identifying the 
c$$$        !> north layer exchanged
c$$$        !
c$$$        !> @date
c$$$        !> 27_05_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param size_y
c$$$        !> extent in the y-direction of the buffer layer arrays
c$$$        !
c$$$        !>@param bf_j_min
c$$$        !> min border in the y-direction identifying the 
c$$$        !> north layer exchanged
c$$$        !--------------------------------------------------------------
c$$$        subroutine get_N_exch_indices(size_y, bf_j_min)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), intent(in)  :: size_y
c$$$          integer(ikind), intent(out) :: bf_j_min
c$$$
c$$$          bf_j_min = size_y-2*bc_size+1
c$$$
c$$$        end subroutine get_N_exch_indices
c$$$
c$$$
c$$$        !> @author
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> get the min border in the y-direction identifying the 
c$$$        !> south layer exchanged
c$$$        !
c$$$        !> @date
c$$$        !> 27_05_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param bf_j_min
c$$$        !> min border in the y-direction identifying the 
c$$$        !> south layer exchanged
c$$$        !--------------------------------------------------------------
c$$$        subroutine get_S_exch_indices(bf_j_min)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), intent(out) :: bf_j_min
c$$$
c$$$          bf_j_min = 1
c$$$
c$$$        end subroutine get_S_exch_indices


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the buffer layer and its neighbor have gridpoints
        !> in common along the x-direction
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param bf_alignment
        !> alignment of the first buffer layer
        !
        !>@param nbf_alignment
        !> alignment of the neighboring buffer layer
        !
        !>@return overlap
        !> logical indicating an overlap of the two buffer layers
        !--------------------------------------------------------------
        function do_grdpts_overlap_along_x_dir(
     $     bf_alignment, nbf_alignment) result(overlap)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          logical                                     :: overlap

          overlap = (min(bf_alignment(1,2), nbf_alignment(1,2)) -
     $              max(bf_alignment(1,1), nbf_alignment(1,1)) +
     $              2*bc_size + 1).gt.0

        end function do_grdpts_overlap_along_x_dir


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying teh borders of the layer
        !> exchanged in the x-direction
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param bf_alignment
        !> alignment of the first buffer layer
        !
        !>@param nbf_alignment
        !> alignment of the neighboring buffer layer
        !
        !>@param bf_i_min
        !> integer identifying the min x-border of the layer exchanged
        !> for the first buffer layer
        !
        !>@param nbf_i_min
        !> integer identifying the min x-border of the layer exchanged
        !> for the neighboring buffer layer
        !
        !>@param bf_copy_size_x
        !> extent of the layer exchanged in the x-direction
        !--------------------------------------------------------------
        subroutine get_exchange_indices(
     $     bf_alignment, nbf_alignment,dir,
     $     bf_min, nbf_min, bf_copy_size)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer                       , intent(in)  :: dir
          integer(ikind)                , intent(out) :: bf_min
          integer(ikind)                , intent(out) :: nbf_min
          integer(ikind)                , intent(out) :: bf_copy_size

          integer(ikind) :: min_border, max_border


          !get the min and max borders along the x-direction as
          !x-component of the general coordinates
          min_border = max(bf_alignment(dir,1), nbf_alignment(dir,1)) - bc_size
          max_border = min(bf_alignment(dir,2), nbf_alignment(dir,2)) + bc_size

          !convert the previous data into local coordinates for the current
          !buffer layer and the neighbor1
          bf_min  = get_local_coord(dir,min_border, bf_alignment)
          nbf_min = get_local_coord(dir,min_border, nbf_alignment)

          !copy size x
          bf_copy_size = max_border-min_border+1

        end subroutine get_exchange_indices


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the local coordinates along the x-direction
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param i_general_coord
        !> x-general coordinate of the grid point
        !
        !>@param bf_alignment
        !> alignment of the buffer layer
        !
        !>@return i_local_coord
        !> x-local coordinate of the grid point
        !--------------------------------------------------------------
        function get_local_coord(
     $     dir,general_coord, bf_alignment)
     $     result(local_coord)

          implicit none

          integer                       , intent(in) :: dir
          integer(ikind)                , intent(in) :: general_coord
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                             :: local_coord

          local_coord = general_coord - (bf_alignment(dir,1) - bc_size - 1)

        end function get_local_coord              


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying the layers to be exchanged
        !> between a buffer layer and the interior domain
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param localization
        !> cardinal coordinate identifying the position of the
        !> buffer layer
        !
        !>@param bf_alignment
        !> alignment of the buffer layer
        !
        !>@param in_send
        !> x- and y-indices for the SW corner of the table send by
        !> the interior domain
        !
        !>@param in_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the interior domain
        !
        !>@param bf_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer
        !
        !>@param bf_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer
        !
        !>@param ex_size
        !> size-x and size-y of the exchanged arrays
        !--------------------------------------------------------------
        subroutine get_sync_indices_with_interior(
     $     localization,
     $     bf_alignment,
     $     bf_size,
     $     in_send,
     $     in_recv,
     $     bf_send,
     $     bf_recv,
     $     ex_size)

          implicit none

          integer                       , intent(in)  :: localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf_size
          integer(ikind), dimension(2)  , intent(out) :: in_send
          integer(ikind), dimension(2)  , intent(out) :: in_recv
          integer(ikind), dimension(2)  , intent(out) :: bf_send
          integer(ikind), dimension(2)  , intent(out) :: bf_recv
          integer(ikind), dimension(2)  , intent(out) :: ex_size


          select case(localization)

            case(N)

               in_send =
     $             [max(1 , bf_alignment(1,1)-bc_size),
     $              ny-2*bc_size+1]

               in_recv =
     $             [in_send(1),
     $              ny-bc_size+1]

               bf_send = 
     $             [in_send(1) - (bf_alignment(1,1)-(bc_size+1)),
     $              bc_size+1]

               bf_recv = 
     $             [bf_send(1),
     $              1]

               ex_size =
     $             [min(nx, bf_alignment(1,2)+bc_size)-in_send(1)+1,
     $              bc_size]

            case(S)

               in_send =
     $             [max(1 , bf_alignment(1,1)-bc_size),
     $              bc_size+1]

               in_recv =
     $             [in_send(1),
     $              1]

               bf_send = 
     $             [in_send(1) - (bf_alignment(1,1)-(bc_size+1)),
     $              bf_size(2)-2*bc_size+1]

               bf_recv = 
     $             [bf_send(1),
     $              bf_size(2)-bc_size+1]

               ex_size =
     $             [min(nx, bf_alignment(1,2)+bc_size)-in_send(1)+1,
     $              bc_size]

            case(E)
               
               in_send =
     $             [nx-2*bc_size+1,
     $              max(1 , bf_alignment(2,1)-bc_size)]

               in_recv =
     $             [nx-bc_size+1,
     $              in_send(2)]

               bf_send = 
     $             [bc_size+1,
     $              in_send(2) - (bf_alignment(2,1)-(bc_size+1))]

               bf_recv = 
     $             [1,
     $              bf_send(2)]

               ex_size =
     $             [bc_size,
     $              min(ny, bf_alignment(2,2)+bc_size)-in_send(2)+1]

            case(W)

               in_send =
     $             [bc_size+1,
     $              max(1 , bf_alignment(2,1)-bc_size)]

               in_recv =
     $             [1,
     $              in_send(2)]

               bf_send = 
     $             [bf_size(1)-2*bc_size+1,
     $              in_send(2) - (bf_alignment(2,1)-(bc_size+1))]

               bf_recv = 
     $             [bf_size(1)-bc_size+1,
     $              bf_send(2)]

               ex_size =
     $             [bc_size,
     $              min(ny, bf_alignment(2,2)+bc_size)-in_send(2)+1]

            case default
               call error_mainlayer_id(
     $              'bf_layer_exchange_module',
     $              'get_sync_indices_with_interior',
     $              localization)

          end select          

        end subroutine get_sync_indices_with_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> exchange the nodes between two arrays
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_nodes
        !> grid points from the interior domain
        !
        !>@param bf_nodes
        !> grid points from the buffer layer
        !
        !>@param in_send
        !> x- and y-indices for the SW corner of the table send by
        !> the interior domain
        !
        !>@param in_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the interior domain
        !
        !>@param bf_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer
        !
        !>@param bf_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer
        !
        !>@param ex_size
        !> size-x and size-y of the exchanged arrays
        !--------------------------------------------------------------
        subroutine sync_nodes(
     $     bf1_nodes,
     $     bf1_send,
     $     bf1_recv,
     $     bf2_nodes,
     $     bf2_send,
     $     bf2_recv,
     $     ex_size)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(inout) :: bf1_nodes
          integer(ikind), dimension(2)    , intent(in)    :: bf1_send
          integer(ikind), dimension(2)    , intent(in)    :: bf1_recv
          real(rkind)   , dimension(:,:,:), intent(inout) :: bf2_nodes
          integer(ikind), dimension(2)    , intent(in)    :: bf2_send
          integer(ikind), dimension(2)    , intent(in)    :: bf2_recv
          integer(ikind), dimension(2)    , intent(in)    :: ex_size


          integer(ikind) :: i,j
          integer        :: k


          if((ex_size(1).gt.0).and.(ex_size(2).gt.0)) then

             do k=1, ne
                do j=1, ex_size(2)
                   do i=1, ex_size(1)

                      bf2_nodes(bf2_recv(1)+(i-1),bf2_recv(2)+(j-1),k) =
     $                     bf1_nodes(bf1_send(1)+(i-1),bf1_send(2)+(j-1),k)

                      bf1_nodes(bf1_recv(1)+(i-1),bf1_recv(2)+(j-1),k) =
     $                     bf2_nodes(bf2_send(1)+(i-1),bf2_send(2)+(j-1),k)

                   end do
                end do
             end do             

          end if

        end subroutine sync_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying the layers to be exchanged
        !> between a buffer layer and its neighbor1
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param bf1_localization
        !> cardinal coordinate identifying the position of the
        !> buffer layer 1
        !
        !>@param bf1_alignment
        !> alignment of the buffer layer 1
        !
        !>@param bf1_size
        !> size of the nodes table for the buffer layer 1
        !
        !>@param bf1_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer 1
        !
        !>@param bf1_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer 1
        !
        !>@param bf2_alignment
        !> alignment of the buffer layer 2
        !
        !>@param bf2_size
        !> size of the nodes table for the buffer layer 2
        !
        !>@param bf2_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer 2
        !
        !>@param bf2_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer 2
        !
        !>@param ex_size
        !> size-x and size-y of the exchanged arrays
        !--------------------------------------------------------------
        subroutine get_sync_indices_with_neighbor1(
     $     bf1_localization,
     $     bf1_alignment,
     $     bf1_size,
     $     bf1_send,
     $     bf1_recv,
     $     bf2_alignment,
     $     bf2_size,
     $     bf2_send,
     $     bf2_recv,
     $     ex_size)

          implicit none

          integer                       , intent(in)  :: bf1_localization
          integer(ikind), dimension(2,2), intent(in)  :: bf1_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf1_size
          integer(ikind), dimension(2)  , intent(out) :: bf1_send
          integer(ikind), dimension(2)  , intent(out) :: bf1_recv
          integer(ikind), dimension(2,2), intent(in)  :: bf2_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf2_size
          integer(ikind), dimension(2)  , intent(out) :: bf2_send
          integer(ikind), dimension(2)  , intent(out) :: bf2_recv
          integer(ikind), dimension(2)  , intent(out) :: ex_size


          !determination of the x-coordinates for the exchanges
          call get_x_sync_indices(
     $         bf1_alignment,
     $         bf1_send(1), bf1_recv(1),
     $         bf2_alignment,
     $         bf2_send(1), bf2_recv(1),
     $         ex_size)          

          select case(bf1_localization)

            case(N)
                              
               !N is bf1
               bf1_send(2) = bc_size+1
               bf1_recv(2) = 1
               
               !W is bf2
               bf2_send(2) = bf2_size(2)-2*bc_size+1
               bf2_recv(2) = bf2_size(2)-bc_size+1

            case(S)

               !S is bf1
               bf1_send(2) = bf1_size(2)-2*bc_size+1
               bf1_recv(2) = bf1_size(2)-bc_size+1
               
               !W is bf2
               bf2_send(2) = bc_size+1
               bf2_recv(2) = 1

            case(E,W)

               !E is bf1
               bf1_send(2) = bc_size+1
               bf1_recv(2) = 1
               
               !S is bf2
               bf2_send(2) = bf2_size(2)-2*bc_size+1
               bf2_recv(2) = bf2_size(2)-bc_size+1
 
            case default
               call error_mainlayer_id(
     $              'bf_layer_exchange_module',
     $              'get_sync_indices_with_interior',
     $              bf1_localization)
          end select          

        end subroutine get_sync_indices_with_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the indices identifying the layers to be exchanged
        !> between a buffer layer and its neighbor1
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param bf1_localization
        !> cardinal coordinate identifying the position of the
        !> buffer layer 1
        !
        !>@param bf1_alignment
        !> alignment of the buffer layer 1
        !
        !>@param bf1_size
        !> size of the nodes table for the buffer layer 1
        !
        !>@param bf1_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer 1
        !
        !>@param bf1_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer 1
        !
        !>@param bf2_alignment
        !> alignment of the buffer layer 2
        !
        !>@param bf2_size
        !> size of the nodes table for the buffer layer 2
        !
        !>@param bf2_send
        !> x- and y-indices for the SW corner of the table send by
        !> the buffer layer 2
        !
        !>@param bf2_recv
        !> x- and y-indices for the SW corner of the table received
        !> by the buffer layer 2
        !
        !>@param ex_size
        !> size-x and size-y of the exchanged arrays
        !--------------------------------------------------------------
        subroutine get_sync_indices_with_neighbor2(
     $     bf1_localization,
     $     bf1_alignment,
     $     bf1_size,
     $     bf1_send,
     $     bf1_recv,
     $     bf2_alignment,
     $     bf2_size,
     $     bf2_send,
     $     bf2_recv,
     $     ex_size)

          implicit none

          integer                       , intent(in)  :: bf1_localization
          integer(ikind), dimension(2,2), intent(in)  :: bf1_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf1_size
          integer(ikind), dimension(2)  , intent(out) :: bf1_send
          integer(ikind), dimension(2)  , intent(out) :: bf1_recv
          integer(ikind), dimension(2,2), intent(in)  :: bf2_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf2_size
          integer(ikind), dimension(2)  , intent(out) :: bf2_send
          integer(ikind), dimension(2)  , intent(out) :: bf2_recv
          integer(ikind), dimension(2)  , intent(out) :: ex_size


          !determination of the x-coordinates for the exchanges
          call get_x_sync_indices(
     $         bf1_alignment,
     $         bf1_send(1), bf1_recv(1),
     $         bf2_alignment,
     $         bf2_send(1), bf2_recv(1),
     $         ex_size)          

          select case(bf1_localization)

            case(N)
                              
               !N is bf1
               bf1_send(2) = bc_size+1
               bf1_recv(2) = 1
               
               !E is bf2
               bf2_send(2) = bf2_size(2)-2*bc_size+1
               bf2_recv(2) = bf2_size(2)-bc_size+1

            case(S)

               !S is bf1
               bf1_send(2) = bf1_size(2)-2*bc_size+1
               bf1_recv(2) = bf1_size(2)-bc_size+1
               
               !E is bf2
               bf2_send(2) = bc_size+1
               bf2_recv(2) = 1

            case(E,W)

               !E,W is bf1
               bf1_send(2) = bf1_size(2)-2*bc_size+1
               bf1_recv(2) = bf1_size(2)-bc_size+1
               
               !N is bf2
               bf2_send(2) = bc_size+1
               bf2_recv(2) = 1
 
            case default
               call error_mainlayer_id(
     $              'bf_layer_exchange_module',
     $              'get_sync_indices_with_interior',
     $              bf1_localization)
          end select          

        end subroutine get_sync_indices_with_neighbor2


        subroutine get_x_sync_indices(
     $     bf1_alignment,
     $     bf1_send_x,
     $     bf1_recv_x,
     $     bf2_alignment,
     $     bf2_send_x,
     $     bf2_recv_x,
     $     ex_size)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf1_alignment
          integer(ikind)                , intent(out) :: bf1_send_x
          integer(ikind)                , intent(out) :: bf1_recv_x
          integer(ikind), dimension(2,2), intent(in)  :: bf2_alignment
          integer(ikind)                , intent(out) :: bf2_send_x
          integer(ikind)                , intent(out) :: bf2_recv_x
          integer(ikind), dimension(2)  , intent(out) :: ex_size

          integer(ikind), dimension(2) :: gen_coords

          gen_coords(1) = max(bf1_alignment(1,1),bf2_alignment(1,1))-bc_size
          gen_coords(2) = min(bf1_alignment(1,2),bf2_alignment(1,2))+bc_size

          bf1_send_x = gen_coords(1) - (bf1_alignment(1,1)-(bc_size+1))
          bf1_recv_x = bf1_send_x
          bf2_send_x = gen_coords(1) - (bf2_alignment(1,1)-(bc_size+1))
          bf2_recv_x = bf2_send_x

          ex_size(1) = gen_coords(2)-gen_coords(1)+1
          ex_size(2) = bc_size

        end subroutine get_x_sync_indices


        subroutine update_alignment_for_exchanges(
     $     bf_mainlayer_id,
     $     bf_alignment)

          implicit none

          integer                       , intent(in)    :: bf_mainlayer_id
          integer(ikind), dimension(2,2), intent(inout) :: bf_alignment


          select case(bf_mainlayer_id)

            case(N,S)
               if(  ((bf_alignment(1,1)-bc_size).le.(align_W+bc_size)).and.
     $              ((bf_alignment(1,1)-bc_size).ge.(align_W))) then
                  bf_alignment(1,1) = align_W+1
               end if

               if(  ((bf_alignment(1,2)+bc_size).ge.(align_E-bc_size)).and.
     $              ((bf_alignment(1,2)+bc_size).le.(align_E))) then
                  bf_alignment(1,2) = align_E-1
               end if


            case(E,W)
               if(  ((bf_alignment(2,1)-bc_size).le.(align_S+bc_size)).and.
     $              ((bf_alignment(2,1)-bc_size).ge.(align_S))) then
                  bf_alignment(2,1) = align_S+1
               end if

               if(  ((bf_alignment(2,2)+bc_size).ge.(align_N-bc_size)).and.
     $              ((bf_alignment(2,2)+bc_size).le.(align_N))) then
                  bf_alignment(2,2) = align_N-1
               end if


            case default

               call error_mainlayer_id(
     $              'bf_layer_exchange_module',
     $              'update_alignment_for_exchanges',
     $              bf_mainlayer_id)

          end select

        end subroutine update_alignment_for_exchanges

      end module bf_layer_exchange_module
