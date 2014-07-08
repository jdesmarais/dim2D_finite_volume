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

        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : ne, bc_size
        use parameters_kind    , only : ikind, rkind
      
        implicit none

        private
        public :: do_grdpts_overlap_along_x_dir,
     $            get_match_indices_for_exchange_with_neighbor1,
     $            get_match_indices_for_exchange_with_neighbor2,
     $            copy_from_bf1_to_bf2

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

          
          do k=1, ne
             do j=bf1_j_min, bf1_j_min+bf_copy_size_y-1
                do i=bf1_i_min, bf1_i_min+bf_copy_size_x-1
                   bf2_nodes(i-bf1_i_min+bf2_i_min,
     $                       j-bf1_j_min+bf2_j_min,
     $                       k) = bf1_nodes(i,j,k)
                end do
             end do
          end do

          do j=bf1_j_min, bf1_j_min+bf_copy_size_y-1
             do i=bf1_i_min, bf1_i_min+bf_copy_size_x-1
                bf2_grdpts_id(
     $               i-bf1_i_min+bf2_i_min,
     $               j-bf1_j_min+bf2_j_min) = bf1_grdpts_id(i,j)
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
     $     localization,
     $     bf_alignment, bf_size_y,
     $     nbf_alignment, nbf_size_y,
     $     bf_i_min, bf_j_min,
     $     nbf_i_min, nbf_j_min,
     $     bf_copy_size_x, bf_copy_size_y)

          implicit none

          integer                       , intent(in)  :: localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind)                , intent(in)  :: bf_size_y
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer(ikind)                , intent(in)  :: nbf_size_y
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
          ! - for the N and S buffer layers, the neighbor1 corresponds
          !   to the neighbor on the W main layer
          ! - for the E and W buffer layers, the neighbor1 corresponds
          !   to the neighbor on the S main layer
          
          !using the previous assumptions, it is possible to determine the
          !borders along the y-direction as local coordinates for the tables
          !(bf_j_min, bf_j_max) and (nbf_j_min, nbf_j_max)
          call get_x_exchange_indices(bf_alignment, nbf_alignment,
     $                                bf_i_min, nbf_i_min,
     $                                bf_copy_size_x)
          
          !get the local min and max borders along the y-direction as
          !local coordinates
          select case(localization)
            case(N,E,W)
               call get_S_exch_indices(bf_j_min)
               call get_N_exch_indices(nbf_size_y, nbf_j_min)
            case(S)
               call get_N_exch_indices(bf_size_y, bf_j_min)
               call get_S_exch_indices(nbf_j_min)
            case default
               print '(''bf_layer_class'')'
               print '(''copy_from_neighbor1'')'
               print '(''localization not recognized'')'
               print '(''localization: '', I2)', localization
               stop 'change loclaization'
          end select

          bf_copy_size_y = 2*bc_size

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
     $     localization,
     $     bf_alignment, bf_size_y,
     $     nbf_alignment, nbf_size_y,
     $     bf_i_min, bf_j_min,
     $     nbf_i_min, nbf_j_min,
     $     bf_copy_size_x, bf_copy_size_y)

          implicit none

          integer                       , intent(in)  :: localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind)                , intent(in)  :: bf_size_y
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer(ikind)                , intent(in)  :: nbf_size_y
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
          call get_x_exchange_indices(bf_alignment, nbf_alignment,
     $                                bf_i_min, nbf_i_min,
     $                                bf_copy_size_x)
          
          !get the local min and max borders along the y-direction as
          !local coordinates
          select case(localization)
            case(N)
               call get_S_exch_indices(bf_j_min)
               call get_N_exch_indices(nbf_size_y, nbf_j_min)
            case(S,E,W)
               call get_N_exch_indices(bf_size_y, bf_j_min)
               call get_S_exch_indices(nbf_j_min)
            case default
               print '(''bf_layer_class'')'
               print '(''copy_from_neighbor1'')'
               print '(''localization not recognized'')'
               print '(''localization: '', I2)', localization
               stop 'change loclaization'
          end select

          bf_copy_size_y = 2*bc_size

        end subroutine get_match_indices_for_exchange_with_neighbor2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the min border in the y-direction identifying the 
        !> north layer exchanged
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param size_y
        !> extent in the y-direction of the buffer layer arrays
        !
        !>@param bf_j_min
        !> min border in the y-direction identifying the 
        !> north layer exchanged
        !--------------------------------------------------------------
        subroutine get_N_exch_indices(size_y, bf_j_min)

          implicit none

          integer(ikind), intent(in)  :: size_y
          integer(ikind), intent(out) :: bf_j_min

          bf_j_min = size_y-2*bc_size+1

        end subroutine get_N_exch_indices


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the min border in the y-direction identifying the 
        !> south layer exchanged
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param bf_j_min
        !> min border in the y-direction identifying the 
        !> south layer exchanged
        !--------------------------------------------------------------
        subroutine get_S_exch_indices(bf_j_min)

          implicit none

          integer(ikind), intent(out) :: bf_j_min

          bf_j_min = 1

        end subroutine get_S_exch_indices


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
        subroutine get_x_exchange_indices(
     $     bf_alignment, nbf_alignment,
     $     bf_i_min, nbf_i_min, bf_copy_size_x)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
          integer(ikind)                , intent(out) :: bf_i_min
          integer(ikind)                , intent(out) :: nbf_i_min
          integer(ikind)                , intent(out) :: bf_copy_size_x

          integer(ikind) :: min_border, max_border


          !get the min and max borders along the x-direction as
          !x-component of the general coordinates
          min_border = max(bf_alignment(1,1), nbf_alignment(1,1)) - bc_size
          max_border = min(bf_alignment(1,2), nbf_alignment(1,2)) + bc_size

          !convert the previous data into local coordinates for the current
          !buffer layer and the neighbor1
          bf_i_min  = get_x_local_coord(min_border, bf_alignment)
          nbf_i_min = get_x_local_coord(min_border, nbf_alignment)

          !copy size x
          bf_copy_size_x = max_border-min_border+1

        end subroutine get_x_exchange_indices


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
        function get_x_local_coord(
     $     i_general_coord, bf_alignment)
     $     result(i_local_coord)

          implicit none

          integer(ikind)                , intent(in) :: i_general_coord
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                             :: i_local_coord

          i_local_coord = i_general_coord - (bf_alignment(1,1) - (bc_size+1))        

        end function get_x_local_coord


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the local coordinates along the y-direction
        !
        !> @date
        !> 27_05_2014 - initial version - J.L. Desmarais
        !
        !>@param j_general_coord
        !> y-general coordinate of the grid point
        !
        !>@param bf_alignment
        !> alignment of the buffer layer
        !
        !>@return j_local_coord
        !> y-local coordinates of the grid point
        !--------------------------------------------------------------
        function get_y_local_coord(
     $     j_general_coord, bf_alignment)
     $     result(j_local_coord)

          implicit none

          integer(ikind)                , intent(in) :: j_general_coord
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                             :: j_local_coord

          j_local_coord = j_general_coord - (bf_alignment(2,1) - (bc_size+1))

        end function get_y_local_coord        

      end module bf_layer_exchange_module