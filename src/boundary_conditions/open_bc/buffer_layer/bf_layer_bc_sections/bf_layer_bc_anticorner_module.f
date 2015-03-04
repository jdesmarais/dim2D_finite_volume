      !> @file
      !> annex subroutines useful when computing the anti-corner
      !> boundary sections in the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to extract grid points
      !> from the interior nodes to complete the grid points needed
      !> when computing fluxes for anti-corners
      !
      !> @date
      ! 04_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_bc_anticorner_module

        use bf_layer_extract_module, only :
     $       get_indices_to_extract_interior_data,
     $       get_indices_to_extract_bf_layer_data,
     $       get_grdpts_id_from_interior

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none

        private
        public :: 
     $       get_coords_from_pattern,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       extract_grdpts_to_compute_anticorner_fluxes

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether gridpoints are needed around the 
        !> central gridpoint to compute the fluxes in the 
        !> x-direction
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param p_model
        !> physical model
        !
        !> @param operator_type
        !> integer identifying the type of space discretization
        !> operator used
        !
        !> @param i
        !> x-index of the central grid point
        !
        !> @param j
        !> y-index of the central grid point
        !
        !> @param size_x
        !> x-extent of the array used to potentially compute the x-fluxes
        !
        !> @param size_y
        !> y-extent of the array used to potentially compute the x-fluxes
        !
        !> @return border_coords
        !> indices identifying the SW and NE corners of the grid points
        !> needed to compute the fluxes
        !
        !> @return cpt_coords
        !> indices identifying the position of the central gridpoint
        !> in the temporary array of gridpoints created around
        !
        !> @return grdpts_needed
        !> logical identifying whether additional gridpoints are needed
        !> to potentially compute the x-fluxes
        !--------------------------------------------------------------
        function are_grdpts_needed_for_flux_x(
     $     p_model,
     $     operator_type,
     $     i,j,
     $     size_x,size_y,
     $     border_coords,
     $     cpt_coords)
     $     result(grdpts_needed)

          implicit none

          type(pmodel_eq)               , intent(in)  :: p_model
          integer                       , intent(in)  :: operator_type
          integer(ikind)                , intent(in)  :: i
          integer(ikind)                , intent(in)  :: j
          integer(ikind)                , intent(in)  :: size_x
          integer(ikind)                , intent(in)  :: size_y
          integer(ikind), dimension(2,2), intent(out) :: border_coords
          integer(ikind), dimension(2)  , intent(out) :: cpt_coords
          logical                                     :: grdpts_needed
          
          integer, dimension(2,2) :: pattern

          pattern = p_model%get_sd_pattern_flux_x(operator_type)

          call get_coords_from_pattern(
     $         pattern,i,j,size_x,size_y,
     $         border_coords,
     $         cpt_coords,
     $         grdpts_needed)          

        end function are_grdpts_needed_for_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether gridpoints are needed around the 
        !> central gridpoint to compute the fluxes in the 
        !> y-direction
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param p_model
        !> physical model
        !
        !> @param operator_type
        !> integer identifying the type of space discretization
        !> operator used
        !
        !> @param i
        !> x-index of the central grid point
        !
        !> @param j
        !> y-index of the central grid point
        !
        !> @param size_x
        !> x-extent of the array used to potentially compute the y-fluxes
        !
        !> @param size_y
        !> y-extent of the array used to potentially compute the y-fluxes
        !
        !> @return border_coords
        !> indices identifying the SW and NE corners of the grid points
        !> needed to compute the fluxes
        !
        !> @return cpt_coords
        !> indices identifying the position of the central gridpoint
        !> in the temporary array of gridpoints created around
        !
        !> @return grdpts_needed
        !> logical identifying whether additional gridpoints are needed
        !> to potentially compute the y-fluxes
        !--------------------------------------------------------------
        function are_grdpts_needed_for_flux_y(
     $     p_model,
     $     operator_type,
     $     i,j,
     $     size_x,size_y,
     $     border_coords,
     $     cpt_coords)
     $     result(grdpts_needed)

          implicit none

          type(pmodel_eq)               , intent(in)  :: p_model
          integer                       , intent(in)  :: operator_type
          integer(ikind)                , intent(in)  :: i
          integer(ikind)                , intent(in)  :: j
          integer(ikind)                , intent(in)  :: size_x
          integer(ikind)                , intent(in)  :: size_y
          integer(ikind), dimension(2,2), intent(out) :: border_coords
          integer(ikind), dimension(2)  , intent(out) :: cpt_coords
          logical                                     :: grdpts_needed
          
          integer, dimension(2,2) :: pattern

          pattern = p_model%get_sd_pattern_flux_y(operator_type)

          call get_coords_from_pattern(
     $         pattern,i,j,size_x,size_y,
     $         border_coords,
     $         cpt_coords,
     $         grdpts_needed)

        end function are_grdpts_needed_for_flux_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the coords of the borders of the grid points needed
        !> as well as the coords of the central gridpoint and if
        !> grid points are needed
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param pattern
        !> integers identifying the extent of the grid points needed
        !> around the central grid point
        !
        !> @param i
        !> x-index of the central grid point
        !
        !> @param j
        !> y-index of the central grid point
        !
        !> @param size_x
        !> x-extent of the array used to potentially compute the y-fluxes
        !
        !> @param size_y
        !> y-extent of the array used to potentially compute the y-fluxes
        !
        !> @return border_coords
        !> indices identifying the SW and NE corners of the grid points
        !> needed to compute the fluxes
        !
        !> @return cpt_coords
        !> indices identifying the position of the central gridpoint
        !> in the temporary array of gridpoints created around
        !
        !> @return grdpts_needed
        !> logical identifying whether additional gridpoints are needed
        !> to potentially compute the y-fluxes
        !--------------------------------------------------------------        
        subroutine get_coords_from_pattern(
     $     pattern, i,j, size_x, size_y,
     $     border_coords,
     $     cpt_coords,
     $     grdpts_needed)

          implicit none

          integer       , dimension(2,2), intent(in)  :: pattern
          integer(ikind)                , intent(in)  :: i
          integer(ikind)                , intent(in)  :: j
          integer(ikind)                , intent(in)  :: size_x
          integer(ikind)                , intent(in)  :: size_y
          integer(ikind), dimension(2,2), intent(out) :: border_coords
          integer(ikind), dimension(2)  , intent(out) :: cpt_coords
          logical                                     :: grdpts_needed

          border_coords(1,1) = i+pattern(1,1)
          border_coords(1,2) = i+pattern(1,2)
          border_coords(2,1) = j+pattern(2,1)
          border_coords(2,2) = j+pattern(2,2)

          cpt_coords(1)   = -pattern(1,1)+1
          cpt_coords(2)   = -pattern(2,1)+1

          grdpts_needed =
     $         (border_coords(1,1).lt.1).or.
     $         (border_coords(1,2).gt.size_x).or.
     $         (border_coords(2,1).lt.1).or.
     $         (border_coords(2,2).gt.size_y)

        end subroutine get_coords_from_pattern


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> gather grid points from the interior and the buffer layer
        !> to have all the grid points needed to compute the fluxes
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !
        !> @param bf_nodes
        !> buffer layer grid points
        !
        !> @param interior_nodes
        !> grid points of the interior domain
        !
        !> @param gen_coords
        !> coordinates of the four borders of the tmp array expressed as
        !> SW_corner = [gen_coords(1,1),gen_coords(2,1)]
        !> NE_corner = [gen_coords(1,2),gen_coords(2,2)]
        !> the coordinates are expressed in the general coordinate reference
        !> frame
        !
        !> @return tmp_nodes
        !> temporary array containing the nodes needed to compute 
        !> the fluxes
        !--------------------------------------------------------------        
        subroutine extract_grdpts_to_compute_anticorner_fluxes(
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_nodes,
     $     interior_nodes,
     $     gen_coords,
     $     tmp_nodes)

          implicit none
          
          integer(ikind), dimension(2,2)     , intent(in)  :: bf_alignment
          integer       , dimension(:,:)     , intent(in)  :: bf_grdpts_id
          real(rkind)   , dimension(:,:,:)   , intent(in)  :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)  :: gen_coords
          real(rkind)   , dimension(:,:,:)   , intent(out) :: tmp_nodes

          
          integer(ikind)                       :: size_x,size_y
          integer(ikind)                       :: i_recv,i_send,j_recv,j_send
          integer(ikind)                       :: i,j
          integer                              :: k
          integer, dimension(:,:), allocatable :: tmp_grdpts_id
          logical                              :: grdpts_available


          ! synchronize the nodes of the interior with the tmp nodes
          call get_indices_to_extract_interior_data(
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)

          do k=1,ne
             do j=1, size_y
                do i=1, size_x
                   
                   tmp_nodes(i_recv+i-1,j_recv+j-1,k) = 
     $                  interior_nodes(i_send+i-1,j_send+j-1,k)

                end do
             end do
          end do

          ! synchronize the nodes of the buffer layer with the tmp_nodes
          call get_indices_to_extract_bf_layer_data(
     $         bf_alignment,
     $         gen_coords,
     $         size_x, size_y,
     $         i_recv, j_recv,
     $         i_send, j_send)

          do k=1,ne
             do j=1, size_y
                do i=1, size_x
                   
                   tmp_nodes(i_recv+i-1,j_recv+j-1,k) = 
     $                  bf_nodes(i_send+i-1,j_send+j-1,k)

                end do
             end do
          end do

          ! determine the grdpts_id for the extraction of the grid-points
          allocate(tmp_grdpts_id(
     $         gen_coords(1,2)-gen_coords(1,1)+1,
     $         gen_coords(2,2)-gen_coords(2,1)+1))

          ! from the interior
          call get_grdpts_id_from_interior(tmp_grdpts_id,gen_coords)

          ! from the buffer layer
          do j=1, size_y
             do i=1, size_x
                
                tmp_grdpts_id(i_recv+i-1,j_recv+j-1) = 
     $               bf_grdpts_id(i_send+i-1,j_send+j-1)

             end do
          end do

          ! check if all grid points needed to compute the
          ! fluxes are available
          grdpts_available = .true.
          do j=1, size(tmp_grdpts_id,2)
             do i=1,size(tmp_grdpts_id,1)
                if(tmp_grdpts_id(i,j).eq.no_pt) then
                   grdpts_available = .false.
                end if
             end do
          end do

          if(.not.grdpts_available) then
             print '(''bf_layer_bc_anticorner_module'')'
             print '(''extract_grdpts_to_compute_anticorner_fluxes'')'
             print '(''not all grid points are available'')'
             print '(''grdpts_id: '')'
             print *, tmp_grdpts_id
             stop ''
          end if          

          deallocate(tmp_grdpts_id)

        end subroutine extract_grdpts_to_compute_anticorner_fluxes

      end module bf_layer_bc_anticorner_module
      
