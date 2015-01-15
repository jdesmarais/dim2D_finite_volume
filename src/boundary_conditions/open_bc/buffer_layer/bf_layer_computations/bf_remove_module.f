      !> @file
      !> module encapsulating the subroutines needed to check whether
      !> a bf_layer object should be removed
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines needed to check whether
      !> a bf_layer object should be removed
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      ! 17_07_2014 - openbc_undermined in p_model - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_remove_module

        use bf_layer_errors_module, only :
     $     error_mainlayer_id

        use parameters_bf_layer, only :
     $       align_N,align_S,
     $       align_E,align_W,
     $       no_pt,
     $       search_dcr

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq
        
        implicit none

        private
        public :: check_if_bf_layer_remains,
     $            get_check_line_param


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the grid points of the buffer layer common
        !> with the interior domain are such that the buffer layer
        !> can be removed without considering the neighbor dependencies
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param bf_localization
        !> cardinal coordinate identifying the position of the buffer
        !> layer
        !
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !
        !> @param bf_match_table
        !> table allowing to convert general coordinates into local
        !> coordinates for the buffer layer
        !
        !> @param bf_grdpts_id
        !> role of the grid points in the buffer layer
        !
        !> @param bf_nodes
        !> grid points of the the buffer layer
        !
        !> @param interior_nodes
        !> grid points of the interior domain
        !
        !> @return bf_remains
        !> logical identifying whether the local removal is approved or
        !> not
        !---------------------------------------------------------------
        function check_if_bf_layer_remains(
     $       bf_localization,
     $       bf_alignment,
     $       bf_match_table,
     $       bf_grdpts_id,
     $       bf_x_map,
     $       bf_y_map,
     $       bf_nodes,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model)
     $       result(bf_remains)

          implicit none

          integer                          , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2)   , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)     , intent(in)  :: bf_match_table
          integer    , dimension(:,:)      , intent(in)  :: bf_grdpts_id
          real(rkind), dimension(:)        , intent(in)  :: bf_x_map
          real(rkind), dimension(:)        , intent(in)  :: bf_y_map
          real(rkind), dimension(:,:,:)    , intent(in)  :: bf_nodes
          real(rkind), dimension(nx)       , intent(in)  :: interior_x_map
          real(rkind), dimension(ny)       , intent(in)  :: interior_y_map
          real(rkind), dimension(nx,ny,ne) , intent(in)  :: interior_nodes
          type(pmodel_eq)                  , intent(in)  :: p_model
          logical                                        :: bf_remains

          
          integer(ikind), dimension(2,2) :: bf_coords
          integer(ikind), dimension(2,2) :: in_coords
          

          !check if the buffer layer has grid points in common with
          !the interior domain
          if(grdpts_common_with_check_layer(bf_alignment)) then
             
             !depending on the buffer layer position, determine
             !the line around which the grid points will be checked
             !to see whether the region undermines the open boundary
             !conditions or not
             call get_check_line_param(
     $            bf_localization,
     $            bf_alignment,
     $            bf_match_table,
     $            size(bf_nodes,1),
     $            size(bf_nodes,2),
     $            bf_coords,
     $            in_coords)

             !check the neighboring points around the line
             call check_line_neighbors(
     $            bf_coords, in_coords,
     $            bf_grdpts_id,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model,
     $            bf_remains)
          
          !if the buffer layer has no grid point in common with the
          !interior domain, it can be removed immediately (the
          !neighboring buffer layers depending on this buffer layer
          !are not taken into account)
          else
             bf_remains = .false.
          end if

        end function check_if_bf_layer_remains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the alignment of the buffer layer compared to
        !> the interior is such that the buffer layer has grid points
        !> in common with the interior domain
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !
        !> @return grdpts_common
        !> logical stating whether the buffer layer has grid points in
        !> common with the interior domain
        !---------------------------------------------------------------
        function grdpts_common_with_check_layer(bf_alignment)
     $     result(grdpts_common)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          logical                                    :: grdpts_common

          integer(ikind) :: common_layer_size_x
          integer(ikind) :: common_layer_size_y


          common_layer_size_x = min(nx+search_dcr,bf_alignment(1,2)+bc_size)-
     $                          max(1-search_dcr, bf_alignment(1,1)-bc_size)+1

          common_layer_size_y = min(ny+search_dcr,bf_alignment(2,2)+bc_size)-
     $                          max(1-search_dcr, bf_alignment(2,1)-bc_size)+1

          grdpts_common = (common_layer_size_x.gt.0).and.
     $                    (common_layer_size_y.gt.0)


        end function grdpts_common_with_check_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the parameters constraining the line around
        !> which the neighboring points will be checked to see
        !> whether the open boundary conditions are undermined
        !> in this region or not
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param bf_localization
        !> cardinal coordinate identifying the position of the buffer
        !> layer
        !
        !> @param bf_alignment
        !> relative position of the buffer layer to the interior domain
        !
        !> @param bf_match_table
        !> table allowing to convert general coordinates into local
        !> coordinates for the buffer layer
        !
        !> @param bf_size_x
        !> extent of the nodes and grdpts arrays of the buffer layer
        !> in the x-direction
        !
        !> @param bf_size_y
        !> extent of the nodes and grdpts arrays of the buffer layer
        !> in the y-direction
        !
        !> @param bf_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the buffer layer grid points
        !
        !> @param in_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the interior domain grid points
        !---------------------------------------------------------------
        subroutine get_check_line_param(
     $     bf_localization,
     $     bf_alignment, bf_match_table,
     $     bf_size_x, bf_size_y,
     $     bf_coords, in_coords)


          implicit none

          integer                       , intent(in)  :: bf_localization
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2)  , intent(in)  :: bf_match_table
          integer(ikind)                , intent(in)  :: bf_size_x
          integer(ikind)                , intent(in)  :: bf_size_y
          integer(ikind), dimension(2,2), intent(out) :: bf_coords
          integer(ikind), dimension(2,2), intent(out) :: in_coords
          

          select case(bf_localization)

            case(N)
               bf_coords(1,1) = max(2, align_W - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(2, align_N - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_size_x-1, align_E + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_size_y-1, align_N + search_dcr - bf_match_table(2))


               in_coords(1,1) = max(2, bf_alignment(1,1) - bc_size)
               in_coords(2,1) = max(2, align_N - search_dcr)

               in_coords(1,2) = min(nx-1, bf_alignment(1,2) + bc_size         )
               in_coords(2,2) = min(ny-1, bf_coords(2,1)-1 + bf_match_table(2))


            case(S)
               bf_coords(1,1) = max(2, align_W - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(2, align_S - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_size_x-1, align_E + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_size_y-1, align_S + search_dcr - bf_match_table(2))


               in_coords(1,1) = max(2, bf_alignment(1,1) - bc_size)
               in_coords(2,1) = max(2, bf_coords(2,2)+1 + bf_match_table(2))

               in_coords(1,2) = min(nx-1, bf_alignment(1,2) + bc_size)
               in_coords(2,2) = min(ny-1, align_S + search_dcr)


            case(E)
               bf_coords(1,1) = max(2, align_E - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(2, align_S - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_size_x-1, align_E + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_size_y-1, align_N + search_dcr - bf_match_table(2))


               in_coords(1,1) = max(2, align_E - search_dcr)
               in_coords(2,1) = max(2, bf_alignment(2,1) - bc_size)

               in_coords(1,2) = min(nx-1, bf_coords(1,1)-1 + bf_match_table(1))
               in_coords(2,2) = min(ny-1, bf_alignment(2,2) + bc_size)


            case(W)
               bf_coords(1,1) = max(2, align_W - search_dcr - bf_match_table(1))
               bf_coords(2,1) = max(2, align_S - search_dcr - bf_match_table(2))

               bf_coords(1,2) = min(bf_size_x-1, align_W + search_dcr - bf_match_table(1))
               bf_coords(2,2) = min(bf_size_y-1, align_N + search_dcr - bf_match_table(2))


               in_coords(1,1) = max(2, bf_coords(1,2)+1 + bf_match_table(1))
               in_coords(2,1) = max(2, bf_alignment(2,1) - bc_size)

               in_coords(1,2) = min(nx-1, align_W + search_dcr)
               in_coords(2,2) = min(ny-1, bf_alignment(2,2) + bc_size)


            case default
               call error_mainlayer_id(
     $              'bf_layer_remove_module.f',
     $              'get_check_line_param',
     $              bf_localization)

          end select
          
        end subroutine get_check_line_param


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check the grid points in the layer determined previously
        !> by get_check_line_param
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param bf_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the buffer layer grid points
        !
        !> @param in_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the interior domain grid points
        !
        !> @param bf_grdpts_id
        !> role of the grid points in the buffer layer
        !
        !> @param bf_nodes
        !> grid points of the the buffer layer
        !
        !> @param interior_nodes
        !> grid points of the interior domain
        !
        !> @param p_model
        !> physical model telling when the open_bc are undermined
        !
        !> @return bf_remains
        !> logical identifying whether the local removal is approved or
        !> not
        !---------------------------------------------------------------
        subroutine check_line_neighbors(
     $     bf_coords, in_coords,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     bf_remains)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_coords
          integer(ikind), dimension(2,2)     , intent(in)    :: in_coords
          integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                    , intent(in)    :: p_model
          logical                            , intent(out)   :: bf_remains


          bf_remains = .false.


          !check the interior points
          call check_layer_interior(
     $         in_coords,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         bf_remains)
          
          
          !check the buffer layer points
          if(.not.bf_remains) then
             call check_layer_bf(
     $            bf_coords,
     $            bf_grdpts_id,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            bf_remains)
          end if


        end subroutine check_line_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the layer of grid points in the interior domain
        !> undermine the open boundary conditions
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param pt_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the interior domain grid points
        !
        !> @param interior_nodes
        !> grid points of the interior domain
        !
        !> @return bf_remains
        !> logical identifying whether the local removal is approved or
        !> not
        !---------------------------------------------------------------
        subroutine check_layer_interior(
     $     pt_coords,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     bf_remains)
        
          implicit none
          
          integer(ikind), dimension(2,2)     , intent(in)    :: pt_coords
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                    , intent(in)    :: p_model
          logical                            , intent(inout) :: bf_remains
          
          
          integer(ikind) :: i,j


          if( ((pt_coords(1,2)-pt_coords(1,1)).gt.0).and.
     $        ((pt_coords(2,2)-pt_coords(2,1)).gt.0) ) then

             if(
     $            (pt_coords(2,1).ge.2).and.(pt_coords(2,1).le.(ny-1)).and.
     $            (pt_coords(2,2).ge.2).and.(pt_coords(2,2).le.(ny-1)).and.
     $            (pt_coords(1,1).ge.2).and.(pt_coords(1,1).le.(nx-1)).and.
     $            (pt_coords(1,2).ge.2).and.(pt_coords(1,2).le.(nx-1))) then
                
                do j=pt_coords(2,1), pt_coords(2,2)
                   do i=pt_coords(1,1), pt_coords(1,2)
                      
                      bf_remains = p_model%are_openbc_undermined(
     $                     interior_x_map(i-1:i+1),
     $                     interior_y_map(j-1:j+1),
     $                     interior_nodes(i-1:i+1,j-1:j+1,:))
                      
                      if(bf_remains) then
                         exit
                      end if
                      
                   end do
                   
                   if(bf_remains) then
                      exit
                   end if
                   
                end do
                
             else
                
                print '(''bf_remove_module'')'
                print '(''check_layer_interior'')'
                print '(''the nodes that should be checked are'')'
                print '(''not all in the interior domain'')'
                print '(''pt_coords(1,1): '',F6.3)', pt_coords(1,1)
                print '(''pt_coords(1,2): '',F6.3)', pt_coords(1,2)
                print '(''pt_coords(2,1): '',F6.3)', pt_coords(2,1)
                print '(''pt_coords(2,2): '',F6.3)', pt_coords(2,2)
                stop ''
                
             end if

          end if

        end subroutine check_layer_interior

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check if the layer of grid points in the buffer layer 
        !> undermine the open boundary conditions
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !> @param pt_coords
        !> integer table identifying the borders of the layer to be
        !> checked in the interior domain grid points
        !
        !> @param grdpts_id
        !> role of the grid points in the buffer layer
        !
        !> @param nodes
        !> grid points of the buffer layer
        !
        !> @return bf_remains
        !> logical identifying whether the local removal is approved or
        !> not
        !---------------------------------------------------------------
        subroutine check_layer_bf(
     $     pt_coords,
     $     grdpts_id,
     $     x_map,
     $     y_map,
     $     nodes,
     $     p_model,
     $     bf_remains)
        
          implicit none
          
          integer(ikind), dimension(2,2)     , intent(in)   :: pt_coords
          integer       , dimension(:,:)     , intent(in)   :: grdpts_id
          real(rkind)   , dimension(:)       , intent(in)   :: x_map
          real(rkind)   , dimension(:)       , intent(in)   :: y_map
          real(rkind)   , dimension(:,:,:)   , intent(in)   :: nodes
          type(pmodel_eq)                    , intent(in)   :: p_model
          logical                            , intent(inout):: bf_remains
          
          
          integer(ikind) :: i,j
          integer(ikind) :: size_x,size_y

          
          if(  ((pt_coords(1,2)-pt_coords(1,1)).gt.0).and.
     $         ((pt_coords(2,2)-pt_coords(2,1)).gt.0) ) then

             size_x = size(x_map,1)
             size_y = size(y_map,1)

             if(
     $            (pt_coords(2,1).ge.2).and.(pt_coords(2,1).le.(size_y-1)).and.
     $            (pt_coords(2,2).ge.2).and.(pt_coords(2,2).le.(size_y-1)).and.
     $            (pt_coords(1,1).ge.2).and.(pt_coords(1,1).le.(size_x-1)).and.
     $            (pt_coords(1,2).ge.2).and.(pt_coords(1,2).le.(size_x-1))) then

                do j=pt_coords(2,1), pt_coords(2,2)
                   do i=pt_coords(1,1), pt_coords(1,2)
                      
                      if(no_nogrdpt(grdpts_id,i,j)) then
                         
                         bf_remains = p_model%are_openbc_undermined(
     $                        x_map(i-1:i+1),
     $                        y_map(j-1:j+1),
     $                        nodes(i-1:i+1,j-1:j+1,:))
                         
                         if(bf_remains) then
                            exit
                         end if
                      end if
                      
                   end do
                
                   if(bf_remains) then
                      exit
                   end if
                   
                end do

             else
                
                print '(''bf_remove_module'')'
                print '(''check_layer_interior'')'
                print '(''the nodes that should be checked are'')'
                print '(''not all in the buffer layer'')'
                print '(''pt_coords(1,1): '',F6.3)', pt_coords(1,1)
                print '(''pt_coords(1,2): '',F6.3)', pt_coords(1,2)
                print '(''pt_coords(2,1): '',F6.3)', pt_coords(2,1)
                print '(''pt_coords(2,2): '',F6.3)', pt_coords(2,2)
                print '(''size_x: '',I6)', size_x
                print '(''size_y: '',I6)', size_y
                stop ''
                
             end if

          end if

        end subroutine check_layer_bf


        function no_nogrdpt(grdpts_id,i,j)

          implicit none

          integer       , dimension(:,:), intent(in) :: grdpts_id
          integer(ikind)                , intent(in) :: i
          integer(ikind)                , intent(in) :: j
          logical                                    :: no_nogrdpt

          integer(ikind) :: i1,j1

          no_nogrdpt = .true.


          do j1=j-1,j+1
             do i1=i-1,i+1
                
                if(grdpts_id(i1,j1).eq.no_pt) then
                   no_nogrdpt = .false.
                   exit
                end if

             end do

             if(.not.no_nogrdpt) then
                exit
             end if

          end do          

        end function no_nogrdpt

      end module bf_remove_module
