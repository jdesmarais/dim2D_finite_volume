      !> @file
      !> module encapsulating the subroutines for the determination of
      !> the boundary layers computed in the interior domain when some
      !> buffer layers are activated
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines for the determination of
      !> the boundary layers computed in the interior domain when some
      !> buffer layers are activated
      !
      !> @date
      !> 29_10_2014 - initial version         - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interior_bc_sections_module
      
        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type

        use parameters_input, only :
     $       nx,
     $       ny,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public ::
     $       ini_interior_bc_sections,
     $       determine_interior_bc_sections,
     $       close_last_bc_section,
     $       set_full_interior_bc_section,
     $       minimize_interior_bc_section,
     $       process_bc_sections_into_bc_procedure


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the parameters when determining the interior
        !> boundary layers
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary layers computed by the interior domain
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !
        !>@param no_bf_common_with_interior
        !> logical indicating whether there are buffer layers that have
        !> grid points in common with the interior domain
        !--------------------------------------------------------------
        subroutine ini_interior_bc_sections(
     $       nb_bc_sections,
     $       min_initialized,
     $       max_initialized,
     $       no_bf_common_with_interior)

          implicit none

          integer, intent(out) :: nb_bc_sections
          logical, intent(out) :: min_initialized
          logical, intent(out) :: max_initialized
          logical, intent(out) :: no_bf_common_with_interior

          nb_bc_sections             = 0
          min_initialized            = .false.
          max_initialized            = .false.
          no_bf_common_with_interior = .true.

        end subroutine ini_interior_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the limits of the current boundary section by
        !> analyzing the alignment of the buffer layer investigated
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param bf_alignment
        !> alignment of the buffer layer in the direction of the
        !> bc_section investigated (ex: for North buffer layers,
        !> the borders of the boundary sections that are relevant are
        !> the borders in the x-direction as the borders in the
        !> y-direction are always known)
        !
        !>@param interior_inf
        !> lower border for the interior boundary sections in the
        !> direction of the bc_section investigated (ex: for N and S,
        !> the lower border is 1 but for E and W borders, the lower
        !> border is bc_size+1)
        !
        !>@param interior_sup
        !> upper border for the interior boundary sections in the
        !> direction of the bc_section investigated (ex: for N and S,
        !> the upper border is nx but for E and W borders, the upper
        !> border is ny-bc_size)
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !
        !>@param no_bf_common_with_interior
        !> logical indicating whether there are buffer layers that have
        !> grid points in common with the interior domain
        !--------------------------------------------------------------
        subroutine determine_interior_bc_sections(
     $     bf_alignment,
     $     interior_inf,
     $     interior_sup,
     $     nb_bc_sections,
     $     bc_sections,
     $     min_initialized,
     $     max_initialized,
     $     no_bf_common_with_interior)

          implicit none

          integer(ikind), dimension(2)               , intent(in)    :: bf_alignment
          integer(ikind)                             , intent(in)    :: interior_inf
          integer(ikind)                             , intent(in)    :: interior_sup
          integer                                    , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(inout) :: min_initialized
          logical                                    , intent(inout) :: max_initialized
          logical                                    , intent(inout) :: no_bf_common_with_interior


          !does the buffer layer investigated has grid points
          !in common with the interior ?
          if((
     $         min(interior_sup,(bf_alignment(2)+bc_size))-
     $         max(interior_inf,(bf_alignment(1)-bc_size))).ge.0) then

             no_bf_common_with_interior = .false.

             
             !if the boundary layer covers the inferior part of
             !the interior domain, the superior limit of the buffer layer
             !can be used as the lower limit of the current bc_section
             if((bf_alignment(1)-bc_size).le.interior_inf) then

                ! |/////////////|
                !   |                |
                if((bf_alignment(2)+bc_size).lt.interior_sup) then
                   call set_as_min(
     $                  nb_bc_sections,
     $                  bc_sections,
     $                  bf_alignment(2)+bc_size+1,
     $                  min_initialized,
     $                  max_initialized)

                ! |////////////////////|
                !   |                |
                end if

             else
                
                !        |/////////////|
                !   |                |
                if(.not.min_initialized) then
                   call set_as_min(
     $                  nb_bc_sections,
     $                  bc_sections,
     $                  interior_inf,
     $                  min_initialized,
     $                  max_initialized)

                end if

                call set_as_max(
     $               nb_bc_sections,
     $               bc_sections,
     $               bf_alignment(1)-bc_size-1,
     $               min_initialized,
     $               max_initialized)

                !        |/////////|
                !   |                |
                if((bf_alignment(2)+bc_size).lt.interior_sup) then
                   call set_as_min(
     $                  nb_bc_sections,
     $                  bc_sections,
     $                  bf_alignment(2)+bc_size+1,
     $                  min_initialized,
     $                  max_initialized)

                !        |/////////////|
                !   |                |   
                end if

             end if

          end if

        end subroutine determine_interior_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> finalize the determination of the interior boundary sections
        !> by deciding the maximum border of the last unfinished
        !> boundary section if any
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param interior_sup
        !> upper border for the interior boundary sections in the
        !> direction of the bc_section investigated (ex: for N and S,
        !> the upper border is nx but for E and W borders, the upper
        !> border is ny-bc_size)
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !--------------------------------------------------------------
         subroutine close_last_bc_section(
     $     nb_bc_sections,
     $     bc_sections,
     $     interior_sup,
     $     min_initialized,
     $     max_initialized)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections
          integer(ikind)                , intent(in)    :: interior_sup
          logical                       , intent(inout) :: min_initialized
          logical                       , intent(inout) :: max_initialized

          if(min_initialized.and.(.not.max_initialized)) then
             
             call set_as_max(
     $         nb_bc_sections,bc_sections,interior_sup,
     $         min_initialized, max_initialized)
                
          end if

        end subroutine close_last_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the minimum border of the current interior boundary
        !> section
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param min
        !> value set as minimum of the current interior boundary section
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !--------------------------------------------------------------
        subroutine set_as_min(
     $     nb_bc_sections, bc_sections, min,
     $     min_initialized, max_initialized)

          implicit none

          integer                                    , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind)                             , intent(in)    :: min
          logical                                    , intent(inout) :: min_initialized
          logical                                    , intent(inout) :: max_initialized

          integer, parameter :: nb_bc_sections_alloc=5
          integer(ikind), dimension(:,:), allocatable :: bc_sections_tmp

          if(allocated(bc_sections)) then
             if(nb_bc_sections+1>size(bc_sections,2)) then

                allocate(bc_sections_tmp(2,nb_bc_sections+nb_bc_sections_alloc))
                bc_sections_tmp(:,1:size(bc_sections,2)) = bc_sections(:,:)
                call MOVE_ALLOC(bc_sections_tmp,bc_sections)

             end if

          else
             allocate(bc_sections(2,nb_bc_sections_alloc))
          end if

          bc_sections(1,nb_bc_sections+1) = min

          min_initialized = .true.

          call check_bc_section(nb_bc_sections,min_initialized,max_initialized)

        end subroutine set_as_min


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the maximum border of the current interior boundary
        !> section
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param max
        !> value set as maximum of the current interior boundary section
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !--------------------------------------------------------------
        subroutine set_as_max(
     $     nb_bc_sections,bc_sections,max,
     $     min_initialized, max_initialized)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections
          integer(ikind)                , intent(in)    :: max
          logical                       , intent(inout) :: min_initialized
          logical                       , intent(inout) :: max_initialized

          bc_sections(2,nb_bc_sections+1) = max

          max_initialized = .true.

          call check_bc_section(nb_bc_sections,min_initialized,max_initialized)

        end subroutine set_as_max 


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the current interior boundary layer is
        !> completed and if so increment the number of boundary
        !> sections, and reinitialize the min_initialized and
        !> max_initialized booleans
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !--------------------------------------------------------------
        subroutine check_bc_section(
     $     nb_bc_sections,
     $     min_initialized,
     $     max_initialized)

          implicit none

          integer, intent(inout) :: nb_bc_sections
          logical, intent(inout) :: min_initialized
          logical, intent(inout) :: max_initialized

          if(min_initialized.and.max_initialized) then
             nb_bc_sections = nb_bc_sections+1
             min_initialized = .false.
             max_initialized = .false.
          end if

        end subroutine check_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> if there are no buffer layers that have grid points in common
        !> with the interior domain, the boundary layer computed by the
        !> interior domain should be set as the full boundary section
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param min_initialized
        !> logical indicating whether the lower border of the current
        !> boundary layer is initialized
        !
        !>@param max_initialized
        !> logical indicating whether the upper border of the current
        !> boundary layer is initialized
        !
        !>@param interior_inf
        !> lower border for the interior boundary sections in the
        !> direction of the bc_section investigated (ex: for N and S,
        !> the lower border is 1 but for E and W borders, the lower
        !> border is bc_size+1)
        !
        !>@param interior_sup
        !> upper border for the interior boundary sections in the
        !> direction of the bc_section investigated (ex: for N and S,
        !> the upper border is nx but for E and W borders, the upper
        !> border is ny-bc_size)
        !--------------------------------------------------------------
        subroutine set_full_interior_bc_section(
     $     nb_bc_sections,
     $     bc_sections,
     $     interior_inf,
     $     interior_sup)

          implicit none

          integer                                    , intent(out)   :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind)                             , intent(in)    :: interior_inf
          integer(ikind)                             , intent(in)    :: interior_sup

          if(allocated(bc_sections)) then
             deallocate(bc_sections)
          end if

          allocate(bc_sections(2,1))
          bc_sections(1,1) = interior_inf
          bc_sections(2,1) = interior_sup

          nb_bc_sections = 1

        end subroutine set_full_interior_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> minimze the interior boundary sections by optimizing the
        !> storage: the number of boundary sections saved in the
        !> array may be smaller than the extent of the array
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !--------------------------------------------------------------
        subroutine minimize_interior_bc_section(
     $     nb_bc_sections,
     $     bc_sections)

          implicit none

          integer                                    , intent(in) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections

          integer(ikind), dimension(:,:), allocatable :: bc_sections_temp

          if(nb_bc_sections.ne.size(bc_sections_temp,2)) then

             allocate(bc_sections_temp(2,nb_bc_sections))
             bc_sections_temp(:,:) = bc_sections(:,1:nb_bc_sections)
             call MOVE_ALLOC(bc_sections_temp,bc_sections)
             
          end if

        end subroutine minimize_interior_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> convert the information on the extent of each boundary
        !> sections into information for the procedure to be applied
        !> on each section
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_bc_sections_N
        !> integer table containing the borders of the north boundary
        !> sections computed by the interior domain (ex:
        !> bc_sections(:,1): min and max borders of the first boundary
        !> section)
        !
        !>@param interior_bc_sections_S
        !> integer table containing the borders of the south boundary
        !> sections computed by the interior domain
        !
        !>@param interior_bc_sections_E
        !> integer table containing the borders of the east boundary
        !> sections computed by the interior domain
        !
        !>@param interior_bc_sections_W
        !> integer table containing the borders of the west boundary
        !> sections computed by the interior domain
        !
        !>@param interior_bc_sections
        !> integer table with the procedures to be applied for each
        !> interior boundary section
        !--------------------------------------------------------------
        subroutine process_bc_sections_into_bc_procedure(
     $     interior_bc_sections_N,
     $     interior_bc_sections_S,
     $     interior_bc_sections_E,
     $     interior_bc_sections_W,
     $     interior_bc_sections)

          implicit none

          integer(ikind), dimension(:,:), allocatable, intent(in)  :: interior_bc_sections_N
          integer(ikind), dimension(:,:), allocatable, intent(in)  :: interior_bc_sections_S
          integer(ikind), dimension(:,:), allocatable, intent(in)  :: interior_bc_sections_E
          integer(ikind), dimension(:,:), allocatable, intent(in)  :: interior_bc_sections_W
          integer(ikind), dimension(:,:), allocatable, intent(out) :: interior_bc_sections

          integer :: corner_nb
          integer :: nb_bc_sections

          !evaluate the total number of bc_sections
          corner_nb = 0          
          call add_corner(corner_nb,interior_bc_sections_S)
          call add_corner(corner_nb,interior_bc_sections_N)

          nb_bc_sections = corner_nb
          call add_edge(nb_bc_sections,interior_bc_sections_N)
          call add_edge(nb_bc_sections,interior_bc_sections_S)
          call add_edge(nb_bc_sections,interior_bc_sections_E)
          call add_edge(nb_bc_sections,interior_bc_sections_W)

          
          if(nb_bc_sections.ne.0) then

             !allocate spaces for the bc_sections
             if(allocated(interior_bc_sections)) then
                if(size(interior_bc_sections,2).ne.nb_bc_sections) then
                   deallocate(interior_bc_sections)
                   allocate(interior_bc_sections(4,nb_bc_sections))
                end if
             else
                allocate(interior_bc_sections(4,nb_bc_sections))
             end if

             !fill the bc_sections with the procedures
             nb_bc_sections = 0

             !process the S_edge
             if(allocated(interior_bc_sections_S)) then
                
                call process_bc_sections_NS(
     $               nb_bc_sections,
     $               interior_bc_sections,
     $               interior_bc_sections_S,
     $               1,
     $               SW_corner_type,
     $               S_edge_type,
     $               SE_corner_type)

             end if

             !process the W_edge
             if(allocated(interior_bc_sections_W)) then

                call process_bc_sections_EW(
     $               nb_bc_sections,
     $               interior_bc_sections,
     $               interior_bc_sections_W,
     $               1,
     $               W_edge_type)

             end if

             !process the E_edge
             if(allocated(interior_bc_sections_E)) then

                call process_bc_sections_EW(
     $               nb_bc_sections,
     $               interior_bc_sections,
     $               interior_bc_sections_E,
     $               nx-1,
     $               E_edge_type)

             end if

             !process the N_edge
             if(allocated(interior_bc_sections_N)) then
                
                call process_bc_sections_NS(
     $               nb_bc_sections,
     $               interior_bc_sections,
     $               interior_bc_sections_N,
     $               ny-bc_size+1,
     $               NW_corner_type,
     $               N_edge_type,
     $               NE_corner_type)

             end if

          else

             if(allocated(interior_bc_sections)) then
                deallocate(interior_bc_sections)
             end if

          end if
          

        end subroutine process_bc_sections_into_bc_procedure


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether there are corner procedures in the
        !> boundary sections of the layer and whether there will
        !> be separated from the edge computations and so be counted
        !> as an extra element for the procedures
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_bc_sections
        !> integer table containing the borders of the south boundary
        !> sections computed by the interior domain
        !
        !>@param corner_nb
        !> number of corners requiring an extra procedure element
        !--------------------------------------------------------------
        subroutine add_corner(corner_nb,interior_bc_sections)
        
          implicit none

          integer                                    , intent(inout) :: corner_nb
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: interior_bc_sections

          
          !shall the SW corner procedure be counted as
          !an additional bc_section element ?
          if(allocated(interior_bc_sections)) then
             if(
     $            (interior_bc_sections(1,1).eq.1).and.
     $            (interior_bc_sections(2,1).ge.(bc_size+1))) then

                corner_nb = corner_nb+1

             end if


          !shall the SE corner procedure be counted as
          !an additional bc_section element ?
             if(
     $            (interior_bc_sections(1,size(interior_bc_sections,2)).le.(nx-bc_size)).and.
     $            (interior_bc_sections(2,size(interior_bc_sections,2)).eq.nx)) then

                corner_nb = corner_nb+1

             end if
          end if

        end subroutine add_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the number of procedure elements by counting the number
        !> of edge procedure elements in the interior boundary sections
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param interior_bc_sections
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain
        !
        !>@param nb_bc_sections
        !> number of procedure elements
        !--------------------------------------------------------------
        subroutine add_edge(nb_bc_sections,interior_bc_sections)

          implicit none

          integer                                    , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: interior_bc_sections

          if(allocated(interior_bc_sections)) then
             nb_bc_sections = nb_bc_sections + size(interior_bc_sections,2)
          end if

        end subroutine add_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> translate the borders of the (North or South) boundary
        !> sections into procedure elements
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param interior_bc_sections
        !> integer table containing the borders of the procedure
        !> elements for the interior domain
        !
        !>@param interior_bc_sections_NS
        !> integer table containing the borders of the boundary
        !> sections
        !
        !>@param j_min
        !> integer identifying the lower border of the boundary
        !> sections along the y-direction (1 for south, ny-bc_size+1
        !> for north)
        !
        !>@param left_corner_type
        !> procedure type for the left corner (SW_corner_type for
        !> south, NW_corner_type for north)
        !
        !>@param edge_type
        !> procedure type for the edge (S_edge_type for south,
        !> N_edge_type for north)
        !
        !>@param right_corner_type
        !> procedure type for the right corner (SE_corner_type for
        !> south, NE_corner_type for north)
        !--------------------------------------------------------------
        subroutine process_bc_sections_NS(
     $     nb_bc_sections,
     $     interior_bc_sections,
     $     interior_bc_sections_NS,
     $     j_min,
     $     left_corner_type,
     $     edge_type,
     $     right_corner_type)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: interior_bc_sections
          integer(ikind), dimension(:,:), intent(in)    :: interior_bc_sections_NS
          integer(ikind)                , intent(in)    :: j_min
          integer                       , intent(in)    :: left_corner_type
          integer                       , intent(in)    :: edge_type
          integer                       , intent(in)    :: right_corner_type

          integer :: nb_sec
          integer :: k_match
          integer :: k


          nb_sec = size(interior_bc_sections_NS,2)

          !special treatment for the left corner
          if(interior_bc_sections_NS(1,1).eq.1) then

             nb_bc_sections=nb_bc_sections+1

             interior_bc_sections(:,nb_bc_sections) = 
     $            [left_corner_type,1,j_min,0]
             
             if(interior_bc_sections_NS(2,1).ge.(bc_size+1)) then

                call add_first_edge(
     $               nb_bc_sections,
     $               interior_bc_sections,
     $               interior_bc_sections_NS,
     $               bc_size+1,
     $               j_min,
     $               edge_type,
     $               right_corner_type)

             end if

          else

             call add_first_edge(
     $            nb_bc_sections,
     $            interior_bc_sections,
     $            interior_bc_sections_NS,
     $            interior_bc_sections_NS(1,1),
     $            j_min,
     $            edge_type,
     $            right_corner_type)
             
          end if

          !interior edge
          k_match = nb_bc_sections-1

          do k=2, size(interior_bc_sections_NS,2)-1

             interior_bc_sections(:,k_match+k) = [
     $            edge_type,
     $            interior_bc_sections_NS(1,k),
     $            j_min,
     $            interior_bc_sections_NS(2,k)]
          end do

          if(nb_sec.ge.3) then
             nb_bc_sections =
     $            nb_bc_sections + 
     $            nb_sec-2
          end if
          
          !special treament for the right corner
          if(nb_sec.ne.1) then

             if(interior_bc_sections_NS(1,nb_sec).le.(nx-bc_size)) then

                if(interior_bc_sections_NS(2,nb_sec).eq.nx) then
                   nb_bc_sections=nb_bc_sections+1
                
                   interior_bc_sections(:,nb_bc_sections) = [
     $                  edge_type,
     $                  interior_bc_sections_NS(1,nb_sec),
     $                  j_min,
     $                  nx-bc_size]

                   nb_bc_sections=nb_bc_sections+1
                
                   interior_bc_sections(:,nb_bc_sections) = [
     $                  right_corner_type,
     $                  nx-bc_size+1,
     $                  j_min,
     $                  0]

                else
                   nb_bc_sections=nb_bc_sections+1
                
                   interior_bc_sections(:,nb_bc_sections) = [
     $                  edge_type,
     $                  interior_bc_sections_NS(1,nb_sec),
     $                  j_min,
     $                  interior_bc_sections_NS(2,nb_sec)]

                end if
             
             else
                if(interior_bc_sections_NS(2,nb_sec).eq.nx) then

                   nb_bc_sections=nb_bc_sections+1
                
                   interior_bc_sections(:,nb_bc_sections) = [
     $                  right_corner_type,
     $                  nx-bc_size+1,
     $                  j_min,
     $                  0]

                end if
             end if

          end if

        end subroutine process_bc_sections_NS


                !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> translate the borders of the (East or West) boundary
        !> sections into procedure elements
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of boundary sections computed by the interior_domain
        !
        !>@param interior_bc_sections
        !> integer table containing the borders of the procedure
        !> elements for the interior domain
        !
        !>@param interior_bc_sections_EW
        !> integer table containing the borders of the boundary
        !> sections
        !
        !>@param i_min
        !> integer identifying the lower border of the boundary
        !> sections along the x-direction (1 for west, nx-bc_size+1
        !> for east)
        !
        !>@param edge_type
        !> procedure type for the edge (W_edge_type for west,
        !> E_edge_type for east)
        !--------------------------------------------------------------
        subroutine process_bc_sections_EW(
     $     nb_bc_sections,
     $     interior_bc_sections,
     $     interior_bc_sections_EW,
     $     i_min,
     $     edge_type)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: interior_bc_sections
          integer(ikind), dimension(:,:), intent(in)    :: interior_bc_sections_EW
          integer(ikind)                , intent(in)    :: i_min
          integer                       , intent(in)    :: edge_type

          integer :: k_match
          integer :: k

          k_match = nb_bc_sections

          do k=1, size(interior_bc_sections_EW,2)

             interior_bc_sections(:,k_match+k) = [
     $            edge_type,
     $            i_min,
     $            interior_bc_sections_EW(1,k),
     $            interior_bc_sections_EW(2,k)]

          end do

          nb_bc_sections = nb_bc_sections +
     $         size(interior_bc_sections_EW,2)

        end subroutine process_bc_sections_EW


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> convert the first boundary element into procedure elements:
        !> depending whether the boundary element contains a right
        !> corner, the boundary element is split into two procedure
        !> elements: an edge and a corner procedure elements
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !>@param nb_bc_sections
        !> number of procedure elements
        !
        !>@param interior_bc_sections
        !> integer table containing the procedure elements
        !
        !>@param interior_bc_sections_NS
        !> integer table containing the borders of the boundary
        !> sections computed by the interior domain (either north
        !> or south)
        !
        !>@param i_min
        !> integer identifying the lower border of the boundary layer
        !> in the x-direction
        !
        !>@param j_min
        !> integer identifying the lower border of the boundary layer
        !> in the y-direction (1 for south, ny-bc_size+1 for north)
        !
        !>@param edge_type
        !> integer identifying the procedure applied for the edge
        !> (S_edge for south, N_edge for north)
        !
        !>@param corner_type
        !> integer identifying the procedure applied for the right
        !> corner (SE_corner for south, NE_corner for north)
        !--------------------------------------------------------------
        subroutine add_first_edge(
     $     nb_bc_sections,
     $     interior_bc_sections,
     $     interior_bc_sections_NS,
     $     i_min,
     $     j_min,
     $     edge_type,
     $     corner_type)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: interior_bc_sections
          integer(ikind), dimension(:,:), intent(in)    :: interior_bc_sections_NS
          integer(ikind)                , intent(in)    :: i_min
          integer(ikind)                , intent(in)    :: j_min
          integer                       , intent(in)    :: edge_type
          integer                       , intent(in)    :: corner_type

          !check if the right corner is also part of
          !the first bc_section element
          if(interior_bc_sections_NS(2,1).eq.nx) then

             if(i_min.le.(nx-bc_size)) then
                nb_bc_sections = nb_bc_sections+1
                interior_bc_sections(:,nb_bc_sections) = [
     $               edge_type,
     $               i_min,
     $               j_min,
     $               nx-bc_size]
             end if
             
             nb_bc_sections = nb_bc_sections+1
             interior_bc_sections(:,nb_bc_sections) = [
     $            corner_type,
     $            nx-bc_size+1,
     $            j_min,
     $            0]
             
          else
             nb_bc_sections = nb_bc_sections+1
             interior_bc_sections(:,nb_bc_sections) = [
     $            edge_type,
     $            i_min,
     $            j_min,
     $            interior_bc_sections_NS(2,1)]
             
          end if

        end subroutine add_first_edge

      end module bf_interior_bc_sections_module
