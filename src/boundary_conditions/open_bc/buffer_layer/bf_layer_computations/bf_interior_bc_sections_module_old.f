      module bf_interior_bc_sections_module

        use parameters_input, only :
     $     bc_size

        use parameters_kind, only :
     $     ikind,
     $     rkind

        implicit none

        private
        public ::
     $       ini_interior_bc_sections,
     $       determine_interior_bc_sections,
     $       close_last_bc_section,
     $       set_full_interior_bc_section


        contains

        subroutine ini_interior_bc_sections(
     $       nb_bc_sections,
     $       min_initialized,
     $       no_bf_common_with_interior)

          implicit none

          integer, intent(out) :: nb_bc_sections
          logical, intent(out) :: min_initialized
          logical, intent(out) :: no_bf_common_with_interior

          nb_bc_sections  = 0
          min_initialized = .false.
          no_bf_common_with_interior = .true.

        end subroutine ini_interior_bc_sections


        !bf_alignment:   alignment in the direction of
        !                the bc_section investigated
        !interior_inf:   minimum in the direction of
        !                the bc_section investigated
        !interior_sup:   maximum in the direction of
        !                the bc_section investigated
        !nb_bc_sections: number of boundary sections
        !                computed by the interior_domain
        !bc_sections:    boundary sections computed by
        !                the interior domain
        subroutine determine_interior_bc_sections(
     $     bf_alignment,
     $     interior_inf,
     $     interior_sup,
     $     nb_bc_sections,
     $     bc_sections,
     $     min_initialized,
     $     no_bf_common_with_interior,
     $     end_of_bc_sections)

          implicit none

          integer(ikind), dimension(2)               , intent(in)    :: bf_alignment
          integer(ikind)                             , intent(in)    :: interior_inf
          integer(ikind)                             , intent(in)    :: interior_sup
          integer                                    , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(inout) :: min_initialized
          logical                                    , intent(inout) :: no_bf_common_with_interior
          logical                                    , intent(out)   :: end_of_bc_sections


          !does the buffer layer investigated has grid points
          !in common with the interior ?
          if((
     $         min(interior_sup,(bf_alignment(2)+bc_size))-
     $         max(interior_inf,(bf_alignment(1)-bc_size))).gt.0) then

             no_bf_common_with_interior = .false.
             end_of_bc_sections = .true.

             !is it the first boundary section of the interior ?
             if(nb_bc_sections.eq.0) then
                
                !if the boundary layer does not cover the inferior part of
                !the interior domain, the inferior limit of the buffer layer
                !can be used as the upper limit of the current bc_section
                if(.not.((bf_alignment(1)-bc_size).le.interior_inf)) then
                   
                   call set_as_min(nb_bc_sections,bc_sections,interior_inf)
                   call set_as_max(nb_bc_sections,bc_sections,bf_alignment(1)-bc_size-1)
                   
                   !the current interior_bc_section is completed and the
                   !number of bc_sections can be incremented
                   nb_bc_sections = nb_bc_sections+1

                end if

                !check the upper bound of the boundary layer
                call check_upper_bound_bc_section(
     $               bf_alignment,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               end_of_bc_sections)
                

             !the boundary section investigated is not the first one
             else

                !the inferior limit of the buffer layer can be used as the
                !upper bound of the current bc_section
                call set_as_max(nb_bc_sections,bc_sections,bf_alignment(1)-bc_size-1)

                !the current interior_bc_section is completed and the
                !number of bc_sections can be incremented
                nb_bc_sections = nb_bc_sections+1

                !check the upper bound of the boundary layer
                call check_upper_bound_bc_section(
     $               bf_alignment,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               end_of_bc_sections)

             end if

          else

             end_of_bc_sections = .false.

          end if

        end subroutine determine_interior_bc_sections


        subroutine check_upper_bound_bc_section(
     $     bf_alignment,
     $     interior_sup,
     $     nb_bc_sections,
     $     bc_sections,
     $     min_initialized,
     $     end_of_bc_sections)

          implicit none
          
          integer(ikind), dimension(2)               , intent(in)    :: bf_alignment
          integer(ikind)                             , intent(in)    :: interior_sup
          integer                                    , intent(in)    :: nb_bc_sections
          integer       , dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(inout) :: min_initialized
          logical                                    , intent(out)   :: end_of_bc_sections


          !is the upper limit of the buffer layer in
          !the interior domain ?
          if((bf_alignment(2)+bc_size).lt.interior_sup) then
             
             !the maximum of the buffer layer is used as the
             !lower limit of the current bc_section
             call set_as_min(nb_bc_sections,bc_sections,bf_alignment(2)+bc_size+1)

             min_initialized    = .true.
             end_of_bc_sections = .false.

          !otherwise the buffer layer covers entirely
          !the interior boundary layer: there are
          !no bc_sections to be computed by the interior
          else

             end_of_bc_sections = .true.

          end if

        end subroutine check_upper_bound_bc_section


        subroutine close_last_bc_section(
     $     nb_bc_sections,
     $     bc_sections,
     $     interior_sup,
     $     min_initialized,
     $     end_of_bc_sections)

          implicit none

          integer                       , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections
          integer(ikind)                , intent(in)    :: interior_sup
          logical                       , intent(in)    :: min_initialized
          logical                       , intent(in)    :: end_of_bc_sections

          if(min_initialized.and.(.not.end_of_bc_sections)) then
             
             call set_as_max(nb_bc_sections,bc_sections,interior_sup)
             nb_bc_sections = nb_bc_sections+1
                
          end if

        end subroutine close_last_bc_section


        subroutine set_as_min(nb_bc_sections,bc_sections,min)

          implicit none

          integer                                    , intent(in)    :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind)                             , intent(in)    :: min

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

        end subroutine set_as_min

      
        subroutine set_as_max(nb_bc_sections,bc_sections,max)

          implicit none

          integer                       , intent(in)    :: nb_bc_sections
          integer(ikind), dimension(:,:), intent(inout) :: bc_sections
          integer(ikind)                , intent(in)    :: max

          bc_sections(2,nb_bc_sections+1) = max

        end subroutine set_as_max 


        subroutine set_full_interior_bc_section(
     $     nb_bc_sections,
     $     bc_sections,
     $     interior_inf,
     $     interior_sup)

          implicit none

          integer                                    , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind)                             , intent(in)    :: interior_inf
          integer(ikind)                             , intent(in)    :: interior_sup

          call set_as_min(nb_bc_sections,bc_sections,interior_inf)
          call set_as_max(nb_bc_sections,bc_sections,interior_sup)
          nb_bc_sections = nb_bc_sections+1

        end subroutine set_full_interior_bc_section

      end module bf_interior_bc_sections_module
