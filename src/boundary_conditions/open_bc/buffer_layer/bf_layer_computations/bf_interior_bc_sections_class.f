      module bf_interior_bc_sections_class

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_kind, only :
     $       ikind

        implicit none

        integer, parameter :: allocation_nb=5


        type :: interior_bc_sections

          integer(ikind), dimension(:,:), allocatable, private :: bc_sections
          integer                                    , private :: nb_sections

          procedure, pass :: ini
          procedure, pass :: process_E_and_W_sections
          procedure, pass :: add_E_edge
          procedure, pass :: add_W_edge
          procedure, pass :: add_EW_edge

          procedure, pass :: add_bc_section

        end type interior_bc_sections

        contains

        subroutine ini(this,nb_bc_sections)

          implicit none

          class(bf_interior_bc_sections_class), intent(inout) :: this

          allocate(this%bc_sections(4,nb_bc_sections))

        end subroutine ini

        subroutine process_E_and_W_sections(
     $     this,
     $     E_sections,
     $     W_sections)

          implicit none

          class(bf_interior_bc_sections), intent(inout) :: this
          integer(ikind), dimension(:,:), intent(inout) :: E_sections
          integer(ikind), dimension(:,:), intent(inout) :: W_sections


          integer :: k1,k2

          k1=1
          k2=1

          do while((k1.le.size(E_sections,2)).and.
     $         (k2.le.size(W_sections,2)))

             !check whether E_sections(:,k1) and
             !W_sections(:,k2) have grid points in
             !common
             common_grdpts = are_grdpts_common(
     $            E_sections(:,k1),
     $            W_sections(:,k2))

             !if there are grdpts in common, the
             !bc_sections are merged
             if(common_grdpts) then
                call merge_E_and_W_sections(
     $               this,
     $               E_sections(:,k1),
     $               W_sections(:,k2))

             !if there are no grdpts in common,
             !one of the two sections should be
             !added to the 
             else

                if



          end do

        end subroutine process_E_and_W_sections


        subroutine add_E_edge(this,borders)

          implicit none

          class(bf_interior_bc_sections), intent(inout) :: this
          integer(ikind), dimension(2)  , intent(in)    :: borders

          integer(ikind), dimension(4) :: bc_section

          bc_section(1) = E_edge_type
          bc_section(2) = nx-1
          bc_section(3) = borders(1)
          bc_section(4) = borders(2)

          call this%add_bc_section(bc_section)

        end subroutine add_E_edge


        subroutine add_W_edge(this,borders)

          implicit none

          class(bf_interior_bc_sections), intent(inout) :: this
          integer(ikind), dimension(2)  , intent(in)    :: borders

          integer(ikind), dimension(4) :: bc_section

          bc_section(1) = W_edge_type
          bc_section(2) = 1
          bc_section(3) = borders(1)
          bc_section(4) = borders(2)

          call this%add_bc_section(bc_section)

        end subroutine add_W_edge


        subroutine add_bc_section(this,bc_section)

          implicit none

          class(bf_interior_bc_sections), intent(inout) :: this
          integer(ikind), dimension(4)  , intent(in)    :: bc_section

          this%nb_sections=this%nb_sections+1

          if(this%nb_sections>size(this%bc_sections,2)) then
             allocate(bc_sections_temp(4,size(this%bc_sections,2)+allocation_nb))
             bc_sections_temp(:,1:size(this%bc_sections,2)) = this%bc_sections(:,:)
             call MOVE_ALLOC(bc_sections_temp,this%bc_sections)
          end if

          this%bc_sections(:,this%nb_sections) = bc_section

        end subroutine add_bc_section


        function are_grdpts_common(bc_section1, bc_section2)
     $     result(grdpts_common)
        
          implicit none

          integer(ikind), dimension(2), intent(in) :: bc_section1
          integer(ikind), dimension(2), intent(in) :: bc_section2
          logical                                  :: grdpts_common
          
          grdpts_common = (
     $         min(bc_section1(2),bc_section2(2))-
     $         max(bc_section1(1),bc_section2(1))).gt.0

        end function are_grdpts_common

      end module bf_interior_bc_sections_class
