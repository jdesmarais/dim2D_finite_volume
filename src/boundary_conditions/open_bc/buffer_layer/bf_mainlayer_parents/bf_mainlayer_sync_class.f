      !> @file
      !> bf_mainlayer_print object augmented with synchronization
      !> procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_mainlayer_print object augmented with synchronization
      !> procedures
      !
      !> @date
      ! 11_04_2014 - initial version      - J.L. Desmarais
      ! 26_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_sync_class
      
        use bf_interior_bc_sections_module, only :
     $     ini_interior_bc_sections,
     $     determine_interior_bc_sections,
     $     close_last_bc_section,
     $     set_full_interior_bc_section,
     $     minimize_interior_bc_section

        use bf_mainlayer_print_class, only :
     $       bf_mainlayer_print

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: bf_mainlayer_sync
        
        
        !> @class bf_mainlayer_sync
        !> class storing bf_sublayers corresponding
        !> to the same cardinal point (N,S,E,W,NE,NW,SE,SW)
        !
        !> @param determine_interior_bc_layers
        !> determine the interior_bc_sections, i.e. the boundary
        !> grid points of the interior domain that are computed
        !> by the interior domain and not exchanged with the
        !> buffer layers
        !
        !> @param sync_nodes_with_interior
        !> synchronise the nodes between the interior domain and
        !> the buffer layers constituing the buffer main layer
        !---------------------------------------------------------------
        type, extends(bf_mainlayer_print) :: bf_mainlayer_sync

          contains

          procedure, pass :: determine_interior_bc_layers
          procedure, pass :: sync_nodes_with_interior

        end type bf_mainlayer_sync


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the extent of the boundary layers computed
        !> by the interior nodes
        !
        !> @date
        !> 28_10_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param interior_bc_sections
        !> extent of the boundary layers computed by the interior
        !> nodes
        !--------------------------------------------------------------
        subroutine determine_interior_bc_layers(this, bc_sections)

          implicit none

          class(bf_mainlayer_sync)                   , intent(in)    :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections

          integer        :: dir
          integer(ikind) :: interior_inf
          integer(ikind) :: interior_sup

          integer :: nb_bc_sections
          logical :: min_initialized
          logical :: max_initialized
          logical :: no_bf_common_with_interior

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: k

          integer(ikind), dimension(2,2) :: bf_alignments
          integer(ikind), dimension(2)   :: bf_alignment


          !initialize the variables for the determination of the
          !interior boundary layers depending on the cardinal
          !coordinate of the main layer
          select case(this%mainlayer_id)
            case(N,S)
               dir          = 1
               interior_inf = 1
               interior_sup = nx

            case(E,W)
               dir          = 2
               interior_inf = bc_size+1
               interior_sup = ny-bc_size
               
            case default
               call error_mainlayer_id(
     $              'bf_mainlayer_class.f',
     $              'determine_interior_bc_sections',
     $              this%mainlayer_id)
               
          end select


          !initialize the interior_bc_sections
          call ini_interior_bc_sections(
     $         nb_bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_interior)


          !initialize the pointer to the first sublayer
          !investigated
          current_sublayer => this%head_sublayer


          !go through the chained list and update the
          !extents of the interior boundary layers
          if(this%nb_sublayers.gt.0) then

             do k=1, this%nb_sublayers

                !determine the alignment of the sublayer
                bf_alignments = current_sublayer%get_alignment_tab()

                !determine the alignment relevant for the interior
                bf_alignment(1) = bf_alignments(dir,1)
                bf_alignment(2) = bf_alignments(dir,2)

                !update the extent of the interior boundary
                !sections
                call determine_interior_bc_sections(
     $               bf_alignment,
     $               interior_inf,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               max_initialized,
     $               no_bf_common_with_interior)

                !get the next sublayer in the mainlayer
                current_sublayer => current_sublayer%get_next()

             end do

             !finalize the interior_bc_sections
             call close_last_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_sup,
     $            min_initialized,
     $            max_initialized)

          end if

          if(no_bf_common_with_interior) then
             call set_full_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_inf,
     $            interior_sup)
          else

            !minimize the extent of the interior boundary
            !layers
             call minimize_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections)

          end if

        end subroutine determine_interior_bc_layers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> exchange grid points between the buffer main layer and
        !> interior domain
        !
        !> @date
        !> 29_10_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param interior_nodes
        !> grid points from the interior domain
        !--------------------------------------------------------------
        subroutine sync_nodes_with_interior(this, interior_nodes)

          implicit none

          class(bf_mainlayer_sync)        , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i


          current_sublayer => this%head_sublayer

          do i=1, this%nb_sublayers

             call current_sublayer%sync_nodes_with_interior(
     $            interior_nodes)

             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine sync_nodes_with_interior


      end module bf_mainlayer_sync_class
      
