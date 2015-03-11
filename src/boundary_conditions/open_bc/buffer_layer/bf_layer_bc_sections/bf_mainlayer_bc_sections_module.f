      !> @file
      !> module encapsulating the subroutines for the determination of
      !> the boundary sections computed in the interior domain when 
      !> main layers are activated
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines for the determination of
      !> the boundary sections computed in the interior domain when 
      !> main layers are activated
      !
      !> @date
      !> 11_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_bc_sections_module

        use bf_interior_bc_sections_module, only :
     $       process_interior_bc_sections_into_bc_procedures

        use bf_mainlayer_pointer_class, only :
     $       bf_mainlayer_pointer

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none

        private
        public ::
     $       determine_interior_bc_sections_from_mainlayers,
     $       update_interior_bc_sections_from_mainlayers
        

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the extent of the boundary sections for the
        !> interior domain from the 4 mainlayers (N,S,E,W)
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !> @param mainlayer_pointers
        !> array with the mainlayers
        !
        !> @param interior_bc_sections_N
        !> extent of the boundary section for the North boundary
        !
        !> @param interior_bc_sections_S
        !> extent of the boundary section for the South boundary
        !
        !> @param interior_bc_sections_E
        !> extent of the boundary section for the East boundary
        !
        !> @param interior_bc_sections_W
        !> extent of the boundary section for the West boundary
        !------------------------------------------------------------
        subroutine determine_interior_bc_sections_from_mainlayers(
     $     mainlayer_pointers,
     $     interior_bc_sections_N,
     $     interior_bc_sections_S,
     $     interior_bc_sections_E,
     $     interior_bc_sections_W)

         implicit none

         type(bf_mainlayer_pointer), dimension(4)   , intent(in)    :: mainlayer_pointers
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_N
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_S
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_E
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_W


         type(bf_mainlayer), pointer :: get_mainlayer


         !North boundary layer
         !---------------------
         !if there are buffer sub-layers saved in the buffer main layer
         !the extents of the interior buffer layers are computed using
         !the buffer main layer
         if(mainlayer_pointers(N)%associated_ptr()) then
            get_mainlayer => mainlayer_pointers(N)%get_ptr()

            call get_mainlayer%update_interior_bc_sections(
     $           interior_bc_sections_N)

         !otherwise, there are no buffer layers and the extent of the
         !north interior boundary layer computed using the boundary
         !conditions is the entire north boundary layer
         else

            allocate(interior_bc_sections_N(2,1))
            interior_bc_sections_N(:,1) = [1,nx]

         end if


         !South boundary layer
         !---------------------
         !if there are buffer sub-layers saved in the buffer main layer
         !the extents of the interior buffer layers are computed using
         !the buffer main layer
         if(mainlayer_pointers(S)%associated_ptr()) then
            get_mainlayer => mainlayer_pointers(S)%get_ptr()

            call get_mainlayer%update_interior_bc_sections(
     $           interior_bc_sections_S)

         !otherwise, there are no buffer layers and the extent of the
         !south interior boundary layer computed using the boundary
         !conditions is the entire south boundary layer
         else

            allocate(interior_bc_sections_S(2,1))
            interior_bc_sections_S(:,1) = [1,nx]

         end if


         !East boundary layer
         !--------------------
         !if there are buffer sub-layers saved in the buffer main layer
         !the extents of the interior buffer layers are computed using
         !the buffer main layer
         if(mainlayer_pointers(E)%associated_ptr()) then
            get_mainlayer => mainlayer_pointers(E)%get_ptr()

            call get_mainlayer%update_interior_bc_sections(
     $           interior_bc_sections_E)

         !otherwise, there are no buffer layers and the extent of the
         !east interior boundary layer computed using the boundary
         !conditions is the entire east boundary layer
         else

            allocate(interior_bc_sections_E(2,1))
            interior_bc_sections_E(:,1) = [1+bc_size,ny-bc_size]

         end if


         !Wast boundary layer
         !--------------------
         !if there are buffer sub-layers saved in the buffer main layer
         !the extents of the interior buffer layers are computed using
         !the buffer main layer
         if(mainlayer_pointers(W)%associated_ptr()) then
            get_mainlayer => mainlayer_pointers(W)%get_ptr()

            call get_mainlayer%update_interior_bc_sections(
     $           interior_bc_sections_W)

         !otherwise, there are no buffer layers and the extent of the
         !west interior boundary layer computed using the boundary
         !conditions is the entire west boundary layer
         else

            allocate(interior_bc_sections_W(2,1))
            interior_bc_sections_W(:,1) = [1+bc_size,ny-bc_size]

         end if

        end subroutine determine_interior_bc_sections_from_mainlayers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the interior bc_sections from the 4 mainlayers
        !> (N,S,E,W)
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !>@param mainlayers
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !>@param interior_bc_sections
        !> extent of the boundary layers computed by the interior
        !> nodes
        !--------------------------------------------------------------
        subroutine update_interior_bc_sections_from_mainlayers(
     $     mainlayer_pointers,
     $     interior_bc_sections)

          implicit none

          type(bf_mainlayer_pointer), dimension(4)   , intent(in)  :: mainlayer_pointers
          integer(ikind), dimension(:,:), allocatable, intent(out) :: interior_bc_sections

          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_N
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_S
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_E
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_W


          !> determine the extents of the boundary sections in
          !> each main layer
          call determine_interior_bc_sections_from_mainlayers(
     $         mainlayer_pointers,
     $         interior_bc_sections_N,
     $         interior_bc_sections_S,
     $         interior_bc_sections_E,
     $         interior_bc_sections_W)
          
          !> deduce the boundary procedures applied at the edges
          !> of the interior domain
          call process_interior_bc_sections_into_bc_procedures(
     $         interior_bc_sections_N,
     $         interior_bc_sections_S,
     $         interior_bc_sections_E,
     $         interior_bc_sections_W,
     $         interior_bc_sections)

          !> clean the allocated arrays
          if(allocated(interior_bc_sections_N)) deallocate(interior_bc_sections_N)
          if(allocated(interior_bc_sections_S)) deallocate(interior_bc_sections_S)
          if(allocated(interior_bc_sections_E)) deallocate(interior_bc_sections_E)
          if(allocated(interior_bc_sections_W)) deallocate(interior_bc_sections_W)
          
        end subroutine update_interior_bc_sections_from_mainlayers

      end module bf_mainlayer_bc_sections_module
