      !> @file
      !> module implementing the object encapsulating data
      !> required to decide which points need to be updated
      !> and whether this leads to the creation or the update
      !> of a buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating data
      !> required to decide which points need to be updated
      !> and whether this leads to the creation or the update
      !> of a buffer layer
      !
      !> @date
      ! 16_04_2014 - initial version      - J.L. Desmarais
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_path_icr_class

        use bf_layer_errors_module, only : error_mainlayer_id
        use bf_interface_class    , only : bf_interface
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_constant   , only : N,S,E,W,
     $                                     x_direction, y_direction,
     $                                     min_border, max_border
        use parameters_input      , only : nx,ny,ne,bc_size
        use parameters_kind       , only : ikind, rkind

        implicit none

        private
        public :: bf_path_icr


        !> @class bf_path_icr
        !> class encapsulating data required to decide which points
        !> need to be updated and whether this leads to the creation
        !> or the update of a buffer layer
        !
        !> @param matching_layer
        !> pointer to the buffer layer that should be updated
        !> if the pointer is not associated, it means that the
        !> buffer layer corresponding to this path should be
        !> created
        !
        !> @param mainlayer_id
        !> cardinal coordinate identifying in which direction the buffer
        !> layers may be updated by the path
        !
        !> @param pts
        !> allocatable table of integers identifying the bc_interior_pts
        !> analyzed that will be used in allocating or updating the 
        !> buffer layer
        !
        !> @param nb_pts
        !> integer identifying the bc_interior_pt contained in this
        !> path, it identifies which bc_interior_pt should be modified
        !> in the interior domain or inside the buffer layer
        !
        !> @param alignment
        !> table identifying the new position of the buffer layer and
        !> whether an existing buffer layer should be reallocated
        !
        !> @param ends
        !> logical identifying whether the path should ends or not
        !> and so whether the path should be processed to allocate
        !> and update an existing buffer layer
        !
        !>@param ini
        !> initialize the path
        !
        !>@param reinitialize
        !> re-initialize the path
        !
        !>@param is_ended
        !> get the ends attribute
        !
        !>@param add_pt_to_path
        !> add the general coordinates of a point
        !> to the current path
        !
        !>@param minimize_path
        !> reallocate the pts attribute to the nb_pts size
        !
        !>@param are_pts_far_away
        !> analyze whether the bc_interior_pt point is far away
        !> from the last point added to the path
        !
        !>@param are_pts_in_same_mainlayer
        !> check whether the bc_interior_pt analyzed is located in
        !> in the same buffer main layer as the grid points stored
        !> in the path
        !
        !>@param update_alignment
        !> update the alignment of the path. This data is
        !> used to decide whether the corresponding buffer layer
        !> will be allocated or reallocated and how its position
        !> may vary compared to the interior domain
        !
        !>@param analyze_path_first_pt
        !> analyze a bc_interior_pt that could be
        !> consider as the first point in the path
        !
        !>@param analyze_path_next_pt
        !> analyze a bc_interior_pt that should be
        !> consider as the next point in the path. The path
        !> is not empty
        !
        !>@param analyze_pt
        !> analyze the bc_interior_pt passed as argument
        !
        !>@param should_bf_sublayers_be_merged
        !> check whether the two buffer layers in clockwise direction
        !> should be merged
        !
        !>@param update_allocation_bf_sublayer
        !> update the allocation of the buffer layers to add
        !> new grid points
        !
        !>@param process_path
        !> process the data encapsulated in the path to update
        !> the buffer layers
        !
        !>@param get_nb_pts
        !> get the nb_pts attribute
        !
        !>@param get_pt
        !> get the general coordinates of the ith grid point
        !> stored in the path
        !
        !>@param print_on_screen
        !> print the content of the path on screen
        !---------------------------------------------------------------
        type :: bf_path_icr

          type(bf_sublayer), pointer :: matching_sublayer

          integer                                    , private :: mainlayer_id
          integer(ikind), dimension(:,:), allocatable, private :: pts
          integer                                    , private :: nb_pts
          integer(ikind), dimension(2,2)             , private :: alignment
          logical                                    , private :: ends

          contains
          
          procedure,   pass          :: ini
          procedure,   pass          :: reinitialize
          procedure,   pass          :: is_ended

          procedure,   pass, private :: add_pt_to_path
          procedure,   pass, private :: minimize_path
                       
          procedure,   pass, private :: are_pts_far_away
          procedure,   pass, private :: are_pts_in_same_mainlayer
          procedure,   pass, private :: update_alignment
                       
          procedure,   pass, private :: analyze_path_first_pt
          procedure,   pass, private :: analyze_path_next_pt
          procedure,   pass          :: analyze_pt

          procedure,   pass, private :: should_bf_sublayers_be_merged
          procedure,   pass, private :: update_allocation_bf_sublayer
          procedure,   pass          :: process_path

          procedure, pass :: get_nb_pts
          procedure, pass :: get_pt
          procedure, pass :: print_on_screen

        end type bf_path_icr


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(bf_path_icr), intent(inout) :: this

          nullify(this%matching_sublayer)
          
          this%ends      = .false.
          this%nb_pts    = 0

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> re-initialize the path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine reinitialize(this)
        
          implicit none

          class(bf_path_icr), intent(inout) :: this

          nullify(this%matching_sublayer)
          
          this%ends      = .false.
          this%nb_pts    = 0

        end subroutine reinitialize

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the ends attribute
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@return is_ended
        !> ends attribute
        !--------------------------------------------------------------
        function is_ended(this)

          implicit none

          class(bf_path_icr), intent(in) :: this
          logical                          :: is_ended

          is_ended = this%ends

        end function is_ended


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the general coordinates of a point
        !> to the current path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param point
        !> general coordinates of the grid point added to the path
        !--------------------------------------------------------------
        subroutine add_pt_to_path(this, point)

          implicit none

          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: point
          
          integer(ikind), dimension(:,:), allocatable :: tmp_table
          integer :: k

          if(this%nb_pts.eq.0) then
             if(.not.allocated(this%pts)) then
                allocate(this%pts(2,10))
             end if
          else
             if(size(this%pts,2)<this%nb_pts+1) then
                allocate(tmp_table(2,this%nb_pts+10))
                do k=1, this%nb_pts
                   tmp_table(1,k) = this%pts(1,k)
                   tmp_table(2,k) = this%pts(2,k)
                end do
                call MOVE_ALLOC(tmp_table,this%pts)
             end if
          end if

          this%nb_pts = this%nb_pts+1

          this%pts(1,this%nb_pts) = point(1)
          this%pts(2,this%nb_pts) = point(2)

        end subroutine add_pt_to_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reallocate the pts attribute to the nb_pts size
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine minimize_path(this)

          implicit none

          class(bf_path_icr)         , intent(inout) :: this
          
          integer(ikind), dimension(:,:), allocatable :: tmp_table
          integer :: k          

          if(this%nb_pts.ne.0) then
             if(this%nb_pts.ne.size(this%pts,2)) then
                allocate(tmp_table(2,this%nb_pts))
                do k=1, this%nb_pts
                   tmp_table(1,k) = this%pts(1,k)
                   tmp_table(2,k) = this%pts(2,k)
                end do
                call MOVE_ALLOC(tmp_table,this%pts)
             end if
          end if

        end subroutine minimize_path
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the bc_interior_pt passed as argument
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !--------------------------------------------------------------
        subroutine analyze_pt(
     $     this,
     $     bc_interior_pt_analyzed,
     $     interface_used)

          implicit none

          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(bf_interface)         , intent(in)    :: interface_used


          integer :: tolerance_pts_same_path
          integer :: tolerance_pts_same_sublayer


          !initialization of the tolerance
          tolerance_pts_same_path     = 2*bc_size+1
          tolerance_pts_same_sublayer = bc_size+1


          !if this is the first point in the path,
          !it should be added to the path
          if(this%nb_pts.eq.0) then

             !< analyze the point considering that
             !> it is the first point in the path
             call this%analyze_path_first_pt(
     $            bc_interior_pt_analyzed,
     $            interface_used,
     $            tolerance_pts_same_sublayer)

          !if this is the second or more, the pt
          !should be analyzed whether it is part
          !of the same path or to be added to
          !another path
          else

             !analyze the point considering that
             !it is the second or more point in
             !the path
             call this%analyze_path_next_pt(
     $            bc_interior_pt_analyzed,
     $            interface_used,
     $            tolerance_pts_same_path,
     $            tolerance_pts_same_sublayer)

          end if

        end subroutine analyze_pt
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> process the data encapsulated in the path to update
        !> the buffer layers
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !--------------------------------------------------------------
        subroutine process_path(
     $     this,
     $     interface_used,
     $     interior_nodes,
     $     dx, dy)
        
          implicit none
        
          class(bf_path_icr)              , intent(inout) :: this
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy

        
          type(bf_sublayer), pointer :: modified_sublayer


          !update the allocation required for the buffer layer
          modified_sublayer => update_allocation_bf_sublayer(
     $         this,
     $         interface_used,
     $         interior_nodes,
     $         dx, dy)

          !update the grid points for the increase
          call interface_used%update_grdpts_after_increase(
     $         modified_sublayer, this%pts(1:2,1:this%nb_pts))

          !reinitialize the path
          call this%reinitialize()

        end subroutine process_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze a bc_interior_pt that could be
        !> consider as the first point in the path
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !>@param tolerance_pts_same_sublayer
        !> tolerance to decide whether a grid point should belong
        !> to an existing sublayer
        !--------------------------------------------------------------        
        subroutine analyze_path_first_pt(
     $     this,
     $     bc_interior_pt_analyzed,
     $     interface_used,
     $     tolerance_pts_same_sublayer)

          implicit none

          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(bf_interface)         , intent(in)    :: interface_used
          integer                     , intent(in)    :: tolerance_pts_same_sublayer

          
          integer(ikind), dimension(2) :: local_coord


          !get the mainlayer id of the path
          this%mainlayer_id = interface_used%get_mainlayer_id(
     $         bc_interior_pt_analyzed)


          !check if there is already an exiting buffer layer in which
          !this bc_interior_pt could fit: if there is a matching layer,
          !the pointer is associated to the sublayer matching the general
          !coordinates. Otherwise, the pointer is nullified
          this%matching_sublayer => interface_used%get_sublayer(
     $         bc_interior_pt_analyzed, local_coord,
     $         tolerance_i=tolerance_pts_same_sublayer,
     $         mainlayer_id_i=this%mainlayer_id)


          !initialize the alignment with the position of the
          !bc_interior_pt
          this%alignment(1,1) = bc_interior_pt_analyzed(1)
          this%alignment(1,2) = bc_interior_pt_analyzed(1)
          this%alignment(2,1) = bc_interior_pt_analyzed(2)
          this%alignment(2,2) = bc_interior_pt_analyzed(2)


          !add the bc_interior_pt_analyzed to the current path
          call this%add_pt_to_path(bc_interior_pt_analyzed)

        end subroutine analyze_path_first_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze a bc_interior_pt that should be
        !> consider as the next point in the path. The path
        !> is not empty
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !>@param tolerance_pts_same_path
        !> tolerance to decide whether a grid point should belong
        !> to the current path
        !
        !>@param tolerance_pts_same_sublayer
        !> tolerance to decide whether a grid point should belong
        !> to an existing sublayer
        !--------------------------------------------------------------        
        subroutine analyze_path_next_pt(
     $     this,
     $     bc_interior_pt_analyzed,
     $     interface_used,
     $     tolerance_pts_same_path,
     $     tolerance_pts_same_sublayer)

          implicit none

          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(bf_interface)         , intent(in)    :: interface_used
          integer                     , intent(in)    :: tolerance_pts_same_path
          integer                     , intent(in)    :: tolerance_pts_same_sublayer

          
          integer(ikind), dimension(2) :: local_coord
          logical :: same_mainlayer
          logical :: pts_far_away
          type(bf_sublayer), pointer :: sublayer_new_pt


          !check whether the bc_interior_pt analyzed belongs to the same
          !mainlayer as the previous path
          same_mainlayer = this%are_pts_in_same_mainlayer(
     $         bc_interior_pt_analyzed,
     $         interface_used)


          !if they are in the same layer, we check whether the points are
          !too far away and are still part of the same sublayer
          if(same_mainlayer) then

             !check if there is already an exiting buffer layer associated
             !with the current path
             if(associated(this%matching_sublayer)) then
                
                
                !check if the current point is not far away from the
                !previous point analyzed
                pts_far_away = this%are_pts_far_away(
     $               bc_interior_pt_analyzed,
     $               tolerance_pts_same_path)
             
             
                !if the points are too far away, we need to check
                !whether they could belong to the same sublayer
                if(pts_far_away) then
                   
                   sublayer_new_pt => interface_used%get_sublayer(
     $                  bc_interior_pt_analyzed,
     $                  local_coord,
     $                  tolerance_i=tolerance_pts_same_sublayer,
     $                  mainlayer_id_i=this%mainlayer_id)
             
                   !if the points are too far away and cannot belong
                   !to the same buffer layer afterwards, the path
                   !should be ended
                   if(.not.associated(sublayer_new_pt, this%matching_sublayer)) then
                      this%ends = .true.
                   end if
             
                end if
             
             !if no sublayer is already matching the current path,
             !we need to check whether the current point belongs to
             !the same path as the previous point and if it could be
             !part of an existing sublayer
             else
             
                !check if the current point is not far away from the
                !previous point analyzed
                pts_far_away = this%are_pts_far_away(
     $               bc_interior_pt_analyzed,
     $               tolerance_pts_same_path)
             
             
                !if both points do not belong to the same path, the
                !current path should be ended
                if(pts_far_away) then
                   this%ends=.true.
             
             
                !if both points belongs to the same path, we should
                !check if this new point could match an existing layer
                !and so if the current path could eventually be part
                !of an existing layer
                else                
                   this%matching_sublayer => interface_used%get_sublayer(
     $                  bc_interior_pt_analyzed,
     $                  local_coord,
     $                  tolerance_i=tolerance_pts_same_sublayer,
     $                  mainlayer_id_i=this%mainlayer_id)
             
                end if
             
             end if

          !if the bc_interior_pt analyzed is not part of the same 
          !mainlayer, we need to cut the path
          else
             this%ends=.true.
          end if


          !if the path is not ended then the path should be
          !updated with the new bc_interior_pt and so it 
          !should update also the alignment
          if(.not.this%ends) then
             call this%update_alignment(bc_interior_pt_analyzed)
             call this%add_pt_to_path(bc_interior_pt_analyzed)
          end if

        end subroutine analyze_path_next_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the bc_interior_pt analyzed is located in
        !> in the same buffer main layer as the grid points stored
        !> in the path
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !--------------------------------------------------------------
        function are_pts_in_same_mainlayer(
     $     this,
     $     bc_interior_pt_analyzed,
     $     interface_used)

          implicit none

          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(bf_interface)         , intent(in)    :: interface_used
          logical                                     :: are_pts_in_same_mainlayer

          integer :: mainlayer_current_pt
          

          !compute the mainlayer of the current point
          mainlayer_current_pt = interface_used%get_mainlayer_id(
     $         bc_interior_pt_analyzed)

          !compare with the mainlayer of the path
          are_pts_in_same_mainlayer =
     $         this%mainlayer_id.eq.mainlayer_current_pt


        end function are_pts_in_same_mainlayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze whether the bc_interior_pt point is far away
        !> from the last point added to the path
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param tolerance
        !> tolerance deciding the maximum distance between the last
        !> point fo teh path and the new point so that the point
        !> belongs to the same path
        !--------------------------------------------------------------
        function are_pts_far_away(
     $     this, bc_interior_pt_analyzed, tolerance)

          implicit none

          class(bf_path_icr)        , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: bc_interior_pt_analyzed
          integer                     , intent(in) :: tolerance
          logical                                  :: are_pts_far_away

          are_pts_far_away = 
     $         (abs(this%pts(1,this%nb_pts) -
     $         bc_interior_pt_analyzed(1)).gt.tolerance).or.
     $         (abs(this%pts(2,this%nb_pts) -
     $         bc_interior_pt_analyzed(2)).gt.tolerance)

        end function are_pts_far_away


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the two buffer layers in clockwise direction
        !> should be merged
        !
        !> @warning bf_layer2%alignment(direction,max) >
        !> this%alignment(direction,min)
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bf_sublayer_ptr
        !> reference to the second bf_sublayer merged. The first one is
        !> the matching bf_sublayer of the path
        !
        !>@param merge
        !> logical stating whether the two bf_sublayers should be
        !> merged
        !--------------------------------------------------------------
        function should_bf_sublayers_be_merged(
     $     this, bf_sublayer_ptr)
     $     result(merge)


          implicit none

          class(bf_path_icr), intent(in) :: this
          type(bf_sublayer)   , intent(in) :: bf_sublayer_ptr
          logical                          :: merge

          integer :: direction
          
          !determine the direction and the border compared
          !between the path and the neighboring layer of the
          !matching buffer layer of the path
          select case(this%mainlayer_id)
            case(N,S)
               direction = x_direction
            case(E,W)
               direction = y_direction
            case default
               call error_mainlayer_id(
     $              'bf_path_icr_class.f',
     $              'should_bf_sublayers_be_merged',
     $              this%mainlayer_id)
          end select          


          !if the distance between the alignment of the path
          !and the neighboring buffer layer of its matching
          !layer is small, the two buffer layers will be merged
          merge = (this%alignment(direction,max_border)).gt.(
     $         bf_sublayer_ptr%get_alignment(direction,min_border)-2*bc_size)

        end function should_bf_sublayers_be_merged


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the allocation of the buffer layers to add
        !> new grid points
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !>@param interior_nodes
        !> table encapsulating the data of the grid points of the
        !> interior domain
        !
        !>@return modified_sublayer
        !> reference to the bf_sublayer whose size is updated
        !--------------------------------------------------------------
        function update_allocation_bf_sublayer(
     $     this, interface_used, interior_nodes, dx, dy)
     $     result(modified_sublayer)

          implicit none

          class(bf_path_icr)              , intent(inout) :: this
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                     , intent(in)    :: dx
          real(rkind)                     , intent(in)    :: dy
          type(bf_sublayer), pointer                      :: modified_sublayer
          
          type(bf_sublayer), pointer     :: neighboring_sublayer
          logical                        :: merge


          !does the current path have a matching sublayer ?
          !if it has one, then the matching sublayer will be reallocated
          if(associated(this%matching_sublayer)) then

             !we need to check whether the reallocation of the matching
             !sublayer with the alignment of the current path will not
             !lead to a merge with the neighboring buffer layer
             neighboring_sublayer => this%matching_sublayer%get_next()

             if(associated(neighboring_sublayer)) then
                merge = this%should_bf_sublayers_be_merged(
     $               neighboring_sublayer)

                if(merge) then
                   modified_sublayer => interface_used%merge_sublayers(
     $                  this%matching_sublayer,
     $                  neighboring_sublayer,
     $                  interior_nodes,
     $                  this%alignment)                   
                else

                   call update_alignment_for_reallocation(this)
                   
                   call interface_used%reallocate_sublayer(
     $                  this%matching_sublayer,
     $                  interior_nodes,
     $                  this%alignment)
                   modified_sublayer => this%matching_sublayer

                end if

             !if there is no potential neighboring sublayer that should be
             !merged, the matching sublayer is simply reallocated according
             !the alignment of the current path
             else
                call update_alignment_for_reallocation(this)

                call interface_used%reallocate_sublayer(
     $               this%matching_sublayer,
     $               interior_nodes,
     $               this%alignment)
                modified_sublayer => this%matching_sublayer

             end if

          !if there is no matching sublayer, a new buffer layer will be
          !allocated with the needed new alignment
          else
             modified_sublayer => interface_used%allocate_sublayer(
     $            this%mainlayer_id,
     $            interior_nodes,
     $            this%alignment,
     $            dx,dy)

          end if
          
        end function update_allocation_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the alignment such that the reallocation process can
        !> only increase the buffer layer sizes
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine update_alignment_for_reallocation(this)

          implicit none

          class(bf_path_icr), intent(inout) :: this

          integer(ikind), dimension(2,2) :: p_alignment

          p_alignment = this%matching_sublayer%get_alignment_tab()

          this%alignment(1,1) = min(p_alignment(1,1),this%alignment(1,1))
          this%alignment(2,1) = min(p_alignment(2,1),this%alignment(2,1))
          this%alignment(1,2) = max(p_alignment(1,2),this%alignment(1,2))
          this%alignment(2,2) = max(p_alignment(2,2),this%alignment(2,2))

        end subroutine update_alignment_for_reallocation


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the alignment of the path. This data is
        !> used to decide whether the corresponding buffer layer
        !> will be allocated or reallocated and how its position
        !> may vary compared to the interior domain
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !--------------------------------------------------------------
        subroutine update_alignment(this, bc_interior_pt_analyzed)
        
          implicit none
        
          class(bf_path_icr)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
        
          this%alignment(1,1) = min(this%alignment(1,1), bc_interior_pt_analyzed(1))
          this%alignment(1,2) = max(this%alignment(1,2), bc_interior_pt_analyzed(1))
          this%alignment(2,1) = min(this%alignment(2,1), bc_interior_pt_analyzed(2))
          this%alignment(2,2) = max(this%alignment(2,2), bc_interior_pt_analyzed(2))
             
        end subroutine update_alignment


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_pts attribute
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@return get_nb_pts
        !> nb_pts attribute
        !--------------------------------------------------------------
        function get_nb_pts(this)

          implicit none

          class(bf_path_icr), intent(in) :: this
          integer                          :: get_nb_pts

          get_nb_pts = this%nb_pts

        end function get_nb_pts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the general coordinates of the ith grid point
        !> stored in the path
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param i
        !> index identifying the grid point stored in the path
        !
        !>@return get_pt
        !> general coordinates of the ith grid point
        !--------------------------------------------------------------
        function get_pt(this,i)

          implicit none

          class(bf_path_icr), intent(in) :: this
          integer             , intent(in) :: i
          integer(ikind), dimension(2)     :: get_pt

          get_pt(1) = this%pts(1,i)
          get_pt(2) = this%pts(2,i)

        end function get_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the path on screen
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_path_icr object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine print_on_screen(this)

          implicit none

          class(bf_path_icr), intent(in) :: this

          print '(''path created'')'
          print '(''nb_pts: '', I2)', this%nb_pts
          print '(''alignment: '', 4I3)', this%alignment
          if(associated(this%matching_sublayer)) then
             print '(''sublayer_matching: associated'')'
             print '(''  + localization: '', I2)',
     $            this%matching_sublayer%get_localization()
             print '(''  + alignment: '', 4I3)',
     $            this%matching_sublayer%get_alignment_tab()
          else
             print '(''sublayer_matching: non associated'')'
          end if
          print '('''')'
              
        end subroutine print_on_screen

      end module bf_path_icr_class
