      !> @file
      !> object gathering the update operations to be applied on
      !> one buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> object gathering the update operations to be applied on
      !> one buffer layer
      !
      !> @date
      !> 29_10_2014 - initial version - J.L. Desmarais
      !> 20_03_2015 - update          - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_path_class

        use bf_layer_errors_module, only :
     $       error_mainlayer_id
       
        use bf_layer_exchange_module, only :
     $       update_alignment_for_exchanges

        use bf_increase_coords_module, only :
     $       get_mainlayer_coord

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        
        private
        public :: icr_path


        !> @class icr_pth
        !> object gathering the update operations to be applied on
        !> one buffer layer
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
        !>@param share_update_operations_with
        !> check whether the two paths have update operations
        !> in common
        !
        !>@param update_alignment
        !> update the alignment of the path. This data is
        !> used to decide whether the corresponding buffer layer
        !> will be allocated or reallocated and how its position
        !> may vary compared to the interior domain
        !
        !>@param are_pts_in_same_mainlayer
        !> check whether the bc_interior_pt analyzed is located in
        !> in the same buffer main layer as the grid points stored
        !> in the path
        !
        !>@param are_pts_in_same_path
        !> analyze whether the bc_interior_pt point analyzed should
        !> be contained in the same path
        !
        !>@param stage_first_grdpt_for_update
        !> submit the first grid-point to the staging area
        !
        !>@param stage_next_grdpt_for_update
        !> submit a grid-point to the staging area. This grid-point
        !> is not the first of the path
        !
        !>@param stage
        !> submit a grid-point to the staging area
        !
        !>@param should_bf_layers_be_merged
        !> check whether the two buffer layers should be merged
        !
        !>@param update_allocation_bf_layer
        !> update the allocation of the buffer layer to have
        !> memory for the new grid points
        !
        !>@param commit
        !> process the data encapsulated in the path to update
        !> the buffer layers
        !
        !>@param merge
        !> merge the update operations from two paths
        !
        !>@param remove
        !> remove the data stored in the path
        !---------------------------------------------------------------
        type :: icr_path

          type(bf_sublayer), pointer :: matching_sublayer

          integer                                     :: mainlayer_id
          integer(ikind), dimension(:,:), allocatable :: pts
          integer                                     :: nb_pts
          integer(ikind), dimension(2,2)              :: alignment
          logical                                     :: ends

          contains
          
          ! procedures for basic interactions with
          ! the main attributes
          procedure, pass :: ini
          procedure, pass :: reinitialize
          procedure, pass :: is_ended
          procedure, pass :: add_pt_to_path
          procedure, pass :: share_update_operations_with

          ! procedure for storing activated bc_interior_pt
          ! in the path, storing changes to be applied on
          ! the buffer layer
          procedure, pass :: update_alignment
          procedure, pass :: are_pts_in_same_mainlayer
          procedure, pass :: are_pts_in_same_path
          procedure, pass :: stage_first_grdpt_for_update
          procedure, pass :: stage_next_grdpt_for_update
          procedure, pass :: stage

          ! procedure for applying changes on the buffer
          ! layer
          procedure, pass :: should_bf_layers_be_merged
          procedure, pass :: update_allocation_bf_layer
          procedure, pass :: commit

          ! merge and finalization procedures
          procedure, pass :: merge
          procedure, pass :: remove

        end type icr_path


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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(icr_path), intent(inout) :: this

          nullify(this%matching_sublayer)
          
          this%ends   = .false.
          this%nb_pts = 0

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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !--------------------------------------------------------------
        subroutine reinitialize(this)
        
          implicit none

          class(icr_path), intent(inout) :: this

          call ini(this)

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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@return is_ended
        !> ends attribute
        !--------------------------------------------------------------
        function is_ended(this)

          implicit none

          class(icr_path), intent(in) :: this
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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param point
        !> general coordinates of the grid point added to the path
        !--------------------------------------------------------------
        subroutine add_pt_to_path(this, gen_coords)

          implicit none

          class(icr_path)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: gen_coords
          
          integer(ikind), dimension(:,:), allocatable :: tmp_table

          if(this%nb_pts.eq.0) then
             if(.not.allocated(this%pts)) then
                allocate(this%pts(2,10))
             end if
          else
             if(size(this%pts,2).lt.(this%nb_pts+1)) then
                allocate(tmp_table(2,this%nb_pts+10))
                tmp_table(:,1:this%nb_pts) = this%pts
                call MOVE_ALLOC(tmp_table,this%pts)
             end if
          end if

          this%nb_pts = this%nb_pts+1

          this%pts(:,this%nb_pts) = gen_coords

        end subroutine add_pt_to_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the two paths point to the same
        !> buffer layer to perform the update operations
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param other_path
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !--------------------------------------------------------------
        function share_update_operations_with(this, other_path)
     $     result(share_updates)

          implicit none

          class(icr_path), intent(in) :: this
          class(icr_path), intent(in) :: other_path
          logical                     :: share_updates

          integer            :: dir
          integer, parameter :: tolerance = bc_size


          if((this%nb_pts.gt.0).and.(other_path%nb_pts.gt.0)) then
             
             if(this%mainlayer_id.eq.other_path%mainlayer_id) then

                select case(this%mainlayer_id)
                  case(N,S)
                     dir = x_direction
                  case(E,W)
                     dir = y_direction
                  case default
                     call error_mainlayer_id(
     $                    'icr_path_class.f',
     $                    'share_update_operations_with',
     $                    this%mainlayer_id)
                end select

                share_updates = (
     $               min(this%alignment(dir,2)+bc_size,
     $                   other_path%alignment(dir,2)+bc_size+tolerance)-
     $               max(this%alignment(dir,1)-bc_size,
     $                   other_path%alignment(dir,1)-bc_size-tolerance)+1).gt.0
             else
                share_updates = .false.
             end if

          else
             share_updates = .false.
          end if
        
        end function share_update_operations_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the alignment of the path to match the new grid-point
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param bf_alignment
        !> integer, dimension(2,2) with the general coords of the
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
        subroutine update_alignment(this,gen_coords,alignment)

          implicit none

          class(icr_path)                         , intent(inout) :: this
          integer(ikind), dimension(2)  , optional, intent(in)    :: gen_coords
          integer(ikind), dimension(2,2), optional, intent(in)    :: alignment

          
          if(present(gen_coords)) then
             this%alignment(1,1) = min(this%alignment(1,1),gen_coords(1))
             this%alignment(1,2) = max(this%alignment(1,2),gen_coords(1))
             this%alignment(2,1) = min(this%alignment(2,1),gen_coords(2))
             this%alignment(2,2) = max(this%alignment(2,2),gen_coords(2))
          end if


          if(present(alignment)) then
             this%alignment(1,1) = min(this%alignment(1,1),alignment(1,1))
             this%alignment(1,2) = max(this%alignment(1,2),alignment(1,2))
             this%alignment(2,1) = min(this%alignment(2,1),alignment(2,1))
             this%alignment(2,2) = max(this%alignment(2,2),alignment(2,2))
          end if


          call update_alignment_for_exchanges(
     $         this%mainlayer_id,
     $         this%alignment)

        end subroutine update_alignment


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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param gen_coords
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !--------------------------------------------------------------
        function are_pts_in_same_mainlayer(
     $     this,
     $     gen_coords)

          implicit none

          class(icr_path)             , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: gen_coords
          logical                                  :: are_pts_in_same_mainlayer

          integer :: mainlayer_current_pt
          

          !compute the mainlayer of the current point
          mainlayer_current_pt = get_mainlayer_coord(gen_coords)

          !compare with the mainlayer of the path
          are_pts_in_same_mainlayer = this%mainlayer_id.eq.mainlayer_current_pt

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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param gen_coords
        !> integer, dimension(2) with the general coordinates of the
        !> bc_interior_pt analyzed
        !
        !>@param tolerance
        !> tolerance deciding the maximum distance between the last
        !> point fo teh path and the new point so that the point
        !> belongs to the same path
        !--------------------------------------------------------------
        function are_pts_in_same_path(this, gen_coords, tolerance)

          implicit none

          class(icr_path)             , intent(in) :: this
          integer(ikind), dimension(2), intent(in) :: gen_coords
          integer                     , intent(in) :: tolerance
          logical                                  :: are_pts_in_same_path

          integer :: dir

          select case(this%mainlayer_id)
            case(N,S)
               dir = x_direction
            case(E,W)
               dir = y_direction
            case default
               call error_mainlayer_id(
     $              'icr_path_class',
     $              'are_pts_in_same_path',
     $              this%mainlayer_id)
          end select

          are_pts_in_same_path = (
     $         min(this%alignment(dir,2),gen_coords(dir)+tolerance)-
     $         max(this%alignment(dir,1),gen_coords(dir)-tolerance)+1).gt.0

        end function are_pts_in_same_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> stage a bc_interior_pt for update. This grid-point
        !> is the first point stored in the path
        !
        !> @date
        !> 19_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param gen_coords
        !> coordinates of the bc_interior_pt analyzed in the general
        !> reference frame
        !
        !>@param interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !>@param tolerance_pts_same_sublayer
        !> tolerance to decide whether a grid point should belong
        !> to an existing sublayer
        !--------------------------------------------------------------        
        subroutine stage_first_grdpt_for_update(
     $     this,
     $     gen_coords,
     $     bf_interface_used,
     $     tolerance_pts_same_sublayer)

          implicit none

          class(icr_path)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: gen_coords
          class(bf_interface_coords)  , intent(in)    :: bf_interface_used
          integer                     , intent(in)    :: tolerance_pts_same_sublayer


          !get the mainlayer id of the path
          this%mainlayer_id = get_mainlayer_coord(gen_coords)


          !check if there is already an exiting buffer layer in which
          !this bc_interior_pt could fit: if there is a matching layer,
          !the pointer is associated to the sublayer matching the general
          !coordinates. Otherwise, the pointer is nullified
          this%matching_sublayer => bf_interface_used%get_bf_layer_from_gen_coords(
     $         gen_coords,
     $         tolerance=tolerance_pts_same_sublayer,
     $         mainlayer_id=this%mainlayer_id)


          !initialize the alignment with the position of the
          !bc_interior_pt
          this%alignment(1,1) = gen_coords(1)
          this%alignment(1,2) = gen_coords(1)
          this%alignment(2,1) = gen_coords(2)
          this%alignment(2,2) = gen_coords(2)


          !update the alignment to take into account potential
          !alignment adjustments for the buffer layers at the
          !interface between main layers
          if(associated(this%matching_sublayer)) then
             call update_alignment(
     $            this,
     $            alignment=this%matching_sublayer%get_alignment_tab())
          else
             call update_alignment(this)
          end if


          !add the bc_interior_pt_analyzed to the current path
          call this%add_pt_to_path(gen_coords)

        end subroutine stage_first_grdpt_for_update


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
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param gen_coords
        !> general coordinates of the bc_interior_pt analyzed
        !
        !>@param bf_interface_used
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
        subroutine stage_next_grdpt_for_update(
     $     this,
     $     gen_coords,
     $     bf_interface_used,
     $     tolerance_pts_same_path,
     $     tolerance_pts_same_sublayer)

          implicit none

          class(icr_path)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: gen_coords
          class(bf_interface_coords)  , intent(in)    :: bf_interface_used
          integer                     , intent(in)    :: tolerance_pts_same_path
          integer                     , intent(in)    :: tolerance_pts_same_sublayer

          logical :: same_mainlayer
          logical :: pts_in_same_path


          !initialize the boolean for the state of the path
          this%ends = .false.


          !check whether the bc_interior_pt analyzed belongs to the same
          !mainlayer as the previous path
          same_mainlayer = this%are_pts_in_same_mainlayer(gen_coords)


          !if the grid-points are not in the same mainlayer,
          !the path is ended
          if(.not.same_mainlayer) then

             this%ends = .true.

          !if they are in the same layer, we check whether
          !the points are too far away
          else

             pts_in_same_path = this%are_pts_in_same_path(
     $            gen_coords,
     $            tolerance_pts_same_path)

             if(.not.pts_in_same_path) then
                
                this%ends = .true.

             else

                !check if there is already an existing buffer layer
                !associated with the current path
                if(.not.associated(this%matching_sublayer)) then
                                   
                   this%matching_sublayer => bf_interface_used%get_bf_layer_from_gen_coords(
     $                  gen_coords,
     $                  tolerance=tolerance_pts_same_sublayer,
     $                  mainlayer_id=this%mainlayer_id)
             
                   !if the new grid-point can be placed in this sublayer
                   !the alignment of the path is updated with the alignment
                   !of the new matching sublayer
                   if(associated(this%matching_sublayer)) then
                      call update_alignment(
     $                     this,
     $                     alignment=this%matching_sublayer%get_alignment_tab())
                   end if
             
                end if

             end if

          end if


          !if the path is not ended then the path should be
          !updated with the new bc_interior_pt and so it 
          !should update also the alignment
          if(.not.this%ends) then
             call this%update_alignment(gen_coords=gen_coords)
             call this%add_pt_to_path(gen_coords)
          end if

        end subroutine stage_next_grdpt_for_update


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze a bc_interior_pt and add stage it for the
        !> the update operations
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param gen_coords
        !> general coordinates of the bc_interior_pt analyzed
        !
        !>@param bf_interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !--------------------------------------------------------------  
        subroutine stage(
     $     this,
     $     gen_coords,
     $     bf_interface_used)

          implicit none

          class(icr_path)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: gen_coords
          class(bf_interface_coords)  , intent(in)    :: bf_interface_used


          integer, parameter :: tolerance_pts_same_sublayer = 2*bc_size
          integer, parameter :: tolerance_pts_same_path     = 3*bc_size


          if(this%nb_pts.eq.0) then
             
             call stage_first_grdpt_for_update(
     $            this,
     $            gen_coords,
     $            bf_interface_used,
     $            tolerance_pts_same_sublayer)

          else

             call stage_next_grdpt_for_update(
     $            this,
     $            gen_coords,
     $            bf_interface_used,
     $            tolerance_pts_same_path,
     $            tolerance_pts_same_sublayer)

          end if

        end subroutine stage


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether two neighboring buffer layers should be merged
        !
        !> @warning bf_layer2%alignment(direction,max) >
        !> this%alignment(direction,min)
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param bf_sublayer_ptr
        !> reference to the second bf_sublayer merged. The first one is
        !> the matching bf_sublayer of the path
        !
        !>@param merge
        !> logical stating whether the two bf_sublayers should be
        !> merged
        !--------------------------------------------------------------
        function should_bf_layers_be_merged(
     $     this, bf_sublayer_ptr)
     $     result(merge)


          implicit none

          class(icr_path)  , intent(in) :: this
          type(bf_sublayer), intent(in) :: bf_sublayer_ptr
          logical                       :: merge


          integer, parameter :: tolerance = bc_size
          integer            :: direction

          
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
     $              'icr_path_class',
     $              'should_bf_layers_be_merged',
     $              this%mainlayer_id)
          end select          


          !if the distance between the alignment of the path
          !and the neighboring buffer layer of its matching
          !layer is small, the two buffer layers will be merged
          merge = (min(this%alignment(direction,2)+bc_size,
     $                 bf_sublayer_ptr%get_alignment(direction,2)+bc_size+tolerance)-
     $             max(this%alignment(direction,1)-bc_size,
     $                 bf_sublayer_ptr%get_alignment(direction,1)-bc_size-tolerance)+1).gt.0

        end function should_bf_layers_be_merged


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the allocation of the buffer layer to add
        !> new grid points
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param bf_interface_used
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
        function update_allocation_bf_layer(
     $     this,
     $     bf_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes)
     $     result(modified_sublayer)

          implicit none

          class(icr_path)                 , intent(inout) :: this
          class(bf_interface_coords)      , intent(inout) :: bf_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(bf_sublayer), pointer                      :: modified_sublayer
          
          type(bf_sublayer), pointer :: neighboring_sublayer
          logical                    :: merge


          !does the current path have a matching sublayer ?
          !if it has one, then the matching sublayer will be reallocated
          if(associated(this%matching_sublayer)) then

             !we need to check whether the reallocation of the matching
             !sublayer with the alignment of the current path will not
             !lead to a merge with the neighboring buffer layer
             neighboring_sublayer => this%matching_sublayer%get_next()

             if(associated(neighboring_sublayer)) then
                merge = this%should_bf_layers_be_merged(
     $               neighboring_sublayer)
             else
                merge = .false.
             end if


             if(merge) then
                modified_sublayer => bf_interface_used%merge_sublayers(
     $               this%matching_sublayer,
     $               neighboring_sublayer,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               this%alignment)
             else
                
                call bf_interface_used%reallocate_sublayer(
     $               this%matching_sublayer,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               this%alignment)
                modified_sublayer => this%matching_sublayer
                
             end if

          !if there is no matching sublayer, a new buffer layer will be
          !allocated with the needed new alignment
          else
             modified_sublayer => bf_interface_used%allocate_sublayer(
     $            this%mainlayer_id,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            this%alignment)

          end if
          
        end function update_allocation_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask the path to perform the update operations for
        !> the update of the buffer layer
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param bf_interface_used
        !> bf_interface object encapsulating the references
        !> to the buffer main layers
        !
        !> @param p_model
        !> physical model
        !
        !> @param t
        !> time 
        !
        !> @param dt
        !> time step
        !
        !> @param interior_x_map
        !> array with x-coordinates of the interior domain
        !
        !> @param interior_y_map
        !> array with y-coordinates of the interior domain
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !--------------------------------------------------------------  
        subroutine commit(
     $     this,
     $     bf_interface_used,
     $     p_model,
     $     t,
     $     dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(icr_path)                    , intent(inout) :: this
          class(bf_interface_coords)         , intent(inout) :: bf_interface_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)                        , intent(in)    :: dt
          real(rkind)   , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes1


          type(bf_sublayer), pointer :: bf_sublayer_ptr


          if(this%nb_pts.gt.0) then

             !1) update the allocation of the buffer layer
             bf_sublayer_ptr => update_allocation_bf_layer(
     $            this,
     $            bf_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes1)
             
          
             !2) update the configuration of the grid points
             !   in the buffer layer
             call bf_interface_used%update_grdpts_id_in_bf_layer(
     $            p_model,
     $            t,
     $            dt,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            bf_sublayer_ptr,
     $            this%pts(:,1:this%nb_pts))

             !3) reinitialize the path
             call reinitialize(this)

          end if

        end subroutine commit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask the path to perform the update operations for
        !> the update of the buffer layer
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !>@param other_path
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !
        !> @param check_for_merge
        !> optional argument controlling if the two paths can be
        !> merged, .false. by default
        !--------------------------------------------------------------  
        subroutine merge(this, other_path, check_for_merge)

          implicit none

          class(icr_path)  , intent(inout) :: this
          class(icr_path)  , intent(in)    :: other_path
          logical, optional, intent(in)    :: check_for_merge

          logical        :: check_for_merge_op
          logical        :: merge_paths
          integer(ikind) :: n1
          integer(ikind) :: n2
          integer(ikind), dimension(:,:), allocatable :: tmp_pts


          if(present(check_for_merge)) then
             check_for_merge_op = check_for_merge
          else
             check_for_merge_op = .true.
          end if

          if(check_for_merge_op) then
             merge_paths = this%share_update_operations_with(other_path)
          else
             merge_paths = .true.
          end if


          if(merge_paths) then

             !update the alignment
             call this%update_alignment(alignment=other_path%alignment)

             !gather the pts
             n1 = this%nb_pts
             n2 = other_path%nb_pts

             if((n1+n2).gt.size(this%pts,2)) then

                allocate(tmp_pts(2,n1+n2))

                tmp_pts(:,1:n1)       = this%pts(:,1:n1)
                tmp_pts(:,n1+1:n1+n2) = other_path%pts(:,1:n2)

                call MOVE_ALLOC(tmp_pts,this%pts)

             else
                
                this%pts(:,n1+1:n1+n2) = other_path%pts(:,1:n2)

             end if
             this%nb_pts = n1+n2
             
             !gather the matching sublayer
             if(.not.associated(this%matching_sublayer).and.
     $            associated(other_path%matching_sublayer)) then

                this%matching_sublayer => other_path%matching_sublayer

             end if

          end if

        end subroutine merge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the data stored in the path
        !
        !> @date
        !> 20_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object gathering the update operations to be applied on
        !> one buffer layer
        !--------------------------------------------------------------  
        subroutine remove(this)

          implicit none

          class(icr_path), intent(inout) :: this

          if(allocated(this%pts)) then
             deallocate(this%pts)
          end if

        end subroutine remove
          
      end module icr_path_class
