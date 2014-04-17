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
      ! 16_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_path_class

        use bf_sublayer_class       , only : bf_sublayer
        use interface_abstract_class, only : interface_abstract
        use parameters_input        , only : nx,ny,bc_size
        use parameters_constant     , only : S_W,S_E,N_W,N_E
        use parameters_kind         , only : ikind

        implicit none

        private
        public :: bf_layer_path


        !> @class bf_layer_path
        !> class encapsulating data required to decide which points
        !> need to be updated and whether this leads to the creation
        !> or the update of a buffer layer
        !
        !> @param leads_to_new_bf_layer
        !> logical indicating whether a new buffer layer should be
        !> created or not
        !
        !> @param matching_layer
        !> pointer to the buffer layer that should be updated
        !> if the pointer is not associated, it means that the
        !> buffer layer corresponding to this path should be
        !> allocated
        !
        !> @param ends
        !> logical identifying whether the path should ends or not
        !> and so whether the path should be interpreted to allocate
        !> and update an existing buffer layer
        !
        !> @param ends_with_corner
        !> logical identifying whether the path ended with a corner
        !> this fact will determine whether the new buffer layer
        !> should be allocated or updated with exchanging gridpoints
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
        !> @param neighbors
        !> table of logical identifying whether neighboring buffer layers
        !> exist and so how to allocate new buffer layers
        !
        !> @param alignment
        !> table identifying the new position of the buffer layer and
        !> whether an existing buffer layer should be reallocated
        !---------------------------------------------------------------
        type :: bf_layer_path

          type(bf_sublayer), pointer :: matching_sublayer

          logical :: ends
          logical :: ends_with_corner
          integer :: corner_id

          integer(ikind), dimension(:,:), allocatable :: pts
          integer                                     :: nb_pts
          logical       , dimension(4)                :: neighbors
          integer(ikind), dimension(2,2)              :: alignment

          contains
          
          procedure, pass          :: ini
          procedure, pass          :: reinitialize
          procedure, pass, private :: add_pt_to_path
          procedure, pass, private :: minimize_path

          procedure, pass, private :: are_pts_far_away
          procedure, pass, private :: update_alignment

          procedure, pass, private :: analyze_path_first_pt
          procedure, pass, private :: analyze_path_next_pt
          procedure, pass          :: analyze_pt

          !procedure, pass          :: process_path

        end type bf_layer_path


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(bf_layer_path), intent(inout) :: this

          nullify(this%matching_sublayer)
          
          this%ends = .false.
          this%ends_with_corner = .false.

          this%nb_pts = 0
          this%neighbors = [.false.,.false.,.false.,.false.]

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !--------------------------------------------------------------
        subroutine reinitialize(this)
        
          implicit none

          class(bf_layer_path), intent(inout) :: this

          nullify(this%matching_sublayer)
          
          this%ends = .false.
          this%ends_with_corner = .false.

          this%nb_pts = 0
          this%neighbors = [.false.,.false.,.false.,.false.]

        end subroutine reinitialize


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the general coordinates of a point
        !> to the current path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param point
        !> general coordinates of the grid point added to the path
        !--------------------------------------------------------------
        subroutine add_pt_to_path(this, point)

          implicit none

          class(bf_layer_path)        , intent(inout) :: this
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
        !> subroutine adding the general coordinates of a point
        !> to the current path
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param point
        !> general coordinates of the grid point added to the path
        !--------------------------------------------------------------
        subroutine minimize_path(this)

          implicit none

          class(bf_layer_path)         , intent(inout) :: this
          
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
        !> subroutine analyzing the bc_interior_pt given
        !
        !> @date
        !> 16_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> interface_abstract class encapsulating the pointers
        !> to the buffer main layers
        !--------------------------------------------------------------
        subroutine analyze_pt(
     $     this,
     $     bc_interior_pt_analyzed,
     $     interface_used)

          implicit none

          class(bf_layer_path)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(interface_abstract)   , intent(in)    :: interface_used


          integer :: corner_id
          integer :: tolerance_pts_same_path
          integer :: tolerance_pts_same_sublayer


          !initialization of the tolerance
          tolerance_pts_same_path     = 2*bc_size+1
          tolerance_pts_same_sublayer = bc_size+1


          !ask whether the pt to be added is a corner pt
          corner_id = is_pt_a_corner_pt(bc_interior_pt_analyzed)

          
          !if it is a corner pt, the path immediately ends
          if(corner_id.ne.0) then
             this%ends = .true.
             this%ends_with_corner = .true.
             this%corner_id=corner_id


          !otherwise, the pt is analyzed to check whether
          !it should be part of the same path or added to
          !another path
          else
             
             !if this is the first point in the path,
             !it should be added to the path
             if(this%nb_pts.eq.0) then

                !< analyze the point considering that
                !> it is the first point in the path
                call this%analyze_path_first_pt(
     $               bc_interior_pt_analyzed,
     $               interface_used,
     $               tolerance_pts_same_sublayer)

             !if this is the second or more, the pt
             !should be analyzed whether it is part
             !of the same path or to be added to
             !another path
             else

                !< analyze the point considering that
                !> it is the second or more point in
                !> the path
                call this%analyze_path_next_pt(
     $               bc_interior_pt_analyzed,
     $               interface_used,
     $               tolerance_pts_same_path,
     $               tolerance_pts_same_sublayer)

             end if

          end if             

        end subroutine analyze_pt
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine analyzing the bc_interior_pt given
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> buffer layer path object encapsulating data determining
        !> whether a buffer layer should be allocated or updated
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> interface_abstract class encapsulating the pointers
        !> to the buffer main layers
        !--------------------------------------------------------------
        !subroutine process_path(
     $  !   this,
     $  !   bc_interior_pt_analyzed,
     $  !   interface_used)
        !
        !  implicit none
        !
        !  class(bf_layer_path)        , intent(inout) :: this
        !  integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
        !  class(interface_abstract)   , intent(in)    :: interface_used
        !
        !
        !  !< only if the path ends, we should consider updating the
        !  !> buffer layers corresponding to the path
        !  if(this%ends) then
        !
        !     !< if the path stoped at a corner, then we should
        !     !> consider allocating the corner buffer layer
        !
        !  end if
        !
        !
        !end subroutine process_path


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine analyzing whether the bc_interior_pt point
        !> analyzed is one of the interior corner of the interior
        !> domain
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_interior_pt
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !--------------------------------------------------------------
        function is_pt_a_corner_pt(bc_interior_pt) result(corner_id)
        
          implicit none

          integer(ikind), dimension(2), intent(in)  :: bc_interior_pt
          integer                                   :: corner_id


          if(bc_interior_pt(1).eq.bc_size) then
             if(bc_interior_pt(2).eq.bc_size) then
                corner_id=S_W
             else
                if(bc_interior_pt(2).eq.(ny-1)) then
                   corner_id=N_W
                else
                   corner_id=0
                end if
             end if

          else
             if(bc_interior_pt(1).eq.(nx-1)) then
                if(bc_interior_pt(2).eq.bc_size) then
                   corner_id=S_E
                else
                   if(bc_interior_pt(2).eq.(ny-1)) then
                      corner_id=S_W
                   else
                      corner_id=0
                   end if
                end if
             else
                corner_id=0
             end if
          end if

        end function is_pt_a_corner_pt        


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine analyzing a bc_interior_pt that should be
        !> consider as the next point in the path. The path
        !> is empty
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> interface_abstract class encapsulating the pointers
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

          class(bf_layer_path)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(interface_abstract)   , intent(in)    :: interface_used
          integer                     , intent(in)    :: tolerance_pts_same_sublayer

          
          integer(ikind), dimension(2) :: local_coord


          !< check if there is already an exiting buffer layer in which
          !> this bc_interior_pt could fit: if there is a matching layer,
          !> the pointer is associated to the sublayer matching the general
          !> coordinates. Otherwise, the pointer is nullified
          this%matching_sublayer => interface_used%get_sublayer(
     $         bc_interior_pt_analyzed, local_coord,
     $         tolerance_pts_same_sublayer)


          !< initialize the alignment with the position of the
          !> bc_interior_pt
          this%alignment(1,1) = bc_interior_pt_analyzed(1)
          this%alignment(1,2) = bc_interior_pt_analyzed(1)
          this%alignment(2,1) = bc_interior_pt_analyzed(2)
          this%alignment(2,2) = bc_interior_pt_analyzed(2)


          !< add the bc_interior_pt_analyzed to the current path
          call this%add_pt_to_path(bc_interior_pt_analyzed)

        end subroutine analyze_path_first_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine analyzing a bc_interior_pt that should be
        !> consider as the next point in the path. The path
        !> is not empty
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> interface_abstract class encapsulating the pointers
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

          class(bf_layer_path)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          class(interface_abstract)   , intent(in)    :: interface_used
          integer                     , intent(in)    :: tolerance_pts_same_path
          integer                     , intent(in)    :: tolerance_pts_same_sublayer

          
          integer(ikind), dimension(2) :: local_coord
          logical :: pts_far_away
          type(bf_sublayer), pointer :: sublayer_new_pt


          !< check if there is already an exiting buffer layer associated
          !> with the current path
          if(associated(this%matching_sublayer)) then
             
             
             !< check if the current point is not far away from the
             !> previous point analyzed
             pts_far_away = this%are_pts_far_away(
     $            bc_interior_pt_analyzed,
     $            tolerance_pts_same_path)

             !< if the points are too far away, we need to check
             !> whether they could belong to the same sublayer
             if(pts_far_away) then
                
                sublayer_new_pt => interface_used%get_sublayer(
     $               bc_interior_pt_analyzed,
     $               local_coord,
     $               tolerance_pts_same_sublayer)              

                !< if the points are too far away and cannot belong
                !> to the same buffer layer afterwards, the path
                !> should be ended
                if(.not.associated(sublayer_new_pt, this%matching_sublayer)) then
                   this%ends = .true.
                end if

             end if

          !< if no sublayer is already matching the current path,
          !> we need to check whether the current point belongs to
          !> the same path as the previous point and if it could be
          !> part of an existing sublayer
          else

             !< check if the current point is not far away from the
             !> previous point analyzed
             pts_far_away = this%are_pts_far_away(
     $            bc_interior_pt_analyzed,
     $            tolerance_pts_same_path)

             !< if both points do not belong to the same path, the
             !> current path should be ended
             if(pts_far_away) then
                this%ends=.true.

             !< if both points belongs to the same path, we should
             !> check if this new point could match an existing layer
             !> and so if the current path could eventually be part
             !> of an existing layer
             else                
                this%matching_sublayer => interface_used%get_sublayer(
     $               bc_interior_pt_analyzed,
     $               local_coord,
     $               tolerance_pts_same_sublayer)

             end if

          end if


          !< if the path is not ended then the path should be
          !> updated with the new bc_interior_pt and so it 
          !> should update also the alignment
          if(.not.this%ends) then
             call this%update_alignment(bc_interior_pt_analyzed)
             call this%add_pt_to_path(bc_interior_pt_analyzed)
          end if

        end subroutine analyze_path_next_pt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine analyzing whether the bc_interior_pt point
        !> analyzed is one of the interior corner of the interior
        !> domain
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !
        !>@param interface_used
        !> interface_abstract class encapsulating the pointers
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
        function are_pts_far_away(this, bc_interior_pt_analyzed, tolerance)

          implicit none

          class(bf_layer_path)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
          integer                     , intent(in)    :: tolerance
          logical                                     :: are_pts_far_away

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
        !> subroutine updating the alignment of the path, that is
        !> used to decide whether the corresponding buffer layer
        !> will be allocated or reallocated and how its position
        !> may vary compared to the interior domain
        !
        !> @date
        !> 17_04_2013 - initial version - J.L. Desmarais
        !
        !>@param bc_interior_pt_analyzed
        !> integer, dimension(2) with the general coords of the
        !> bc_interior_pt analyzed
        !--------------------------------------------------------------        
        subroutine update_alignment(this, bc_interior_pt_analyzed)
        
          implicit none
        
          class(bf_layer_path)        , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: bc_interior_pt_analyzed
        
          this%alignment(1,1) = min(this%alignment(1,1), bc_interior_pt_analyzed(1))
          this%alignment(1,2) = max(this%alignment(1,2), bc_interior_pt_analyzed(1))
          this%alignment(2,1) = min(this%alignment(2,1), bc_interior_pt_analyzed(2))
          this%alignment(2,2) = max(this%alignment(2,2), bc_interior_pt_analyzed(2))
             
        end subroutine update_alignment

      end module bf_layer_path_class
