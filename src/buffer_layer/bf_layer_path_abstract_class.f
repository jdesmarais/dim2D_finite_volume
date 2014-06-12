      !> @file
      !> module implementing the object encapsulating data
      !> required to decide which points need to be updated
      !> and whether this leads to the creation or the update
      !> of a buffer layer but its intrinsic procedures are
      !> not implemented
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating data
      !> required to decide which points need to be updated
      !> and whether this leads to the creation or the update
      !> of a buffer layer, but its intrinsic procedures are
      !> not implemented
      !
      !> @date
      ! 18_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_path_icr_abstract_class

        use bf_sublayer_class       , only : bf_sublayer
        use parameters_kind         , only : ikind

        implicit none

        private
        public :: bf_path_icr_abstract


        !> @class bf_path_icr_abstract
        !> class encapsulating data required to decide which points
        !> need to be updated and whether this leads to the creation
        !> or the update of a buffer layer, nut the intrinsic
        !> procedures are not implemented
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
        type :: bf_path_icr_abstract

          type(bf_sublayer), pointer :: matching_sublayer

          logical :: ends
          logical :: ends_with_corner
          integer :: corner_id

          integer                                     :: mainlayer
          integer(ikind), dimension(:,:), allocatable :: pts
          integer                                     :: nb_pts
          logical       , dimension(4)                :: neighbors
          integer(ikind), dimension(2,2)              :: alignment

        end type bf_path_icr_abstract

      end module bf_path_icr_abstract_class
