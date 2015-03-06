      !> @file
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_interface_class        

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_interface_newgrdpt_class, only :
     $       nbf_interface_newgrdpt

        use sbf_list_class, only :
     $       sbf_list

        implicit none


        private
        public :: nbf_interface

        
        !>@class nbf_interface
        !> nbf_interface_newgrdpt augmented with procedures to remove
        !> buffer layers
        !
        !>@param get_nbf_layers_sharing_grdpts_with
        !> add to the list of sublayer pointers the neighboring 
        !> buffer layers that shares grid points in the
        !> x-direction with the current buffer layer
        !
        !>@param bf_layer_depends_on_neighbors
        !> test whether the bf_sublayer is sharing grid points with
        !> its neighboring buffer layers
        !
        !>@param does_a_neighbor_remains
        !> test whether one of the bf_sublayer neighbors is remaining
        !--------------------------------------------------------------
        type, extends(nbf_interface_newgrdpt) :: nbf_interface

          contains

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

        end type nbf_interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add to the list of sublayer pointers the neighboring 
        !> buffer layers that shares grid points in the
        !> x-direction with the current buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer_i investigated
        !
        !>@param bf_sublayer_i
        !> reference to the bf_sublayer whose neighbors are investigated
        !
        !>@param bf_sublayer_list
        !> list of the bf_sublayer objects that share grid points in the
        !> x-direction with the current buffer layer
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !--------------------------------------------------------------
        subroutine get_nbf_layers_sharing_grdpts_with(
     $     this,
     $     nbf_type,
     $     bf_sublayer_i,
     $     bf_sublayer_list,
     $     bf_mainlayer_id)

          implicit none

          class(nbf_interface)      , intent(in)    :: this
          integer                   , intent(in)    :: nbf_type
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          type(sbf_list)            , intent(inout) :: bf_sublayer_list
          integer         , optional, intent(in)    :: bf_mainlayer_id


          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          call this%nbf_links(mainlayer_id,nbf_type)%get_nbf_layers_sharing_grdpts_with(
     $         bf_sublayer_i, bf_sublayer_list)

        end subroutine get_nbf_layers_sharing_grdpts_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> test whether the bf_sublayer is sharing grid points with
        !> its neighboring buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer investigated
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !
        !>@param dependent
        !> logical stating whether the buffer layer is sharing
        !> grid points with the neighboring buffer layers
        !--------------------------------------------------------------
        function bf_layer_depends_on_neighbors(
     $     this, nbf_type, bf_sublayer_i, bf_mainlayer_id)
     $     result(dependent)

          implicit none

          class(nbf_interface)      , intent(in) :: this
          integer                   , intent(in) :: nbf_type
          type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
          integer         , optional, intent(in) :: bf_mainlayer_id
          logical                                :: dependent

          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          dependent = this%nbf_links(mainlayer_id,nbf_type)%bf_layer_depends_on_neighbors(
     $         bf_sublayer_i)

        end function bf_layer_depends_on_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> test whether one of the bf_sublayer neighbors is remaining
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer investigated
        !
        !>@param bf_sublayer_id
        !> bf_sublayer 
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !
        !>@param a_neighbor_remains
        !> logical stating whether the buffer layer cannot be removed
        !> because a neighboring buffer layer should remain
        !--------------------------------------------------------------
        function does_a_neighbor_remains(
     $     this, nbf_type, bf_sublayer_i, bf_mainlayer_id)
     $     result(a_neighbor_remains)

          implicit none

          class(nbf_interface)      , intent(in)    :: this
          integer                   , intent(in)    :: nbf_type
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          integer         , optional, intent(in)    :: bf_mainlayer_id
          logical                                   :: a_neighbor_remains


          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          a_neighbor_remains = this%nbf_links(mainlayer_id,nbf_type)%does_a_neighbor_remains(
     $         bf_sublayer_i)

        end function does_a_neighbor_remains

      end module nbf_interface_class
