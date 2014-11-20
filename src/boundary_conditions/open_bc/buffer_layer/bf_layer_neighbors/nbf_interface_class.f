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

        use bf_interior_bc_sections_module, only :
     $       ini_interior_bc_sections,
     $       close_last_bc_section,
     $       set_full_interior_bc_section,
     $       minimize_interior_bc_section

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_layer_newgrdpt_procedure_module, only :
     $       get_newgrdpt_procedure,
     $       get_interior_data_for_newgrdpt,
     $       are_intermediate_newgrdpt_data_needed,
     $       get_x_map_for_newgrdpt,
     $       get_y_map_for_newgrdpt
      
        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_element_class, only :
     $       nbf_element

        use nbf_list_class, only :
     $       nbf_list

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size,debug
        
        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sbf_list_class, only :
     $       sbf_list

        implicit none


        private
        public :: nbf_interface

        
        !>@class nbf_interface
        !> object encapsulting links to buffer layers at the edge between
        !> different main layers
        !
        !>@param nbf_links
        !> array of pointers to the buffer layers at the edge of different
        !> main layers
        !> ex: nbf_links(N,1) : links to the buffer layers considered
        !>                      neighbors of type 1 by the north main
        !>                      layer
        !>     nbf_links(N,2) : links to the buffer layers considered
        !>                      neighbors of type 2 by the north main
        !>                      layer
        !
        !>@param ini
        !> initialize the array of links by initializing the
        !> nbf_list containing the links
        !
        !>@param link_neighbor1_to_bf_sublayer
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 1 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !>@param link_neighbor2_to_bf_sublayer
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 2 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !>@param update_link_from_neighbor1_to_bf_sublayer
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 1 by nbf_sublayer1 have their links updated
        !
        !>@param update_link_from_neighbor2_to_bf_sublayer
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 2 by nbf_sublayer1 have their links updated
        !
        !>@param remove_link_from_neighbor1_to_bf_sublayer
        !> remove the links existing from neighbor1 to nbf_sublayer
        !
        !>@param remove_link_from_neighbor2_to_bf_sublayer
        !> remove the links existing from neighbor2 to nbf_sublayer
        !
        !>@param update_grdpts_from_neighbors
        !> ask all the buffer layers that have grid points in common
        !> with the current main layer to send data to the current
        !> buffer layer
        !
        !>@param update_neighbor_grdpts
        !> ask all the main layers that have grid points in common
        !> with the current main layer to receive data from the current
        !> buffer layer
        !
        !>@param sync_interface_nodes
        !> synchronize the nodes located at the interface between
        !> buffer main layers
        !
        !>@param define_integration_borders
        !> define the x_borders,y_borders and N_bc_sections and
        !> S_bc_sections for the buffer layer depending on its
        !> neighboring buffer layers
        !
        !>@param update_neighbors_integration_borders
        !> ask the buffer layers neighboring the recently modified
        !> buffer layer to update thier integration borders if they
        !> effectively share grid points with the current buffer
        !> layer
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
        !
        !>@param print_on_screen
        !> print the links between bf_sublayers on screen
        !--------------------------------------------------------------
        type nbf_interface

          type(nbf_list), dimension(4,2), private :: nbf_links

          contains

          procedure, pass :: ini

          procedure, pass :: link_neighbor1_to_bf_sublayer
          procedure, pass :: link_neighbor2_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor2_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor2_to_bf_sublayer

          procedure, pass :: update_grdpts_from_neighbors
          procedure, pass :: update_neighbor_grdpts

          procedure, pass :: sync_interface_nodes
          procedure, pass :: define_integration_borders
          procedure, pass :: update_neighbors_integration_borders
          procedure, pass :: update_integration_borders

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

          procedure, pass :: compute_newgrdpt

          procedure, pass :: print_on_screen

        end type nbf_interface


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the array of links by initializing the
        !> nbf_list containing the links
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
        subroutine ini(this)
        
          implicit none

          class(nbf_interface), intent(inout) :: this

          integer :: i,j

          do j=1, size(this%nbf_links,2)
             do i=1, size(this%nbf_links,1)
                call this%nbf_links(i,j)%ini()
             end do
          end do

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 1 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------        
        subroutine link_neighbor1_to_bf_sublayer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%add_link_in_list(
     $         nbf_sublayer)

        end subroutine link_neighbor1_to_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 2 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------        
        subroutine link_neighbor2_to_bf_sublayer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%add_link_in_list(
     $         nbf_sublayer)

        end subroutine link_neighbor2_to_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 1 by nbf_sublayer1 have their links updated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer1
        !> bf_sublayer reference before the update
        !
        !>@param nbf_sublayer2
        !> bf_sublayer reference after the update
        !--------------------------------------------------------------
        subroutine update_link_from_neighbor1_to_bf_sublayer(
     $     this, nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer1%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_link_from_neighbor1_to_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 2 by nbf_sublayer1 have their links updated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer1
        !> bf_sublayer reference before the update
        !
        !>@param nbf_sublayer2
        !> bf_sublayer reference after the update
        !--------------------------------------------------------------
        subroutine update_link_from_neighbor2_to_bf_sublayer(
     $     this, nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer1%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_link_from_neighbor2_to_bf_sublayer
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the links existing from neighbor1 to nbf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer reference removed
        !--------------------------------------------------------------
        subroutine remove_link_from_neighbor1_to_bf_sublayer(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%remove_link_from_list(
     $         nbf_sublayer)

        end subroutine remove_link_from_neighbor1_to_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the links existing from neighbor2 to nbf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer reference removed
        !--------------------------------------------------------------
        subroutine remove_link_from_neighbor2_to_bf_sublayer(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%remove_link_from_list(
     $         nbf_sublayer)

        end subroutine remove_link_from_neighbor2_to_bf_sublayer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask all the buffer layers that have grid points in common
        !> with the current main layer to send data to the current
        !> buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer object updated from its references
        !--------------------------------------------------------------
        subroutine update_grdpts_from_neighbors(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: nbf_sublayer

          integer :: mainlayer_id
          
          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor1()) then

             call this%nbf_links(mainlayer_id,1)%copy_from_neighbors_to_bf_layer(
     $            1, nbf_sublayer)
          end if
          

          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_from_neighbors_to_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_grdpts_from_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask all the main layers that have grid points in common
        !> with the current main layer to receive data from the current
        !> buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer object updated from its references
        !--------------------------------------------------------------
        subroutine update_neighbor_grdpts(this, nbf_sublayer)

          implicit none

          class(nbf_interface), intent(inout) :: this
          type(bf_sublayer)   , intent(in)    :: nbf_sublayer

          integer :: mainlayer_id          


          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor1()) then
             call this%nbf_links(mainlayer_id,1)%copy_to_neighbors_from_bf_layer(
     $            1, nbf_sublayer)
          end if


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_to_neighbors_from_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_neighbor_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes located at the interface between
        !> buffer main layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
        subroutine sync_interface_nodes(this)

          implicit none

          class(nbf_interface), intent(inout) :: this


          !synchronize the nodes at the SW interface
          !-----------------------------------------
          !nbf_links(W,1) are references to buffer
          !layers of the S layer that may have grid
          !points in common with the W main layer
          !-----------------------------------------
          !The W layer is a neighbor of type 1 for
          !the S layer
          !-----------------------------------------
          call this%nbf_links(W,1)%sync_nodes_with_neighbor1(
     $         this%nbf_links(S,1))


          !synchronize the nodes at the SE interface
          !-----------------------------------------
          !nbf_links(E,1) are references to buffer
          !layers of the S layer that may have grid
          !points in common with the E main layer
          !-----------------------------------------
          !The E layer is a neighbor of type 2 for
          !the S layer
          !-----------------------------------------
          call this%nbf_links(E,1)%sync_nodes_with_neighbor2(
     $         this%nbf_links(S,2))


          !synchronize the nodes at the NW interface
          !-----------------------------------------
          !nbf_links(W,2) are references to buffer
          !layers of the N layer that may have grid
          !points in common with the W main layer
          !-----------------------------------------
          !The W layer is a neighbor of type 1 for
          !the N layer
          !-----------------------------------------
          call this%nbf_links(W,2)%sync_nodes_with_neighbor1(
     $         this%nbf_links(N,1))


          !synchronize the nodes at the NE interface
          !-----------------------------------------
          !nbf_links(E,2) are references to buffer
          !layers of the N layer that may have grid
          !points in common with the E main layer
          !-----------------------------------------
          !The E layer is a neighbor of type 2 for
          !the N layer
          !-----------------------------------------
          call this%nbf_links(E,2)%sync_nodes_with_neighbor2(
     $         this%nbf_links(N,2))

        end subroutine sync_interface_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes located at the interface between
        !> buffer main layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
        subroutine define_integration_borders(this,buffer_layer)

          implicit none

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: buffer_layer

          integer                                     :: localization
          integer(ikind), dimension(2)                :: sizes
          integer(ikind), dimension(2)                :: x_borders
          integer(ikind), dimension(2)                :: y_borders
          integer(ikind), dimension(:,:), allocatable :: bc_sections


          localization = buffer_layer%get_localization()
          sizes        = buffer_layer%get_sizes()


          select case(localization)
      
            case(N)
               call buffer_layer%set_x_borders([1,sizes(1)])
               call buffer_layer%set_y_borders([bc_size+1,sizes(2)])

               if(buffer_layer%can_exchange_with_neighbor1().or.
     $              buffer_layer%can_exchange_with_neighbor2()) then

                  call determine_bc_sections_for_NorS(
     $                 buffer_layer,
     $                 this%nbf_links(N,1),
     $                 this%nbf_links(N,2),
     $                 bc_sections)

                  if(allocated(bc_sections)) then
                     call buffer_layer%set_S_bc_sections(bc_sections)
                  else
                     call buffer_layer%remove_S_bc_sections()
                  end if

               else
                  call buffer_layer%remove_S_bc_sections()                  
               end if

            case(S)
               call buffer_layer%set_x_borders([1,sizes(1)])
               call buffer_layer%set_y_borders([1,sizes(2)-bc_size])

               if(buffer_layer%can_exchange_with_neighbor1().or.
     $              buffer_layer%can_exchange_with_neighbor2()) then

                  call determine_bc_sections_for_NorS(
     $                 buffer_layer,
     $                 this%nbf_links(S,1),
     $                 this%nbf_links(S,2),
     $                 bc_sections)

                  if(allocated(bc_sections)) then
                     call buffer_layer%set_N_bc_sections(bc_sections)
                  else
                     call buffer_layer%remove_N_bc_sections()
                  end if

               else
                  call buffer_layer%remove_N_bc_sections()
               end if

            case(E,W)
               if(localization.eq.E) then
                  x_borders = [bc_size+1,sizes(1)]
               else
                  x_borders = [1,sizes(1)-bc_size]
               end if
               y_borders = [bc_size+1,sizes(2)-bc_size]

               
               !if the buffer layer exchanges with the south
               !buffer layer, the south boundary sections should
               !be determined
               if(buffer_layer%can_exchange_with_neighbor1()) then

                  call determine_bc_sections_for_EorW(
     $                 localization,
     $                 buffer_layer,
     $                 this%nbf_links(localization,1),
     $                 bc_sections)

                  if(allocated(bc_sections)) then
                     call buffer_layer%set_S_bc_sections(bc_sections)
                  else
                     call buffer_layer%remove_S_bc_sections()
                  end if
                  
               !if the buffer layer does not exchange with the
               !south buffer layer, the south boundary layer
               !should be computed by the buffer layer
               else
                  y_borders(1) = 1
                  call buffer_layer%remove_S_bc_sections()
               end if


               !if the buffer layer exchanges with the north
               !buffer layer, the north boundary sections should
               !be determined
               if(buffer_layer%can_exchange_with_neighbor2()) then

                  call determine_bc_sections_for_EorW(
     $                 localization,
     $                 buffer_layer,
     $                 this%nbf_links(localization,2),
     $                 bc_sections)

                  if(allocated(bc_sections)) then
                     call buffer_layer%set_N_bc_sections(bc_sections)
                  else
                     call buffer_layer%remove_N_bc_sections()
                  end if

               !if the buffer layer does not exchange with the
               !north buffer layer, the south boundary layer
               !should be computed by the buffer layer
               else                
                  y_borders(2) = sizes(2)
                  call buffer_layer%remove_N_bc_sections()
               end if

               call buffer_layer%set_x_borders(x_borders)
               call buffer_layer%set_y_borders(y_borders)


            case default
               call error_mainlayer_id(
     $              'nbf_interface_class.f',
     $              'define_integration_borders',
     $              localization)

          end select

        end subroutine define_integration_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the bc_sections for a north or south
        !> buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !> @param buffer_layer
        !> buffer layer whose bc_sections are determined
        !
        !> @param bf_neighbors1
        !> buffer layer of type neighbors1 that can potentially
        !> exchange grid points with the buffer layer
        !
        !> @param bf_neighbors2
        !> buffer layer of type neighbors2 that can potentially
        !> exchange grid points with the buffer layer
        !
        !> @param bc_sections
        !> bc_sections for the buffer layer
        !--------------------------------------------------------------
        subroutine determine_bc_sections_for_NorS(
     $     buffer_layer,
     $     bf_neighbors1,
     $     bf_neighbors2,
     $     bc_sections)

          implicit none

          type(bf_sublayer)                          , intent(inout) :: buffer_layer
          type(nbf_list)                             , intent(in)    :: bf_neighbors1
          type(nbf_list)                             , intent(in)    :: bf_neighbors2
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: bc_sections


          integer(ikind) :: bf_inf
          integer(ikind) :: bf_sup

          integer(ikind) :: interior_inf
          integer(ikind) :: interior_sup

          integer :: nb_bc_sections
          logical :: min_initialized
          logical :: max_initialized
          logical :: no_bf_common_with_bf_layer


          !0) determine the borders identifying the
          !x-coordinates of the SW and SE corners
          !of the buffer layer
          bf_inf = buffer_layer%get_alignment(1,1)-bc_size
          bf_sup = buffer_layer%get_alignment(1,2)+bc_size

          
          !1) initialize the bc_sections
          call ini_interior_bc_sections(
     $         nb_bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_bf_layer)


          !2) check whether the neighbors of type 1 (W)
          !   have grid points in common with the
          !   current buffer layer
          !   (to prevent the inclusion of the interior
          !   boundary layer in the bc_sections, the
          !   interior_sup is initialized at min(0,bf_sup))
          if(buffer_layer%can_exchange_with_neighbor1()) then

             if(bf_inf.le.0) then

                interior_inf = bf_inf
                interior_sup = min(0,bf_sup)
                
                call bf_neighbors1%update_bc_sections(
     $               interior_inf,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               max_initialized,
     $               no_bf_common_with_bf_layer)
                
             end if
          end if


          !3) check whether the neighbors of type 2 (E)
          !   have grid points in common with the
          !   current buffer layer
          !   (to prevent the inclusion of the interior
          !   boundary layer in the bc_sections, the
          !   interior_inf is initialized at max(bf_inf,nx+1))
          if(buffer_layer%can_exchange_with_neighbor2()) then

             if(bf_sup.ge.(nx+1)) then
                interior_inf = max(bf_inf,nx+1)
                interior_sup = bf_sup

                call bf_neighbors2%update_bc_sections(
     $               interior_inf,
     $               interior_sup,
     $               nb_bc_sections,
     $               bc_sections,
     $               min_initialized,
     $               max_initialized,
     $               no_bf_common_with_bf_layer)

             end if
          end if


          !4) close the last bc_section
          call close_last_bc_section(
     $         nb_bc_sections,
     $         bc_sections,
     $         interior_sup,
     $         min_initialized,
     $         max_initialized)


          !5) check if there were no buffer layers
          !   in common with the buffer layer to
          !   initialize the bc_sections with the
          !   entire boundary layer
          if(no_bf_common_with_bf_layer) then
             
             nb_bc_sections = 0

             if(bf_inf.le.0) then
                nb_bc_sections = nb_bc_sections+1
             end if

             if(bf_sup.ge.(nx+1)) then
                nb_bc_sections = nb_bc_sections+1
             end if

             if(nb_bc_sections.gt.0) then
                allocate(bc_sections(2,nb_bc_sections))

                if(bf_inf.le.0) then
                   bc_sections(:,1) =
     $                  [bf_inf,min(bf_sup,0)]
                end if

                if(bf_sup.ge.(nx+1)) then
                   bc_sections(:,size(bc_sections,2)) =
     $                  [max(bf_inf,nx+1),bf_sup]
                end if
             end if

          else

          !6) adapt the size of the bc_sections
             call minimize_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections)
             
          end if

          !7) convert the general cooordinates saved
          !   in the bc_sections into local coordinates
          !   for the buffer layer
          call convert_bc_sections_to_local_coords(
     $         nb_bc_sections,
     $         bc_sections,
     $         bf_inf)

        end subroutine determine_bc_sections_for_NorS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the bc_sections for a north or south
        !> buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !> @param bf_localization
        !> cardinal coordinate identifying the buffer main layer
        !> containing the buffer layer, buffer_layer
        !
        !> @param buffer_layer
        !> buffer layer whose bc_sections are determined
        !
        !> @param bf_neighbors
        !> buffer layer of type neighbors that can potentially
        !> exchange grid points with the buffer layer
        !
        !> @param bc_sections
        !> bc_sections for the buffer layer
        !--------------------------------------------------------------
        subroutine determine_bc_sections_for_EorW(
     $     bf_localization,
     $     buffer_layer,
     $     bf_neighbors,
     $     bc_sections)

          implicit none

          integer                                    , intent(in)    :: bf_localization
          type(bf_sublayer)                          , intent(inout) :: buffer_layer
          type(nbf_list)                             , intent(in)    :: bf_neighbors
          integer(ikind), dimension(:,:), allocatable, intent(out)   :: bc_sections


          integer(ikind) :: bf_inf
          integer(ikind) :: bf_sup

          integer(ikind) :: interior_inf
          integer(ikind) :: interior_sup

          integer :: nb_bc_sections
          logical :: min_initialized
          logical :: max_initialized
          logical :: no_bf_common_with_bf_layer


          !0) determine the borders identifying the
          !x-coordinates of the SW and SE corners
          !of the buffer layer
          bf_inf = buffer_layer%get_alignment(1,1)-bc_size
          bf_sup = buffer_layer%get_alignment(1,2)+bc_size


          !1) initialize the bc_sections
          call ini_interior_bc_sections(
     $         nb_bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_bf_layer)

          !2) initialize the boundary of the domain
          !   studied
          if(bf_localization.eq.W) then
             
             interior_inf = bf_inf
             interior_sup = bc_size

          else

             interior_inf = nx-1
             interior_sup = bf_sup

          end if

          !3) determine the bc_sections if any
          call bf_neighbors%update_bc_sections(
     $         interior_inf,
     $         interior_sup,
     $         nb_bc_sections,
     $         bc_sections,
     $         min_initialized,
     $         max_initialized,
     $         no_bf_common_with_bf_layer)
          
          !4) close the last bc_section
          call close_last_bc_section(
     $         nb_bc_sections,
     $         bc_sections,
     $         interior_sup,
     $         min_initialized,
     $         max_initialized)

          !5) check if there were no buffer layers
          !   in common with the buffer layer to
          !   initialize the bc_sections with the
          !   entire boundary layer
          if(no_bf_common_with_bf_layer) then

             call set_full_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections,
     $            interior_inf,
     $            interior_sup)

          !6) adapt the size of the bc_sections
          else

             call minimize_interior_bc_section(
     $            nb_bc_sections,
     $            bc_sections)

          end if

          !7) convert the general cooordinates saved
          !   in the bc_sections into local coordinates
          !   for the buffer layer
          call convert_bc_sections_to_local_coords(
     $         nb_bc_sections,
     $         bc_sections,
     $         bf_inf)          

        end subroutine determine_bc_sections_for_EorW


        subroutine convert_bc_sections_to_local_coords(
     $     nb_bc_sections,
     $     bc_sections,
     $     bf_inf)

          implicit none

          integer                                    , intent(in)    :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          integer(ikind)                             , intent(in)    :: bf_inf

          integer :: k

          if(nb_bc_sections.gt.0) then

             do k=1, nb_bc_sections
                bc_sections(1,k) = bc_sections(1,k)-(bf_inf-1)
                bc_sections(2,k) = bc_sections(2,k)-(bf_inf-1)
             end do

          end if

        end subroutine convert_bc_sections_to_local_coords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask the buffer layers sharing grid points with the current
        !> buffer layer to update the integration borders of the layer
        !> (x_borders, y_borders, N_bc_sections, S_bc_sections)
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param buffer_layer
        !> buffer layer which has been recently modified and whose
        !> modifications have a consequence on the neighboring buffer
        !> layers
        !--------------------------------------------------------------
        subroutine update_neighbors_integration_borders(this,buffer_layer)

          implicit none

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: buffer_layer

          integer                    :: localization
          type(nbf_element), pointer :: nbf_current_ele
          integer                    :: k
          type(bf_sublayer), pointer :: nbf_sublayer


          localization = buffer_layer%get_localization()
          

          !if the buffer layer can potentially share grid points with
          !its neighboring buffer layer of type 1, the neighboring
          !buffer layers of type 1 are updated only if they effectively
          !share grid points with the current buffer layer
          !----------------------------------------------------------
          if(buffer_layer%can_exchange_with_neighbor1()) then

             nbf_current_ele => this%nbf_links(localization,1)%get_head()

             do k=1,  this%nbf_links(localization,1)%get_nb_elements()

                if(nbf_current_ele%shares_grdpts_along_x_dir_with(buffer_layer)) then
                   nbf_sublayer => nbf_current_ele%get_ptr()

                   call define_integration_borders(this, nbf_sublayer)

                end if

                nbf_current_ele => nbf_current_ele%get_next()
             end do
          end if

          !if the buffer layer can potentially share grid points with
          !its neighboring buffer layer of type 2, the neighboring
          !buffer layers of type 2 are updated only if they effectively
          !share grid points with the current buffer layer
          !----------------------------------------------------------
          if(buffer_layer%can_exchange_with_neighbor2()) then

             nbf_current_ele => this%nbf_links(localization,2)%get_head()

             do k=1,  this%nbf_links(localization,2)%get_nb_elements()

                if(nbf_current_ele%shares_grdpts_along_x_dir_with(buffer_layer)) then
                   nbf_sublayer => nbf_current_ele%get_ptr()

                   call define_integration_borders(this, nbf_sublayer)

                end if

                nbf_current_ele => nbf_current_ele%get_next()
             end do
          end if

        end subroutine update_neighbors_integration_borders


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the integration borders of the current buffer layer
        !> and ask the neighboring buffer layers that are effectively
        !> sharing grid points with the current buffer layer to update
        !> their integration borders
        !
        !> @date
        !> 31_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param buffer_layer
        !> buffer layer recently modified whose integration borders
        !> are determined and whose neighboring buffer layers effectively
        !> sharing grid points with the current buffer layer are asked to
        !> update their integration borders
        !--------------------------------------------------------------
        subroutine update_integration_borders(this,buffer_layer)

          implicit none

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: buffer_layer


          call define_integration_borders(this,buffer_layer)
          call update_neighbors_integration_borders(this,buffer_layer)

        end subroutine update_integration_borders


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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the new grid point of the buffer layer
        !
        !> @date
        !> 19_11_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !> @param bf_sublayer_ptr
        !> buffer layer whose new grid point is computed
        !
        !> @param i1
        !> x-index of the new grid point
        !
        !> @param j1
        !> y-index of the new grid point
        !
        !> @param interior_nodes0
        !> array with the grid point data at t=t-dt
        !
        !> @param interior_nodes1
        !> array with the grid point data at t=t
        !--------------------------------------------------------------
        subroutine compute_newgrdpt(
     $     this,
     $     bf_sublayer_ptr,
     $     p_model,t,dt,
     $     i1,j1,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1)

          implicit none

          class(nbf_interface)            , intent(in)    :: this
          type(bf_sublayer)               , intent(inout) :: bf_sublayer_ptr
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          integer(ikind)                  , intent(in)    :: i1
          integer(ikind)                  , intent(in)    :: j1
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2)   :: gen_coords
          logical                        :: tmp_needed
          integer(ikind), dimension(2,2) :: gen_borders

          integer    , dimension(2*bc_size+1,2*bc_size+1)    :: tmp_grdpts_id0
          real(rkind), dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes0
          real(rkind), dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes1

          integer                                :: localization
          type(bf_newgrdpt)                      :: bf_newgrdpt_used
          integer                                :: procedure_type
          integer                                :: gradient_type
          integer(ikind), dimension(2,2)         :: bf_align
          real(rkind)   , dimension(2*bc_size+1) :: tmp_x_map
          real(rkind)   , dimension(2*bc_size+1) :: tmp_y_map
          real(rkind)   , dimension(ne)          :: new_grdpt
          integer                                :: k


          !localization of the buffer layer
          localization = bf_sublayer_ptr%get_localization()

          !compute the general coordinates of the new grid point
          match_table   = bf_sublayer_ptr%get_general_to_local_coord_tab()
          gen_coords(1) = match_table(1) + i1
          gen_coords(2) = match_table(2) + j1
          
          !determine the extent of the data needed
          gen_borders(1,1) = gen_coords(1)-bc_size
          gen_borders(1,2) = gen_coords(1)+bc_size
          gen_borders(2,1) = gen_coords(2)-bc_size
          gen_borders(2,2) = gen_coords(2)+bc_size

          !determine whether the new grid point is at the interface
          !between buffer layers or at the interface with the interior
          !and so requires to create intermediate arrays gathering data
          !from the interior and the neighboring buffer layers
          tmp_needed = are_intermediate_newgrdpt_data_needed(
     $         localization,
     $         gen_borders)

          if(.not.tmp_needed) then
             tmp_needed = .not.(bf_sublayer_ptr%does_previous_timestep_exist())
          end if


          !if intermediate data are required, data are gathered from
          !the interior, the buffer layer and its neighbors
          if(tmp_needed) then

             !gather the data from the interior domain
             call get_interior_data_for_newgrdpt(
     $            interior_nodes0,
     $            interior_nodes1,
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            gen_borders)

             !gather the data from the current buffer
             !layer
             call bf_sublayer_ptr%get_data_for_newgrdpt(
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            gen_borders)

             !gather the data from the potential neighbors
             !of the buffer layer
             !- gather data from neighbor of type 1
             if(gen_borders(1,1).le.0) then
                
                call this%nbf_links(localization,1)%get_data_for_newgrdpt(
     $               tmp_grdpts_id0,
     $               tmp_nodes0,
     $               tmp_nodes1,
     $               gen_borders)

             end if

             !- gather data from neighbor of type 2
             if(gen_borders(1,2).ge.(nx+1)) then

                call this%nbf_links(localization,2)%get_data_for_newgrdpt(
     $               tmp_grdpts_id0,
     $               tmp_nodes0,
     $               tmp_nodes1,
     $               gen_borders)

             end if

             
             !determine the procedure and gradient type
             call get_newgrdpt_procedure(
     $            bc_size+1,bc_size+1,
     $            tmp_grdpts_id0,
     $            procedure_type,
     $            gradient_type)


             !compute the new grdpt
             ! - construct the bf_align array to localize
             !   the buffer layer
             bf_align(1,1) = gen_coords(1)
             bf_align(1,2) = gen_coords(1)
             bf_align(2,1) = gen_coords(2)
             bf_align(2,2) = gen_coords(2)

             ! - construct the bf_x_map array
             tmp_x_map = get_x_map_for_newgrdpt(interior_x_map, gen_borders)

             ! - construct the bf_y_map array
             tmp_y_map = get_y_map_for_newgrdpt(interior_y_map, gen_borders)

             ! - compute the new grdpt
             new_grdpt = bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t, dt,
     $            bf_align, tmp_x_map, tmp_y_map, tmp_nodes0,
     $            bf_align, tmp_x_map, tmp_y_map, tmp_nodes1,
     $            bc_size+1, bc_size+1,
     $            procedure_type, gradient_type)


             !set the new grdpt in the buffer layer
             do k=1,ne
                call bf_sublayer_ptr%set_nodes_pt(i1,j1,k,new_grdpt(k))
             end do


          !otherwise, the grid point can be directly computed
          !from the buffer layer
          else

              call bf_sublayer_ptr%compute_newgrdpt(
     $            p_model,
     $            t,dt,
     $            i1,j1)

          end if

        end subroutine compute_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the links between bf_sublayers on screen
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
        subroutine print_on_screen(this)

          implicit none

          class(nbf_interface), intent(in) :: this

          integer     , dimension(4,2) :: neighbors
          character(1), dimension(4)   :: bf_layer_char
          integer                      :: i,j

          neighbors(N,1) = W
          neighbors(N,2) = E
          neighbors(S,1) = W
          neighbors(S,2) = E
          neighbors(E,1) = S
          neighbors(E,2) = N
          neighbors(W,1) = S
          neighbors(W,2) = N

          bf_layer_char = ['N','S','E','W']          

          do j=1, size(this%nbf_links,2)
             do i=1, size(this%nbf_links,1)
                print '(A1,'' --> '',A1)',
     $               bf_layer_char(neighbors(i,j)),
     $               bf_layer_char(i)
                call this%nbf_links(i,j)%print_on_screen()
             end do
          end do

        end subroutine print_on_screen

      end module nbf_interface_class
