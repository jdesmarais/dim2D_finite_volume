      !> @file
      !> module implementing the object encapsulating the buffer layers
      !> around the interior domain and the relations to synchronize the
      !> data between them
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating the buffer layers
      !> around the interior domain and the relations to synchronize the
      !> data between them
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_class

        use netcdf

        use bc_operators_class, only :
     $       bc_operators
        
        use bf_interior_bc_sections_module, only :
     $       process_bc_sections_into_bc_procedure

        use bf_layer_errors_module, only :
     $       error_mainlayer_id,
     $       error_incompatible_neighbor

        use bf_restart_module, only :
     $       get_restart_alignment

        use bf_sublayer_class, only :
     $       bf_sublayer

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_mainlayer_pointer_class, only :
     $       bf_mainlayer_pointer

        use bf_layer_nf90_operators_module, only :
     $       print_interior_grdpts_id_on_netcdf

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use nbf_interface_newgrdpt_class, only :
     $       nbf_interface_newgrdpt

        use nf90_operators_module, only :
     $       nf90_close_file

        use nf90_operators_read_module, only :
     $       nf90_open_file_for_reading,
     $       nf90_get_varid,
     $       nf90_get_var_model_nopt,
     $       nf90_read_borders

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W,
     $       BF_SUCCESS

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction,
     $       interior,
     $       left,
     $       right

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       debug

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sbf_list_class, only :
     $       sbf_list

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators

        implicit none

        private
        public :: bf_interface
        

        !>@class bf_interface
        !> class encapsulating the buffer layers around the interior
        !> domain and subroutines to synchronize the data between them
        !
        !>@param mainlayer_pointers
        !> table with reference to the bf_mainlayers objects gathering
        !> the bf_sublayer corresponding to a specific cardinal point
        !
        !>@param border_interface
        !> object encapsulating references to the sublayers at the edge
        !> between the main layers and subroutines to synchronize the data
        !> between them
        !
        !>@param ini
        !> initialize the buffer/interior domain interface
        !
        !>@param get_mainlayer
        !> get the bf_mainlayer object corresponding to the cardinal point
        !
        !>@param allocate_sublayer
        !> allocate a bf_sublayer and insert it in the corresponding
        !> buffer mainlayer
        !
        !>@param reallocate_sublayer
        !> reallocate a buffer sublayer and check whether the neighboring
        !> buffer layer changed and so if the neighboring links should be
        !> updated
        !
        !>@param merge_sublayers
        !> merge the content of two sublayers and update
        !> the links to the border buffer layers
        !
        !>@param remove_sublayer
        !> remove a sublayer from the buffer main layers
        !
        !>@param update_bf_grdpts_after_increase
        !> compute the new grid points of the bf_sublayer after its
        !> increase and synchronize the neighboring buffer layers
        !
        !> @param determine_interior_bc_layers
        !> determine the boundary sections in the interior domain
        !
        !> @param determine_interior_bc_procedures
        !> determine the boundary sections and prcoedures for the
        !> interior domain
        !
        !> @param sync_nodes_with_interior
        !> synchronize the nodes at the edge between the interior
        !> domain and the buffer layers
        !
        !> @param sync_nodes_at_mainlayer_interfaces
        !> synchronize the nodes the edge between the buffer main
        !> layers
        !
        !> @param sync_nodes_at_domain_interfaces
        !> synchronize the nodes both at the edge between the
        !> interior and the buffer layers and at the edge between
        !> the buffer main layers
        !
        !>@param get_mainlayer_id
        !> get the cardinal coordinate corresponding to
        !> the general coordinates
        !
        !>@param get_sublayer
        !> get the sublayer corresponding to the general coordinates
        !> asked within the tolerance
        !
        !>@param get_nodes
        !> extract the governing variables at a general coordinates
        !> asked by the user
        !
        !>@param update_grdpts_from_neighbors
        !> if the buffer sublayer passed as argument has grid points
        !> in common with buffer layers from other main layers, the
        !> grid points in common are updated from the neighboring buffer
        !> layers
        !
        !>@param update_neighbor_grdpts
        !< if the buffer sublayer passed as argument has grid points
        !> in common with the buffer layers from other main layers, the
        !> grid points in common are updated in the neighboring buffer
        !> layers from the current buffer layer
        !
        !>@param shares_with_neighbor1
        !> check if the alignment of the future sublayer makes it a
        !> potential neighboring buffer layer of type 1 and if so
        !> update its alignment if it is shifted by one grid
        !> point from the border as it makes exchanges easier later the
        !> sublayer is here tested as a neighbor buffer layer of type 1
        !> this function echoes bf_layer_class%shares_grdpts_with_neighbor1
        !
        !>@param shares_with_neighbor2
        !> check if the alignment of the future sublayer makes it a
        !> potential neighboring buffer layer of type 2 and if so
        !> update its alignment if it is shifted by one grid
        !> point from the border as it makes exchanges easier later the
        !> sublayer is here tested as a neighbor buffer layer of type 1
        !> this function echoes bf_layer_class%shares_grdpts_with_neighbor1
        !
        !>@param get_nbf_layers_sharing_grdpts_with
        !> determine the neighboring buffer layers sharing grid points
        !> with the current buffer layer
        !
        !>@param bf_layer_depends_on_neighbors
        !> determine whether a buffer layer depends on its neighboring
        !> buffer layers
        !
        !>@param does_a_neighbor_remains
        !> check if the neighboring bf_layers from bf_sublayer_i can
        !> all be removed
        !
        !>@param resolve_bc_overlap_conflicts
        !> check if the buffer layers at the interface between main
        !> layers (N<->E, N<->W, S<->E, S<->W) have boundary layers
        !> that cannot be computed due to a lack of grid points (after
        !> the neighboring boundary layers updated their grid points w/o
        !> updating its allocation)
        !
        !>@param print_binary
        !> print the content of the interface on external binary files
        !
        !>@param print_netcdf
        !> print the content of the interface on external netcdf files
        !
        !> @param allocate_before_timeInt
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param deallocate_after_timeInt
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
        !---------------------------------------------------------------
        type :: bf_interface

          type(bf_mainlayer_pointer), dimension(4), private :: mainlayer_pointers
          type(nbf_interface_newgrdpt)            , private :: border_interface

          contains

          !buffer layer interactions management
          procedure, pass :: ini
          procedure, pass :: restart
          procedure, pass :: get_mainlayer

          procedure, pass :: allocate_sublayer
          procedure, pass :: reallocate_sublayer
          procedure, pass :: merge_sublayers
          procedure, pass :: remove_sublayer

          !procedure computation of new grid points
          procedure, pass :: update_bf_grdpts_after_increase

          !interior boundary procedures
          procedure, pass :: determine_interior_bc_layers
          procedure, pass :: determine_interior_bc_procedures

          !synchronization between the domains
          procedure, pass :: sync_nodes_with_interior
          procedure, pass :: sync_nodes_at_mainlayer_interfaces
          procedure, pass :: sync_nodes_at_domain_interfaces

          procedure, nopass :: get_mainlayer_id
          procedure, pass   :: get_sublayer
          procedure, pass   :: get_nodes

          procedure, pass, private :: update_grdpts_from_neighbors
          procedure, pass, private :: update_neighbor_grdpts
          
          procedure, nopass, private :: shares_with_neighbor1
          procedure, nopass, private :: shares_with_neighbor2

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

          procedure, pass :: resolve_bc_overlap_conflicts

          !i/o management
          procedure, pass :: print_binary
          procedure, pass :: print_netcdf

          !for time integration
          procedure, pass :: allocate_before_timeInt
          procedure, pass :: deallocate_after_timeInt
          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step
          procedure, pass :: update_integration_borders

        end type bf_interface

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer/interior domain interface
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        ! the data between them
        !--------------------------------------------------------------
        subroutine ini(this,interior_x_map,interior_y_map)

          implicit none

          class(bf_interface)       , intent(inout) :: this
          real(rkind), dimension(nx), intent(in)    :: interior_x_map
          real(rkind), dimension(ny), intent(in)    :: interior_y_map

          integer :: i
          
          do i=1, size(this%mainlayer_pointers,1)
             call this%mainlayer_pointers(i)%ini()             
          end do

          call this%border_interface%ini()

          !create a netcdf file with the grdpts_id for the interior
          !computational domain
          call print_interior_grdpts_id_on_netcdf(
     $         'interior_grdpts_id.nc',
     $         interior_x_map,
     $         interior_y_map)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer/interior domain interface
        !
        !> @date
        !> 15_12_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        ! the data between them
        !--------------------------------------------------------------
        subroutine restart(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     nb_bf_layers,
     $     p_model,
     $     timestep)

          implicit none

          class(bf_interface)                     , intent(inout) :: this
          real(rkind)        , dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind)        , dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind)        , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer            , dimension(4)       , intent(in)    :: nb_bf_layers
          type(pmodel_eq)                         , intent(in)    :: p_model
          integer                                 , intent(in)    :: timestep


          character(len=1) , dimension(4)                  :: bf_mainlayer_id
          integer                                          :: k
          integer                                          :: i
          integer                                          :: t_index_format
          integer                                          :: bf_index_format
          character(len=19)                                :: filename_format
          character(len=25)                                :: bf_filename
          integer                                          :: ncid
          integer          , dimension(3)                  :: coordinates_id
          integer          , dimension(ne)                 :: data_id
          integer                                          :: grdptsid_id
          real(rkind)      , dimension(2)                  :: x_borders
          real(rkind)      , dimension(2)                  :: y_borders
          integer          , dimension(2)                  :: sizes
          integer          , dimension(2,2)                :: bf_alignment
          type(bf_sublayer), pointer                       :: bf_sublayer_restart
          real(rkind)                                      :: time
          real(rkind)      , dimension(:)    , allocatable :: bf_x_map
          real(rkind)      , dimension(:)    , allocatable :: bf_y_map
          real(rkind)      , dimension(:,:,:), allocatable :: bf_nodes
          integer          , dimension(:,:)  , allocatable :: bf_grdpts_id

          bf_mainlayer_id(N) = 'N'
          bf_mainlayer_id(S) = 'S'
          bf_mainlayer_id(E) = 'E'
          bf_mainlayer_id(W) = 'W'


          !initialize the nbf_interface + pointers to main layers
          call this%ini(interior_x_map,interior_y_map)

          !initialize the buffer layers
          !1) extract the number of buffer layers for each localization
          !2) extract the buffer layer interior_map and the nodes
          do k=1,4
             do i=1, nb_bf_layers(k)

                !number of character to represent the timestep index
                if(timestep.eq.0) then
                   t_index_format = 1
                else
                   t_index_format  = floor(log10(real(timestep))+1.0)
                end if

                !number of character to represent the bf_layer index
                if(bf_index_format.eq.0) then
                   bf_index_format = 1
                else
                   bf_index_format = floor(log10(real(i))+1.0)
                end if

                !format for the filename
                write(filename_format,
     $               '(''(A1,A1,I'',I1,'',A1,I'',I1,'',A3)'')')
     $               bf_index_format,
     $               t_index_format

                !filename
                write(bf_filename, filename_format)
     $               bf_mainlayer_id(k),'_',
     $               i,'_',
     $               timestep,
     $               '.nc'

                !extract [x_min,x_max]x[y_min,y_max]
                call nf90_open_file_for_reading(
     $               trim(bf_filename),
     $               ncid)

                call nf90_get_varid(
     $               ncid,
     $               p_model,
     $               coordinates_id,
     $               data_id,
     $               grdptsid_id=grdptsid_id)

                !read borders
                call nf90_read_borders(
     $               ncid,
     $               coordinates_id, 
     $               x_borders,
     $               y_borders,
     $               sizes)

                !determine the alignment of the buffer layer
                bf_alignment = get_restart_alignment(
     $               interior_x_map,
     $               interior_y_map,
     $               x_borders,
     $               y_borders)

                !allocate the new buffer layer
                bf_sublayer_restart => allocate_sublayer(
     $               this,
     $               k,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               bf_alignment)

                !determine the x_map,y_map and nodes for
                !the buffer layer
                allocate(bf_x_map(sizes(1)))
                allocate(bf_y_map(sizes(2)))
                allocate(bf_nodes(sizes(1),sizes(2),ne))
                allocate(bf_grdpts_id(sizes(1),sizes(2)))

                call nf90_get_var_model_nopt(
     $               ncid,
     $               coordinates_id,
     $               data_id,
     $               time, bf_nodes, bf_x_map, bf_y_map,
     $               [1,sizes(1),sizes(2)],
     $               grdptsid_id=grdptsid_id,grdpts_id=bf_grdpts_id)

                call bf_sublayer_restart%set_x_map(bf_x_map)
                call bf_sublayer_restart%set_y_map(bf_y_map)
                call bf_sublayer_restart%set_nodes(bf_nodes)
                call bf_sublayer_restart%set_grdpts_id(bf_grdpts_id)

                !close file
                call nf90_close_file(ncid)

                !update the integration borders
                call this%update_integration_borders(bf_sublayer_restart)

                !update the neighboring buffer layers
                call this%update_neighbor_grdpts(bf_sublayer_restart)

             end do
          end do
          

        end subroutine restart


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the bf_mainlayer object corresponding to the cardinal point
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param mainlayer_id
        !> cardinal coordinate for the bf_mainlayer asked
        !
        !>@return get_mainlayer
        !> reference to the bf_mainlayer corresponding to the cardinal
        !> coordinate
        !--------------------------------------------------------------
        function get_mainlayer(this, mainlayer_id)
       
         implicit none
   
         class(bf_interface), intent(in) :: this
         integer            , intent(in) :: mainlayer_id
         type(bf_mainlayer) , pointer    :: get_mainlayer

         if(debug) then
            if((mainlayer_id.lt.1).or.(mainlayer_id.gt.4)) then
               call error_mainlayer_id(
     $              'bf_interface_class',
     $              'get_mainlayer',
     $              mainlayer_id)
            end if
         end if

         if(this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
         else
            nullify(get_mainlayer)
         end if

       end function get_mainlayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if the alignment of the future sublayer makes it a
       !> potential neighboring buffer layer of type 1 and if so
       !> update its alignment if it is shifted by one grid
       !> point from the border as it makes exchanges easier later the
       !> sublayer is here tested as a neighbor buffer layer of type 1
       !> this function echoes bf_layer_class%shares_grdpts_with_neighbor1
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param mainlayer_id
       !> cardinal coordinate for the bf_mainlayer asked
       !
       !>@param bf_final_alignment
       !> position of the future buffer layer
       !
       !>@return share_grdpts
       !> logical stating whether the future sublayer is a potential
       !> buffer layer at the edge between main layers
       !--------------------------------------------------------------
       function shares_with_neighbor1(mainlayer_id, bf_final_alignment)
     $     result(share_grdpts)

         implicit none

         integer                       , intent(in)    :: mainlayer_id
         integer(ikind), dimension(2,2), intent(inout) :: bf_final_alignment
         logical                                       :: share_grdpts


         select case(mainlayer_id)
           case(N,S)
              share_grdpts = bf_final_alignment(1,1).le.(align_W+bc_size+1)
           case(E,W)
              share_grdpts = bf_final_alignment(2,1).le.(align_S+bc_size+1)
              if(bf_final_alignment(2,1).eq.(align_S+bc_size)) then
                 bf_final_alignment(2,1)=align_S+1
              end if
           case default
              call error_mainlayer_id(
     $             'bf_layer_class.f',
     $             'share_grdpts_with_neighbor1',
     $             mainlayer_id)
         end select

       end function shares_with_neighbor1


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if the alignment of the future sublayer makes it a
       !> potential neighboring buffer layer of type 2 and if so
       !> update its alignment if it is shifted by one grid
       !> point from the border as it makes exchanges easier later the
       !> sublayer is here tested as a neighbor buffer layer of type 1
       !> this function echoes bf_layer_class%shares_grdpts_with_neighbor1
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param mainlayer_id
       !> cardinal coordinate for the bf_mainlayer asked
       !
       !>@param bf_final_alignment
       !> position of the future buffer layer
       !
       !>@return share_grdpts
       !> logical stating whether the future sublayer is a potential
       !> buffer layer at the edge between main layers
       !--------------------------------------------------------------
       function shares_with_neighbor2(mainlayer_id, bf_final_alignment)
     $     result(share_grdpts)

         implicit none

         integer                       , intent(in)    :: mainlayer_id
         integer(ikind), dimension(2,2), intent(inout) :: bf_final_alignment
         logical                                       :: share_grdpts


         select case(mainlayer_id)
           case(N,S)
              share_grdpts = bf_final_alignment(1,2).ge.(align_E-bc_size-1)
           case(E,W)
              share_grdpts = bf_final_alignment(2,2).ge.(align_N-bc_size-1)
              if(bf_final_alignment(2,2).eq.(align_N-bc_size)) then
                 bf_final_alignment(2,2)=align_N-1
              end if
           case default
              call error_mainlayer_id(
     $             'bf_layer_class.f',
     $             'share_grdpts_with_neighbor1',
     $             mainlayer_id)
         end select

       end function shares_with_neighbor2


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> allocate a bf_sublayer and insert it in the corresponding
       !> buffer mainlayer
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param mainlayer_id
       !> cardinal coordinate for the new bf_sublayer
       !
       !>@param nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@param alignment
       !> table of integers characterizing the
       !> correspondance between the interior grid points
       !> and the buffer layer elements
       !
       !>@return added_sublayer
       !> reference to the newly allocated bf_sublayer
       !--------------------------------------------------------------
       function allocate_sublayer(
     $     this,
     $     mainlayer_id,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(added_sublayer)
        
          class(bf_interface)             , intent(inout) :: this
          integer                         , intent(in)    :: mainlayer_id
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer, dimension(2,2)         , intent(inout) :: alignment

          type(bf_sublayer), pointer                      :: added_sublayer


          logical :: share_with_neighbor1
          logical :: share_with_neighbor2


          !0) debug : check the mainlayer id
          if(debug) then
             if((mainlayer_id.lt.1).or.(mainlayer_id.gt.4)) then
                call error_mainlayer_id(
     $              'bf_interface_class',
     $              'allocate_sublayer',
     $              mainlayer_id)
             end if
          end if

          !1) check if the new sublayer is at the interface between several
          !   main layers and update the alignment in case it is just one grid
          !   point away from a corner: exchanges are made easier
          share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
          share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)


          !2) check if the mainlayer corresponding to the cardinal point is
          !   indeed allocated: if the memory space is not allocated, the space
          !   in memory is first allocated, the pointer identifying the mainlayer
          !   is initialized and the main layer itself is initialized
          if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
             call this%mainlayer_pointers(mainlayer_id)%ini_mainlayer(mainlayer_id)
          end if


          !3) now that we are sure that space is allocated for the main layer,
          !   the sublayer can be integrated to the mainlayer and the buffer
          !   layer can be initialized using the nodes, alignment and neighbors
          !   arguments
          added_sublayer => this%mainlayer_pointers(mainlayer_id)%add_sublayer(
     $         interior_x_map, interior_y_map, interior_nodes, alignment)
          call added_sublayer%set_neighbor1_share(share_with_neighbor1)
          call added_sublayer%set_neighbor2_share(share_with_neighbor2)

          !4) if the sublayer newly allocated is indeed a buffer layer at the
          !   interface between main layers and of type 1, the neighboring
          !   buffer layers should know that they will exchange data with this
          !   buffer layer.          
          !   The same is true for an interface buffer layer of type 2          
          if(share_with_neighbor1) then
             call this%border_interface%link_neighbor1_to_bf_sublayer(
     $            added_sublayer)
          end if
          if(share_with_neighbor2) then
             call this%border_interface%link_neighbor2_to_bf_sublayer(
     $            added_sublayer)
          end if
          
          !5) Morover, the new buffer layer has been allocated without considering
          !   the other buffer layers already allocated which could share grid
          !   points with this buffer layer. These grid points are now copy from
          !   the other buffer layers to this buffer layer
          call this%border_interface%update_grdpts_from_neighbors(added_sublayer)


          !6) The integration borders are computed for the new buffer layer
          !   and the neighboring buffer layers are asked to update their own
          !   integration borders
          call this%update_integration_borders(added_sublayer)

       end function allocate_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> reallocate a buffer sublayer and check whether the neighboring
       !> buffer layer changed and so if the neighboring links should be
       !> updated
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_r
       !> reference to the bf_sublayer to be reallocated
       !
       !>@param interior_x_map
       !> array encapsulating the coordinates along the x-axis for
       !> the interior computational domain
       !
       !>@param interior_y_map
       !> array encapsulating the coordinates along the y-axis for
       !> the interior computational domain
       !
       !>@param interior_nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@param alignment
       !> table of integers characterizing the
       !> correspondance between the interior grid points
       !> and the buffer layer
       !--------------------------------------------------------------
       subroutine reallocate_sublayer(
     $     this,
     $     bf_sublayer_r,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         type(bf_sublayer), pointer      , intent(inout) :: bf_sublayer_r
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
         integer    , dimension(2,2)     , intent(inout) :: alignment

         integer :: mainlayer_id
         logical :: share_with_neighbor1
         logical :: share_with_neighbor2

         mainlayer_id = bf_sublayer_r%get_localization()


         !1) check if the reallocated sublayer is at the interface between
         !   several main layers and update the alignment in case it is just
         !   one grid point away from a corner: exchanges are made easier
         share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
         share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)

         
         !2) reallocate the buffer sublayer
         call bf_sublayer_r%reallocate_bf_layer(
     $        interior_x_map,
     $        interior_y_map,
     $        interior_nodes,
     $        alignment)

         !3) check if the links to the neighbor1 should be updated
         ! if the bf_sublayer_r was exchanging with neighbor1 before
         if(bf_sublayer_r%can_exchange_with_neighbor1()) then
            
            !and the bf_sublayer_r can no longer exchange
            !the previous links should be removed
            if(.not.share_with_neighbor1) then
               call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         ! if the bf_sublayer_r was not exchanging with neighbor1 before
         else

            !and the bf_sublayer_r can now exchange, links should be added
            !to this bf_sublayer_r
            if(share_with_neighbor1) then
               call this%border_interface%link_neighbor1_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         end if

         !4) check if the links to the neighbor2 should be updated
         ! if the bf_sublayer_r was exchanging with neighbor2 before
         if(bf_sublayer_r%can_exchange_with_neighbor2()) then
            
            !and the bf_sublayer_r can no longer exchange
            !the previous links should be removed
            if(.not.share_with_neighbor2) then
               call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         ! if the bf_sublayer_r was not exchanging with neighbor2 before
         else

            !and the bf_sublayer_r can now exchange, links should be added
            !to this bf_sublayer_r
            if(share_with_neighbor2) then
               call this%border_interface%link_neighbor2_to_bf_sublayer(
     $              bf_sublayer_r)
            end if

         end if

         !5) update the status of the bf_sublayer_r for its neighbors
         call bf_sublayer_r%set_neighbor1_share(share_with_neighbor1)
         call bf_sublayer_r%set_neighbor2_share(share_with_neighbor2)

         !6) update the grid points new allocated using the neighboring sublayers
         call this%border_interface%update_grdpts_from_neighbors(bf_sublayer_r)
         
         !7) The integration borders are computed for the reallocated buffer
         !   layer and the neighboring buffer layers are asked to update their
         !   own integration borders
         call this%update_integration_borders(bf_sublayer_r)

       end subroutine reallocate_sublayer
       

       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> merge the content of two sublayers and update
       !> the links to the border buffer layers
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer1
       !> reference to the first bf_sublayer merged
       !
       !>@param bf_sublayer2
       !> reference to the second bf_sublayer merged
       !
       !>@param nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@param alignment
       !> table of integers characterizing the
       !> correspondance between the interior grid points
       !> and the buffer layer
       !--------------------------------------------------------------
       function merge_sublayers(
     $     this,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     alignment)
     $     result(merged_sublayer)

          implicit none

          class(bf_interface)                        , intent(inout) :: this
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer1
          type(bf_sublayer), pointer                 , intent(inout) :: bf_sublayer2
          real(rkind)      , dimension(nx)           , intent(in)    :: interior_x_map
          real(rkind)      , dimension(ny)           , intent(in)    :: interior_y_map
          real(rkind)      , dimension(nx,ny,ne)     , intent(in)    :: interior_nodes
          integer(ikind)   , dimension(2,2), optional, intent(inout) :: alignment
          type(bf_sublayer), pointer                                 :: merged_sublayer

          integer :: mainlayer_id
          logical :: share_with_neighbor1
          logical :: share_with_neighbor2


          mainlayer_id = bf_sublayer1%get_localization()
   

          if(present(alignment)) then


             !1) check if the sublayer resulting from the merge is at the interface
             !   between several main layers and update the alignment in case it is just
             !   one grid point away from a corner: exchanges are made easier
             share_with_neighbor1 = shares_with_neighbor1(mainlayer_id, alignment)
             share_with_neighbor2 = shares_with_neighbor2(mainlayer_id, alignment)


             !2) check if the links to neighbor1 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor1()) then
                
                !if the bf_sublayer2 was also linked to neighbor1 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor1()) then
                   call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor1 sublayers
             else

                !and the bf_sublayer2 was linked to this neighbor, the links should
                !be updated from bf_sublayer2 to bf_sublayer1 and the status of the
                !neighbor1 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor1()) then

                   call this%border_interface%update_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor1_share(.true.)


                !and if the bf_sublayer2 was also not linked to neighbor1
                else
                   
                   !and if the final alignment is such that the merged sublayer
                   !will exchange with neighbor1, links should be created to
                   !bf_sublayer1 and the state of the neighbor1 for bf_sublayer1
                   !should be updated
                   if(share_with_neighbor1) then
                      
                      call this%border_interface%link_neighbor1_to_bf_sublayer(
     $                     bf_sublayer1)

                      call bf_sublayer1%set_neighbor1_share(.true.)

                   end if
                end if
             end if


             !3) check if the links to neighbor2 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor2()) then
                
                !if the bf_sublayer2 was also linked to neighbor2 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor2()) then
                   call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor2 sublayers
             else

                !and on the contrary the bf_sublayer2 was linked to this neighbor
                !the links should be updated from bf_sublayer2 to bf_sublayer1
                !and the status of the neighbor2 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor2()) then

                   call this%border_interface%update_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor2_share(.true.)


                !and if the bf_sublayer2 was also not linked to neighbor2
                else
                   
                   !and if the final alignment is such that the merged sublayer
                   !will exchange with neighbor2, links should be created to
                   !bf_sublayer1 and the state of the neighbor2 for bf_sublayer1
                   !should be updated
                   if(share_with_neighbor2) then
                      
                      call this%border_interface%link_neighbor2_to_bf_sublayer(
     $                     bf_sublayer1)

                      call bf_sublayer1%set_neighbor2_share(.true.)

                   end if
                end if
             end if


             !4) merge the content of the sublayers
             merged_sublayer => this%mainlayer_pointers(mainlayer_id)%merge_sublayers(
     $            bf_sublayer1,
     $            bf_sublayer2,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            alignment)

          else

             !2) check if the links to neighbor1 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor1()) then
                
                !if the bf_sublayer2 was also linked to neighbor1 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor1()) then
                   call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor1 sublayers
             else

                !and the bf_sublayer2 was linked to this neighbor, the links should
                !be updated from bf_sublayer2 to bf_sublayer1 and the status of the
                !neighbor1 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor1()) then

                   call this%border_interface%update_link_from_neighbor1_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor1_share(.true.)

                end if
             end if


             !3) check if the links to neighbor2 should be updated
             !if the bf_sublayer1 was already linked to neighbor1 sublayers
             if(bf_sublayer1%can_exchange_with_neighbor2()) then
                
                !if the bf_sublayer2 was also linked to neighbor2 sublayers
                !as it will be merged, these links should be removed
                if(bf_sublayer2%can_exchange_with_neighbor2()) then
                   call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer2)

                end if

             !if the bf_sublayer1 was not linked to neighbor2 sublayers
             else

                !and on the contrary the bf_sublayer2 was linked to this neighbor
                !the links should be updated from bf_sublayer2 to bf_sublayer1
                !and the status of the neighbor2 for bf_sublayer1 should be updated
                if(bf_sublayer2%can_exchange_with_neighbor2()) then

                   call this%border_interface%update_link_from_neighbor2_to_bf_sublayer(
     $                  bf_sublayer1, bf_sublayer2)

                   call bf_sublayer1%set_neighbor2_share(.true.)

                end if
             end if


             !4) merge the content of the sublayers
             merged_sublayer => this%mainlayer_pointers(mainlayer_id)%merge_sublayers(
     $            bf_sublayer1,
     $            bf_sublayer2,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes)

          end if


          !5) update the grid points that have been newly allocated
          !   with grid points from the neighboring sublayers
          call this%border_interface%update_grdpts_from_neighbors(bf_sublayer1)


c$$$          !6) if the merge is such that grid points between the two
c$$$          !   sublayers should be updated to prevent a line of inconsistent
c$$$          !   boundary points
c$$$          !    ________  ________        __________________
c$$$          !   |        ||        |      |                  |
c$$$          !   |        ||        |  ->  |                  |
c$$$          !   |        ||        |      |                  |
c$$$          print '(''bf_interface_class'')'
c$$$          print '(''merge_sublayers'')'
c$$$          stop 'not implemented yet'

         !7) The integration borders are computed for the buffer layer resulting
         !   from the merge and the neighboring buffer layers are asked to update
         !   their own integration borders
         call this%update_integration_borders(bf_sublayer1)

       end function merge_sublayers


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> remove a sublayer from the buffer main layers
       !
       !> @date
       !> 27_06_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param sublayer_ptr
       !> reference to the bf_sublayer removed
       !
       !>@param bf_mainlayer_id
       !> cardinal coordinate of the buffer layer removed
       !--------------------------------------------------------------
       subroutine remove_sublayer(
     $     this,
     $     sublayer_ptr,
     $     interior_x_map,
     $     interior_y_map,
     $     bf_mainlayer_id)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         type(bf_sublayer), pointer      , intent(inout) :: sublayer_ptr
         real(rkind)      , dimension(nx), intent(in)    :: interior_x_map
         real(rkind)      , dimension(ny), intent(in)    :: interior_y_map
         integer, optional               , intent(in)    :: bf_mainlayer_id

         
         integer     :: mainlayer_id
         real(rkind) :: dx_s
         real(rkind) :: dy_s

         dx_s = interior_x_map(2)-interior_x_map(1)
         dy_s = interior_y_map(2)-interior_y_map(1)


         !> identify the mainlayer to which the sublayer belongs
         if(present(bf_mainlayer_id)) then
            mainlayer_id = bf_mainlayer_id
         else
            mainlayer_id = sublayer_ptr%get_localization()
         end if

         !> remove the sublayer from the table identifying the
         !> neighboring buffer layers
         if(sublayer_ptr%can_exchange_with_neighbor1()) then
            call this%border_interface%remove_link_from_neighbor1_to_bf_sublayer(
     $           sublayer_ptr)
         end if
         if(sublayer_ptr%can_exchange_with_neighbor2()) then
            call this%border_interface%remove_link_from_neighbor2_to_bf_sublayer(
     $           sublayer_ptr)
         end if

         !> remove the sublayer from the main layer
         call this%mainlayer_pointers(mainlayer_id)%remove_sublayer(sublayer_ptr)

       end subroutine remove_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> get the cardinal coordinate corresponding to
       !> the general coordinates
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param general_coord
       !> integer table giving the general coordinates
       !
       !>@return mainlayer_id
       !> main layer cardinal coordinates
       !--------------------------------------------------------------
       function get_mainlayer_id(general_coord) result(mainlayer_id)

         implicit none

         integer(ikind), dimension(2), intent(in) :: general_coord
         integer                                  :: mainlayer_id

         if(general_coord(2).le.align_S) then
            mainlayer_id = S

         else
            if(general_coord(2).lt.(align_N)) then

               if(general_coord(1).le.align_W) then
                  mainlayer_id = W

               else
                  if(general_coord(1).lt.align_E) then
                     mainlayer_id = interior

                  else
                     mainlayer_id = E

                  end if
               end if
                     
            else
               mainlayer_id = N

            end if
         end if

       end function get_mainlayer_id


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> get the sublayer corresponding to the general coordinates
       !> asked within the tolerance
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param general_coord
       !> table giving the general coordinates of the point analyzed
       !
       !>@param local_coord
       !> table giving the local coordinates of the point analyzed
       !> in the corresponding bf_sublayer
       !
       !>@param tolerance_i
       !> integer indicating how far the gridpoint can be from the
       !> closest sublayer to be considered inside
       !
       !>@param mainlayer_id_i
       !> cardinal coordinate cooresponding to the general coordinates
       !
       !>@return sublayer
       !> reference to the bf_sublayer matching the general coordinates
       !> of the grid point
       !--------------------------------------------------------------
       function get_sublayer(
     $    this,
     $    general_coord,
     $    local_coord,
     $    tolerance_i,
     $    mainlayer_id_i)
     $    result(sublayer)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         integer(ikind), dimension(2), intent(in)  :: general_coord
         integer(ikind), dimension(2), intent(out) :: local_coord
         integer       , optional    , intent(in)  :: tolerance_i
         integer       , optional    , intent(in)  :: mainlayer_id_i
         type(bf_sublayer), pointer                :: sublayer


         integer                     :: direction_tested
         integer                     :: mainlayer_id
         type(bf_mainlayer), pointer :: mainlayer
         integer                     :: tolerance
         logical                     :: grdpt_in_sublayer


         !< identification of the main layer
         if(present(mainlayer_id_i)) then
            mainlayer_id = mainlayer_id_i
         else
            mainlayer_id = get_mainlayer_id(general_coord)
         end if

         !< if the general coordinates match the interior,
         !> no sublayer matches the general coordinates
         !> and the sublayer pointer is nullified
         if(mainlayer_id.eq.interior) then
            nullify(sublayer)

         !< otherwise, the mainlayers are analyzed
         else

            !< check that the main layer exists
            !< if it does not exist, no sublayer can be saved inside
            !< and the pointer to the sublayer is nullified
            if(.not.this%mainlayer_pointers(mainlayer_id)%associated_ptr()) then
               nullify(sublayer)
                 
            !< if the main layer exists, the sublayers saved inside are
            !< checked to decide whether the grid point asked belongs to
            !< one of them or not
            else
               mainlayer => this%mainlayer_pointers(mainlayer_id)%get_ptr()
            
               !< check if sublayers are saved inside the mainlayer
               !> if no sublayers are saved inside the mainlayer,
               !> no existing sublayer can match the general coord
               !> and so the pointer to sublayer is nullified
               if(.not.associated(mainlayer%get_head_sublayer())) then
                  nullify(sublayer)
            
               !< otherwise, the sublayer corresponding to the general
               !> coordinates is searched by going through the different
               !> element of the doubled chained list
               else
                  sublayer => mainlayer%get_head_sublayer()
            
              	  !< processing the tolerance for matching a sublayer
            	  !> if no tolerence is provided, the default option is 0
                  if(.not.present(tolerance_i)) then
                     tolerance=0
                  else
                     tolerance=tolerance_i
                  end if
                  
                  !< if the mainlayer investigated is N,S,E or W, there can
                  !> be sublayers to be investigated
                  select case(mainlayer_id)
                    case(N,S)
                       direction_tested = x_direction
                    case(E,W)
                       direction_tested = y_direction
                  end select 

                  !check if the grid point belongs to the current sublayer
                  grdpt_in_sublayer =
     $                 (  general_coord(direction_tested).ge.
     $                 (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                 .and.(
     $                 general_coord(direction_tested).le.
     $                 (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  !go through the different sublayers
                  do while(.not.grdpt_in_sublayer)
            	        
            	     !if no matching sublayer can be found
            	     !nullify the corresponding pointer
                     if(.not.associated(sublayer%get_next())) then
                        nullify(sublayer)
                        exit
                     end if
                   
                     sublayer => sublayer%get_next()
                     grdpt_in_sublayer =
     $                    (general_coord(direction_tested).ge.
     $                    (sublayer%get_alignment(direction_tested,1)-bc_size-tolerance))
     $                    .and.(
     $                    general_coord(1).le.
     $                    (sublayer%get_alignment(direction_tested,2)+bc_size+tolerance))
            	
                  end do
            
               end if
            end if
         end if

         !< if a sublayer matching the general coordinates was found
         !> compute the local coordinates in this sublayer
         if(associated(sublayer)) then
            local_coord = sublayer%get_local_coord(general_coord)
         end if           
         
       end function get_sublayer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> extract the governing variables at a general coordinates
       !> asked by the user
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param g_coord
       !> table giving the general coordinates of the point analyzed
       !
       !>@param interior_nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !
       !>@return var
       !> governing variables at the grid point asked
       !--------------------------------------------------------------
       function get_nodes(this, g_coords, interior_nodes) result(var)

         implicit none
         
         class(bf_interface)             , intent(in) :: this
         integer(ikind), dimension(2)    , intent(in) :: g_coords
         real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
         real(rkind), dimension(ne)                   :: var

         integer                      :: mainlayer_id
         type(bf_sublayer), pointer   :: sublayer
         integer(ikind), dimension(2) :: l_coords
         
         mainlayer_id = this%get_mainlayer_id(g_coords)
         if((g_coords(1).ge.1).and.
     $        (g_coords(1).le.nx).and.
     $        (g_coords(2).ge.1).and.
     $        (g_coords(2).le.ny)) then
            var = interior_nodes(g_coords(1),g_coords(2),:)
         else
            sublayer => this%get_sublayer(
     $               g_coords, l_coords, mainlayer_id_i=mainlayer_id)
            if(associated(sublayer)) then
               var = sublayer%get_nodes(l_coords)
            else
               print '(''bf_interface_class'')'
               print '(''get_nodes'')'
               print '(''cannot get sublayer'')'
               stop 'check way to get nodes'
            end if               
         end if

       end function get_nodes    


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> if the buffer sublayer passed as argument has grid points
       !> in common with buffer layers from other main layers, the
       !> grid points in common are updated from the neighboring buffer
       !> layers
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param nbf_sublayer
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !--------------------------------------------------------------
       subroutine update_grdpts_from_neighbors(this, nbf_sublayer)

         implicit none

         class(bf_interface), intent(in)    :: this
         type(bf_sublayer)  , intent(inout) :: nbf_sublayer

         call this%border_interface%update_grdpts_from_neighbors(
     $        nbf_sublayer)

       end subroutine update_grdpts_from_neighbors


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !< if the buffer sublayer passed as argument has grid points
       !> in common with the buffer layers from other main layers, the
       !> grid points in common are updated in the neighboring buffer
       !> layers from the current buffer layer
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param nbf_sublayer
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !--------------------------------------------------------------
       subroutine update_neighbor_grdpts(this, nbf_sublayer)

         implicit none

         class(bf_interface), intent(inout) :: this
         type(bf_sublayer)  , intent(inout) :: nbf_sublayer

         call this%border_interface%update_neighbor_grdpts(nbf_sublayer)

       end subroutine update_neighbor_grdpts


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> compute the new grid points of the bf_sublayer after its
       !> increase and synchronize the neighboring buffer layers
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_i
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !
       !>@param selected_grdpts
       !> list of the general coordinates of the grid points to be
       !> computed
       !--------------------------------------------------------------
       subroutine update_bf_grdpts_after_increase(
     $     this,
     $     bf_sublayer_updated,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     selected_grdpts)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         type(bf_sublayer)               , intent(inout) :: bf_sublayer_updated
         type(pmodel_eq)                 , intent(in)    :: p_model
         real(rkind)                     , intent(in)    :: t
         real(rkind)                     , intent(in)    :: dt
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1
         integer(ikind), dimension(:,:)  , intent(in)    :: selected_grdpts

         !compute the new grid points after the increase
         call this%border_interface%update_bf_grdpts_after_increase(
     $        bf_sublayer_updated,
     $        p_model,
     $        t,dt,
     $        interior_x_map,
     $        interior_y_map,
     $        interior_nodes0,
     $        interior_nodes1,
     $        selected_grdpts)

         !update the neighboring buffer layers
         call this%update_neighbor_grdpts(bf_sublayer_updated)

       end subroutine update_bf_grdpts_after_increase


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> determine the extent of the interior boundary
       !> layers
       !
       !> @date
       !> 29_10_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param interior_bc_sections_N
       !> extent of the boundary sections for the north interior
       !> boundary layer
       !
       !>@param interior_bc_sections_S
       !> extent of the boundary sections for the south interior
       !> boundary layer
       !
       !>@param interior_bc_sections_E
       !> extent of the boundary sections for the east interior
       !> boundary layer
       !
       !>@param interior_bc_sections_W
       !> extent of the boundary sections for the west interior
       !> boundary layer
       !--------------------------------------------------------------
       subroutine determine_interior_bc_layers(
     $     this,
     $     interior_bc_sections_N,
     $     interior_bc_sections_S,
     $     interior_bc_sections_E,
     $     interior_bc_sections_W)

         implicit none

         class(bf_interface)                        , intent(in)    :: this
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_N
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_S
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_E
         integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_W

          type(bf_mainlayer) , pointer :: get_mainlayer


         !North boundary layer
         !---------------------
         !if there are buffer sub-layers saved in the buffer main layer
         !the extents of the interior buffer layers are computed using
         !the buffer main layer
         if(this%mainlayer_pointers(N)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(N)%get_ptr()

            call get_mainlayer%determine_interior_bc_layers(
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
         if(this%mainlayer_pointers(S)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(S)%get_ptr()

            call get_mainlayer%determine_interior_bc_layers(
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
         if(this%mainlayer_pointers(E)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(E)%get_ptr()

            call get_mainlayer%determine_interior_bc_layers(
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
         if(this%mainlayer_pointers(W)%associated_ptr()) then
            get_mainlayer => this%mainlayer_pointers(W)%get_ptr()

            call get_mainlayer%determine_interior_bc_layers(
     $           interior_bc_sections_W)

         !otherwise, there are no buffer layers and the extent of the
         !west interior boundary layer computed using the boundary
         !conditions is the entire west boundary layer
         else

            allocate(interior_bc_sections_W(2,1))
            interior_bc_sections_W(:,1) = [1+bc_size,ny-bc_size]

         end if

       end subroutine determine_interior_bc_layers


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> determine the extent of the interior boundary
       !> layers
       !
       !> @date
       !> 29_10_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bc_procedures
       !> array identifying the boundary procedures to be applied in
       !> the boundary layers of the interior domain
       !--------------------------------------------------------------
       subroutine determine_interior_bc_procedures(
     $     this,
     $     bc_procedures)

         implicit none

         class(bf_interface)                        , intent(inout) :: this
         integer(ikind), dimension(:,:), allocatable, intent(out)   :: bc_procedures

         
         integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_N
         integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_S
         integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_E
         integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_W


         call determine_interior_bc_layers(
     $        this,
     $        interior_bc_sections_N,
     $        interior_bc_sections_S,
     $        interior_bc_sections_E,
     $        interior_bc_sections_W)

         call process_bc_sections_into_bc_procedure(
     $        interior_bc_sections_N,
     $        interior_bc_sections_S,
     $        interior_bc_sections_E,
     $        interior_bc_sections_W,
     $        bc_procedures)

         deallocate(interior_bc_sections_N)
         deallocate(interior_bc_sections_S)
         deallocate(interior_bc_sections_E)
         deallocate(interior_bc_sections_W)

       end subroutine determine_interior_bc_procedures


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> exchange the grid points common between the buffer
       !> layers and the interior domain
       !
       !> @date
       !> 29_10_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param interior_nodes
       !> grid points from the interior domain
       !--------------------------------------------------------------
       subroutine sync_nodes_with_interior(
     $     this,
     $     interior_nodes)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

         integer :: i
         integer, dimension(4) :: exch_order

         exch_order = [E,W,N,S]

         do i=1, size(exch_order,1)

            if(this%mainlayer_pointers(exch_order(i))%associated_ptr()) then
               
               call this%mainlayer_pointers(exch_order(i))%sync_nodes_with_interior(
     $              interior_nodes)

            end if
         end do

       end subroutine sync_nodes_with_interior


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> synchronize the nodes located at the interface between
       !> buffer main layers
       !
       !> @date
       !> 30_10_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !--------------------------------------------------------------
       subroutine sync_nodes_at_mainlayer_interfaces(this)

         implicit none

         class(bf_interface), intent(inout) :: this
                  
         call this%border_interface%sync_interface_nodes()

       end subroutine sync_nodes_at_mainlayer_interfaces


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> synchronize the nodes located at the interface between
       !> interior and the buffer layers and the nodes located
       !> at the interface between the buffer main layers
       !
       !> @date
       !> 30_10_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param interior_nodes
       !> grid points from the interior domain
       !--------------------------------------------------------------
       subroutine sync_nodes_at_domain_interfaces(this, interior_nodes)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes
                  
         call this%sync_nodes_with_interior(interior_nodes)
         call this%sync_nodes_at_mainlayer_interfaces()

       end subroutine sync_nodes_at_domain_interfaces


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> determine the neighboring buffer layers sharing grid points
       !> with the current buffer layer
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_i
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !
       !>@param nbf1_list
       !> list of the neighboring buffer layers of type 1 sharing grid
       !> points with bf_sublayer_i
       !
       !>@param nbf2_list
       !> list of the neighboring buffer layers of type 2 sharing grid
       !> points with bf_sublayer_i
       !
       !>@param bf_mainlayer_id
       !> cardinal coordinate of bf_sublayer_i
       !--------------------------------------------------------------
       subroutine get_nbf_layers_sharing_grdpts_with(
     $     this, bf_sublayer_i, nbf1_list, nbf2_list, bf_mainlayer_id)

         implicit none
         
         class(bf_interface)       , intent(in)    :: this
         type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
         type(sbf_list)            , intent(inout) :: nbf1_list
         type(sbf_list)            , intent(inout) :: nbf2_list
         integer         , optional, intent(in)    :: bf_mainlayer_id


         !determine inside the list of neighboring buffer layer of
         !type 1 which ones share grid points with the current
         !buffer layer
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            if(present(bf_mainlayer_id)) then
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              1, bf_sublayer_i, nbf1_list, bf_mainlayer_id)
            else
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              1, bf_sublayer_i, nbf1_list)
            end if
         end if


         !determine inside the list of neighboring buffer layer of
         !type 2 which ones share grid points with the current
         !buffer layer
         if(bf_sublayer_i%can_exchange_with_neighbor2()) then
            if(present(bf_mainlayer_id)) then
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              2, bf_sublayer_i, nbf2_list, bf_mainlayer_id)
            else
               call this%border_interface%get_nbf_layers_sharing_grdpts_with(
     $              2, bf_sublayer_i, nbf2_list)
            end if
         end if         

       end subroutine get_nbf_layers_sharing_grdpts_with


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> determine whether a buffer layer depends on its neighboring
       !> buffer layers
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_i
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !
       !>@param bf_mainlayer_id
       !> cardinal coordinate of bf_sublayer_i
       !
       !>@return dependent
       !> logical stating whether the bf_sublayer_i depends on neighboring
       !> buffer layers
       !--------------------------------------------------------------
       function bf_layer_depends_on_neighbors(
     $     this, bf_sublayer_i, bf_mainlayer_id)
     $     result(dependent)

         implicit none

         class(bf_interface)       , intent(in) :: this
         type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
         integer         , optional, intent(in) :: bf_mainlayer_id
         logical                                :: dependent
       
         
         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 1
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            
            if(present(bf_mainlayer_id)) then
               dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $              1, bf_sublayer_i, bf_mainlayer_id)
            else
               dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $              1, bf_sublayer_i)
            end if
         else
            dependent = .false.
         end if
         

         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 2
         if(.not.dependent) then
            
            if(bf_sublayer_i%can_exchange_with_neighbor2()) then
               
               if(present(bf_mainlayer_id)) then
                  dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $                 2, bf_sublayer_i, bf_mainlayer_id)
               else
                  dependent = this%border_interface%bf_layer_depends_on_neighbors(
     $                 2, bf_sublayer_i)
               end if

            else
               dependent = .false.
            end if
            
         end if

       end function bf_layer_depends_on_neighbors


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if the neighboring bf_layers from bf_sublayer_i can
       !> all be removed
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param bf_sublayer_i
       !> bf_sublayer exchanging data with the neighboring buffer layers
       !
       !>@param bf_mainlayer_id
       !> cardinal coordinate of bf_sublayer_i
       !
       !>@return a_neighbor_remains
       !> logical stating whether the bf_sublayer_i should not be removed
       !> because one of its neighbors remains
       !--------------------------------------------------------------
       function does_a_neighbor_remains(
     $     this, bf_sublayer_i, bf_mainlayer_id)
     $     result(a_neighbor_remains)

         implicit none

         class(bf_interface)       , intent(in) :: this
         type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
         integer         , optional, intent(in) :: bf_mainlayer_id
         logical                                :: a_neighbor_remains


         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 1
         if(bf_sublayer_i%can_exchange_with_neighbor1()) then
            
            if(present(bf_mainlayer_id)) then
               a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $              1, bf_sublayer_i, bf_mainlayer_id)
            else
               a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $              1, bf_sublayer_i)
            end if

         else
            a_neighbor_remains = .false.

         end if
         

         !determine if the buffer layer is sharing grid points
         !with its neighborign buffer layers of type 2
         if(.not.a_neighbor_remains) then
            
            if(bf_sublayer_i%can_exchange_with_neighbor2()) then
               
               if(present(bf_mainlayer_id)) then
                  a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $                 2, bf_sublayer_i, bf_mainlayer_id)
               else
                  a_neighbor_remains = this%border_interface%does_a_neighbor_remains(
     $                 2, bf_sublayer_i)
               end if

            else
               a_neighbor_remains = .false.
            end if
            
         end if

       end function does_a_neighbor_remains


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if the buffer layers at the interface between main
       !> layers (N<->E, N<->W, S<->E, S<->W) have boundary layers
       !> that cannot be computed due to a lack of grid points (after
       !> the neighboring boundary layers updated their grid points
       !> w/o updating its allocation)
       !
       !> @date
       !> 18_12_2014 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param interior_x_map
       !> array encapsulating the coordinates along the x-axis for
       !> the interior computational domain
       !
       !>@param interior_y_map
       !> array encapsulating the coordinates along the y-axis for
       !> the interior computational domain
       !
       !>@param interior_nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !--------------------------------------------------------------
       subroutine resolve_bc_overlap_conflicts(
     $     this,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes


         integer                     :: k
         type(bf_mainlayer), pointer :: mainlayer
         integer                     :: nb_bf_layers
         type(bf_sublayer) , pointer :: bf_sublayer_ptr
         integer                     :: l

         
         !loop over the main layers
         do k=1, size(this%mainlayer_pointers,1)

            !extract the main layer cooresponding to the
            !cardinal coordinate k
            mainlayer => this%get_mainlayer(k)
            
            if(associated(mainlayer)) then

               !extract the number of buffer layers
               !inside the mainlayer
               nb_bf_layers = mainlayer%get_nb_sublayers()

               !loop over the buffer layers contained in
               !the main layer
               bf_sublayer_ptr => mainlayer%get_head_sublayer()

               do l=1, nb_bf_layers

                  !resolve the boundary overlap conflicts b/w
                  !the interface buffer layer if the current
                  !buffer layer is indeed at the interface 
                  !b/w main layers
                  call resolve_bc_overlap_conflicts_for_bf_layer(
     $                 this,
     $                 k,
     $                 bf_sublayer_ptr,
     $                 interior_x_map,
     $                 interior_y_map,
     $                 interior_nodes)

                  bf_sublayer_ptr => bf_sublayer_ptr%get_next()

               end do

            end if

         end do

       end subroutine resolve_bc_overlap_conflicts


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if whether the buffer layer is at the interface b/w
       !> main layers (N<->E, N<->W, S<->E, S<->W) and whether the
       !> buffer layer has boundary layers that cannot be computed due
       !> to a lack of grid points (after the neighboring boundary
       !> layers updated their grid points w/o updating its allocation).
       !> If there is a lack of grid points, the buffer layer is
       !> reallocated to have all the necessary grid points from its
       !> neighbors
       !
       !> @date
       !> 18_12_2014 - initial version - J.L. Desmarais
       !
       !> @param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !> @param bf_localization
       !> cardinal coordinate identifying to which buffer main layer
       !> the 'bf_layer_investigated' belongs
       !> 
       !> @param bf_sublayer_ptr
       !> buffer layer whose grid points at the interface b/w main
       !> layers are checked to see whether additional grid points are
       !> needed to compute the boundary conditions
       !
       !>@param interior_x_map
       !> array encapsulating the coordinates along the x-axis for
       !> the interior computational domain
       !
       !>@param interior_y_map
       !> array encapsulating the coordinates along the y-axis for
       !> the interior computational domain
       !
       !>@param interior_nodes
       !> table encapsulating the data of the grid points of the
       !> interior domain
       !--------------------------------------------------------------
       subroutine resolve_bc_overlap_conflicts_for_bf_layer(
     $     this,
     $     bf_localization,
     $     bf_sublayer_ptr,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes)

         implicit none

         class(bf_interface)             , intent(inout) :: this
         integer                         , intent(in)    :: bf_localization
         type(bf_sublayer), pointer      , intent(inout) :: bf_sublayer_ptr
         real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
         real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes


         integer(ikind), dimension(2,2) :: bf_alignment_update
         logical                        :: should_be_reallocated
         logical                        :: update_alignment
         integer(ikind)                 :: x_border
         integer(ikind), dimension(2)   :: sizes


         !get the current alignment of the buffer layer
         bf_alignment_update = bf_sublayer_ptr%get_alignment_tab()
         

         !first assume that the buffer layer should not
         !be reallocated
         should_be_reallocated = .false.

         
         !get the new alignmment of the buffer layer if modified
         select case(bf_localization)
           case(N)

              !check whether there are overlap conflicts with the
              !W neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor1(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [1,bc_size+1],
     $             left,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,1) = x_border+bc_size
                 should_be_reallocated = .true.
              end if


              !check whether there are overlap conflicts with the
              !E neighboring buffer layers
              sizes = bf_sublayer_ptr%get_sizes()

              call resolve_bc_overlap_conflict_with_neighbor2(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [sizes(1),bc_size+1],
     $             right,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,2) = x_border-bc_size
                 should_be_reallocated = .true.
              end if
           

           case(S)

              !get the extent of the arrays in the buffer layer
              sizes = bf_sublayer_ptr%get_sizes()

              !check whether there are overlap conflicts with the
              !W neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor1(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [1,sizes(2)-bc_size],
     $             left,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,1) = x_border+bc_size
                 should_be_reallocated = .true.
              end if


              !check whether there are overlap conflicts with the
              !E neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor2(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [sizes(1),sizes(2)-bc_size],
     $             right,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,2) = x_border-bc_size
                 should_be_reallocated = .true.
              end if


           case(E)

              !get the extent of the arrays in the buffer layer
              sizes = bf_sublayer_ptr%get_sizes()

              !check whether there are overlap conflicts with the
              !S neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor1(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [sizes(1),bc_size+1],
     $             right,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,2) = x_border-bc_size
                 should_be_reallocated = .true.
              end if


              !check whether there are overlap conflicts with the
              !N neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor2(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [sizes(1),sizes(2)-bc_size],
     $             right,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,2) =
     $                max(bf_alignment_update(1,2),x_border-bc_size)
                 should_be_reallocated = .true.
              end if


           case(W)

              !check whether there are overlap conflicts with the
              !S neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor1(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [1,bc_size+1],
     $             left,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,1) = x_border+bc_size
                 should_be_reallocated = .true.
              end if


              !check whether there are overlap conflicts with the
              !N neighboring buffer layers
              call resolve_bc_overlap_conflict_with_neighbor2(
     $             this,
     $             bf_sublayer_ptr,
     $             bf_localization,
     $             [1,sizes(2)-bc_size],
     $             left,
     $             update_alignment,
     $             x_border)

              if(update_alignment) then
                 bf_alignment_update(1,1) =
     $                max(bf_alignment_update(1,1),x_border+bc_size)
                 should_be_reallocated = .true.
              end if

           case default
              call error_mainlayer_id(
     $             'bf_interface_class',
     $             'resolve_bc_overlap_conflicts_for_bf_layer',
     $             bf_localization)

         end select


         !if the buffer layer should be reallocated to get
         !the grid points needed, the bf_interface updates
         !its allocation
         if(should_be_reallocated) then
            call reallocate_sublayer(
     $           this,
     $           bf_sublayer_ptr,
     $           interior_x_map,
     $           interior_y_map,
     $           interior_nodes,
     $           bf_alignment_update)
         end if

       end subroutine resolve_bc_overlap_conflicts_for_bf_layer


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if whether the buffer layer is at the interface b/w
       !> main layers (N<->E, N<->W, S<->E, S<->W) and whether the
       !> buffer layer has boundary layers that cannot be computed due
       !> to a lack of grid points (after the neighboring boundary
       !> layers updated their grid points w/o updating its allocation).
       !> If there is a lack of grid points, the new x_border is given
       !> in the general frame
       !
       !> @date
       !> 18_12_2014 - initial version - J.L. Desmarais
       !
       !> @param bf_sublayer_ptr
       !> buffer layer whose grid points at the interface b/w main
       !> layers are checked to see whether additional grid points are
       !> needed to compute the boundary conditions
       !
       !> @param bf_localization
       !> cardinal coordinate identifying to which buffer main layer
       !> the 'bf_layer_investigated' belongs
       !
       !>@param l_coords
       !> coordinates of the grid point, expressed in the local frame of
       !> the buffer layer, and located at the edge of the buffer layer
       !> 'bf_sublayer_ptr' to decide whether there are grid points
       !> missing for the computation of the boundary conditions
       !
       !>@param side
       !> direction in which the new x_border should be searched
       !
       !>@param update_alignment
       !> logical stating whether the alignment of the buffer layer
       !> 'bf_layer_ptr' should be updated to incorporate the missing
       !> grid points
       !
       !>@param x_border
       !> integer stating what should be the new border of the buffer
       !> layer (general coordinate of the edge of the buffer layer not
       !> its new alignment)
       !--------------------------------------------------------------
       subroutine resolve_bc_overlap_conflict_with_neighbor1(
     $     this,
     $     bf_sublayer_ptr,
     $     bf_localization,
     $     l_coords,
     $     side,
     $     update_alignment,
     $     x_border)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         type(bf_sublayer)           , intent(in)  :: bf_sublayer_ptr
         integer                     , intent(in)  :: bf_localization
         integer(ikind), dimension(2), intent(in)  :: l_coords
         logical                     , intent(in)  :: side
         logical                     , intent(out) :: update_alignment
         integer(ikind)              , intent(out) :: x_border         

         integer :: neighbor_type

         update_alignment = .false.
         x_border = 1

         !if the buffer layer exchanges with the neighboring
         !buffer layers of type 1
         if(bf_sublayer_ptr%can_exchange_with_neighbor1()) then

            neighbor_type = 1

            call resolve_bc_overlap_conflict_with_neighbor(
     $           this,
     $           bf_sublayer_ptr,
     $           bf_localization,
     $           l_coords,
     $           side,
     $           neighbor_type,
     $           update_alignment,
     $           x_border)

         end if

       end subroutine resolve_bc_overlap_conflict_with_neighbor1


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if whether the buffer layer is at the interface b/w
       !> main layers (N<->E, N<->W, S<->E, S<->W) and whether the
       !> buffer layer has boundary layers that cannot be computed due
       !> to a lack of grid points (after the neighboring boundary
       !> layers updated their grid points w/o updating its allocation).
       !> If there is a lack of grid points, the new x_border is given
       !> in the general frame
       !
       !> @date
       !> 18_12_2014 - initial version - J.L. Desmarais
       !
       !> @param bf_sublayer_ptr
       !> buffer layer whose grid points at the interface b/w main
       !> layers are checked to see whether additional grid points are
       !> needed to compute the boundary conditions
       !
       !> @param bf_localization
       !> cardinal coordinate identifying to which buffer main layer
       !> the 'bf_layer_investigated' belongs
       !
       !>@param l_coords
       !> coordinates of the grid point, expressed in the local frame of
       !> the buffer layer, and located at the edge of the buffer layer
       !> 'bf_sublayer_ptr' to decide whether there are grid points
       !> missing for the computation of the boundary conditions
       !
       !>@param side
       !> direction in which the new x_border should be searched
       !
       !>@param update_alignment
       !> logical stating whether the alignment of the buffer layer
       !> 'bf_layer_ptr' should be updated to incorporate the missing
       !> grid points
       !
       !>@param x_border
       !> integer stating what should be the new border of the buffer
       !> layer (general coordinate of the edge of the buffer layer not
       !> its new alignment)
       !--------------------------------------------------------------
       subroutine resolve_bc_overlap_conflict_with_neighbor2(
     $     this,
     $     bf_sublayer_ptr,
     $     bf_localization,
     $     l_coords,
     $     side,
     $     update_alignment,
     $     x_border)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         type(bf_sublayer)           , intent(in)  :: bf_sublayer_ptr
         integer                     , intent(in)  :: bf_localization
         integer(ikind), dimension(2), intent(in)  :: l_coords
         logical                     , intent(in)  :: side
         logical                     , intent(out) :: update_alignment
         integer(ikind)              , intent(out) :: x_border         

         integer       :: neighbor_type

         update_alignment = .false.
         x_border = 1

         !if the buffer layer exchanges with the neighboring buffer layers
         !of type 2
         if(bf_sublayer_ptr%can_exchange_with_neighbor2()) then

            neighbor_type = 2

            call resolve_bc_overlap_conflict_with_neighbor(
     $           this,
     $           bf_sublayer_ptr,
     $           bf_localization,
     $           l_coords,
     $           side,
     $           neighbor_type,
     $           update_alignment,
     $           x_border)

         end if

       end subroutine resolve_bc_overlap_conflict_with_neighbor2


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> check if whether the buffer layer is at the interface b/w
       !> main layers (N<->E, N<->W, S<->E, S<->W) and whether the
       !> buffer layer has boundary layers that cannot be computed due
       !> to a lack of grid points (after the neighboring boundary
       !> layers updated their grid points w/o updating its allocation).
       !> If there is a lack of grid points, the new x_border is given
       !> in the general frame
       !
       !> @date
       !> 18_12_2014 - initial version - J.L. Desmarais
       !
       !> @param bf_sublayer_ptr
       !> buffer layer whose grid points at the interface b/w main
       !> layers are checked to see whether additional grid points are
       !> needed to compute the boundary conditions
       !
       !> @param bf_localization
       !> cardinal coordinate identifying to which buffer main layer
       !> the 'bf_layer_investigated' belongs
       !
       !>@param l_coords
       !> coordinates of the grid point, expressed in the local frame of
       !> the buffer layer, and located at the edge of the buffer layer
       !> 'bf_sublayer_ptr' to decide whether there are grid points
       !> missing for the computation of the boundary conditions
       !
       !>@param side
       !> direction in which the new x_border should be searched
       !
       !>@param neighbor_type
       !> type of the neighboring buffer layers tested
       !
       !>@param update_alignment
       !> logical stating whether the alignment of the buffer layer
       !> 'bf_layer_ptr' should be updated to incorporate the missing
       !> grid points
       !
       !>@param x_border
       !> integer stating what should be the new border of the buffer
       !> layer (general coordinate of the edge of the buffer layer not
       !> its new alignment)
       !--------------------------------------------------------------
       subroutine resolve_bc_overlap_conflict_with_neighbor(
     $     this,
     $     bf_sublayer_ptr,
     $     bf_localization,
     $     l_coords,
     $     side,
     $     neighbor_type,
     $     update_alignment,
     $     x_border)

         implicit none

         class(bf_interface)         , intent(in)  :: this
         type(bf_sublayer)           , intent(in)  :: bf_sublayer_ptr
         integer                     , intent(in)  :: bf_localization
         integer(ikind), dimension(2), intent(in)  :: l_coords
         integer                     , intent(in)  :: neighbor_type
         logical                     , intent(in)  :: side
         logical                     , intent(out) :: update_alignment
         integer(ikind)              , intent(out) :: x_border

         integer       , dimension(2) :: match_table
         integer(ikind), dimension(2) :: start_grdpt
         logical                      :: err

         !if there is a bc_interior_pt at the left edge of
         !the buffer layer, it will not be possible to compute
         !the boundary condition and so its alignment should be
         !updated                 
         if(bf_sublayer_ptr%is_bc_interior_pt(l_coords(1),l_coords(2))) then

            !get the table converting local to general coords
            match_table = bf_sublayer_ptr%get_general_to_local_coord_tab()

            !convert the local coordinates of the grdpt [i,j]
            !into general coordinates
            start_grdpt(1) = l_coords(1)+match_table(1)
            start_grdpt(2) = l_coords(2)+match_table(2)
            
            !ask the neighboring buffer layer to have the
            !x_border identifying how far the buffer layer
            !should be reallocated to the left
            !the coordinate obtained is expressed in the 
            !general frame (interior + bf_layers)
            x_border = this%border_interface%ask_neighbors_for_bc_overlap(
     $           bf_localization,
     $           neighbor_type,
     $           start_grdpt,
     $           side,
     $           err)

            if(.not.(err.eqv.BF_SUCCESS)) then
               print '(''bf_interface_class'')'
               print '(''resolve_bc_overlap_conflict_with_neighbor'')'
               print '(''*****************************************'')'
               print '(''did not manage to get x_border'')'
               print '(''*****************************************'')'
               print '(''bf_localization: '',I2)', bf_localization
               print '(''l_coords: '',2I2)', l_coords
               print '(''neighbor_type: '',I2)', neighbor_type
               print '(''side: '',I2)', side
               print '(''match_table: '',2I5)', match_table
               print '(''start_grdpt: '',2I5)', start_grdpt
               print '()'
               stop ''
            end if

            update_alignment = .true.
            
         end if

       end subroutine resolve_bc_overlap_conflict_with_neighbor

       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the content of the interface on external binary files
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param suffix_nodes
       !> suffix for the name of the files storing the nodes content
       !
       !>@param suffix_grdid
       !> suffix for the name of the files storing the grdpts_id content
       !
       !>@param suffix_sizes
       !> suffix for the name of the files storing the profile of the 
       !> nodes and grdpts_id arrays
       !
       !>@param suffix_nb_sublayers_max
       !> suffix for the name of the files storing the maximum number of
       !> sublayers per main buffer layer
       !--------------------------------------------------------------
       subroutine print_binary(
     $     this,
     $     suffix_x_map,
     $     suffix_y_map,
     $     suffix_nodes,
     $     suffix_grdid,
     $     suffix_sizes,
     $     suffix_nb_sublayers_max,
     $     timedev)

         implicit none

         class(bf_interface), intent(in) :: this
         character(*)       , intent(in) :: suffix_x_map
         character(*)       , intent(in) :: suffix_y_map
         character(*)       , intent(in) :: suffix_nodes
         character(*)       , intent(in) :: suffix_grdid
         character(*)       , intent(in) :: suffix_sizes
         character(*)       , intent(in) :: suffix_nb_sublayers_max
         logical, optional  , intent(in) :: timedev

         logical :: timedev_op

         integer           :: i
         integer           :: nb_sublayers_max
         character(len=18) :: filename_format
         character(len=28) :: nb_sublayers_filename

         if(present(timedev)) then
            timedev_op = timedev
         else
            timedev_op = .false.
         end if
                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files and 
         !determine the maximum number of sublayers
         nb_sublayers_max = 0
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_binary(
     $              suffix_x_map,
     $              suffix_y_map,
     $              suffix_nodes,
     $              suffix_grdid,
     $              suffix_sizes,
     $              timedev=timedev_op)

               nb_sublayers_max = max(
     $              nb_sublayers_max,
     $              this%mainlayer_pointers(i)%get_nb_sublayers())
            end if            
         end do


         !print the maximum number of sublayers in an output
         !binary file
         write(filename_format,
     $        '(''(A12,A'',I2,'')'')')
     $        len(suffix_nb_sublayers_max)

         write(nb_sublayers_filename, filename_format)
     $        'sublayers_nb',
     $        suffix_nb_sublayers_max

         call print_nb_sublayers_max(
     $        nb_sublayers_filename, nb_sublayers_max)

        end subroutine print_binary


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the maxmimum number of sublayers per main layer
       !
       !> @date
       !> 11_04_2013 - initial version - J.L. Desmarais
       !
       !>@param filename
       !> name of the file for storing the maximum number of sublayer
       !> per mainlayer
       !
       !>@param nb_sublayers
       !> maximum number of sublayers per main layer
       !--------------------------------------------------------------
       subroutine print_nb_sublayers_max(filename, nb_sublayers)

          implicit none

          character(*), intent(in) :: filename
          integer     , intent(in) :: nb_sublayers

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers_max


       !> @author
       !> Julien L. Desmarais
       !
       !> @brief
       !> print the content of the interface on external netcdf files
       !
       !> @date
       !> 11_07_2013 - initial version - J.L. Desmarais
       !
       !>@param this
       !> bf_interface object encapsulating the buffer layers
       !> around the interior domain and subroutines to synchronize
       !> the data between them
       !
       !>@param timestep_written
       !> integer identifying the timestep written
       !
       !>@param name_var
       !> table with the short name for the governing variables saved
       !> in the netcdf file
       !
       !>@param longname_var
       !> table with the long name for the governing variables saved
       !> in the netcdf file
       !
       !>@param unit_var
       !> table with the units of the governing variables saved
       !> in the netcdf file
       !
       !>@param time
       !> time corresponding to the data for the grdpts and the nodes
       !--------------------------------------------------------------
       subroutine print_netcdf(
     $     this,
     $     timestep_written,
     $     name_var,
     $     longname_var,
     $     unit_var,
     $     time)

         implicit none

         class(bf_interface)        , intent(in) :: this
         integer                    , intent(in) :: timestep_written
         character(*), dimension(ne), intent(in) :: name_var
         character(*), dimension(ne), intent(in) :: longname_var
         character(*), dimension(ne), intent(in) :: unit_var
         real(rkind)                , intent(in) :: time

         integer           :: i
                  

         !go through the buffer main layers and
         !print the content of each buffer layer
         !in seperate binary output files
         do i=1, size(this%mainlayer_pointers,1)

            if(this%mainlayer_pointers(i)%associated_ptr()) then
               
               call this%mainlayer_pointers(i)%print_netcdf(
     $              timestep_written,
     $              name_var,
     $              longname_var,
     $              unit_var,
     $              time)

            end if
         end do

        end subroutine print_netcdf


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !--------------------------------------------------------------
        subroutine allocate_before_timeInt(this)

          implicit none

          class(bf_interface), intent(inout) :: this

          integer :: i

          !go through the buffer main layers and
          !allocate the intermediate variables
          !needed to perform the integration step
          do i=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(i)%associated_ptr()) then
                
                call this%mainlayer_pointers(i)%allocate_before_timeInt()

             end if
          end do

        end subroutine allocate_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each main buffer layer contained in this bf_interface
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !--------------------------------------------------------------
        subroutine deallocate_after_timeInt(this)

          implicit none

          class(bf_interface), intent(inout) :: this

          integer :: i

          !go through the buffer main layers and
          !deallocate the intermediate variables
          !needed to perform the integration step
          do i=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(i)%associated_ptr()) then
                
                call this%mainlayer_pointers(i)%deallocate_after_timeInt()

             end if
          end do

        end subroutine deallocate_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the main layers
        !> contained in this bf_interface
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used)

          implicit none

          class(bf_interface), intent(inout) :: this
          type(td_operators) , intent(in)    :: td_operators_used
          real(rkind)        , intent(in)    :: t
          type(sd_operators) , intent(in)    :: s
          type(pmodel_eq)    , intent(in)    :: p_model
          type(bc_operators) , intent(in)    :: bc_used

          integer :: i

          !go through the buffer main layers and
          !compute the time derivatives
          do i=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(i)%associated_ptr()) then
                
                call this%mainlayer_pointers(i)%compute_time_dev(
     $               td_operators_used,
     $               t,s,p_model,bc_used)

             end if
          end do

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each main buffer layer contained in this bf_interface
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param dt
        !> integration time step
        !
        !>@param integration_step_nopt
        !> procedure performing the time integration
        !
        !>@param full
        !> logical to enforce whether the full domain is computed
        !> discarding the x_borders and y_borders supplied
        !> (important for the first integration step when the
        !> previous integration step is saved temporary in another
        !> array)
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt, full)

          implicit none

          class(bf_interface), intent(inout) :: this
          real(rkind)        , intent(in)    :: dt
          procedure(timeInt_step_nopt)       :: integration_step_nopt
          logical            , intent(in)    :: full

          integer :: i

          !go through the buffer main layers and
          !compute the integration step
          do i=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(i)%associated_ptr()) then
                
                call this%mainlayer_pointers(i)%compute_integration_step(
     $               dt, integration_step_nopt, full)

             end if
          end do

        end subroutine compute_integration_step          

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the integration borders of the buffer layer
        !
        !> @date
        !> 04_11_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param added_sublayer
        !> buffer layer whose integration borders are updated
        !--------------------------------------------------------------
        subroutine update_integration_borders(this,added_sublayer)

          implicit none

          class(bf_interface), intent(inout) :: this
          type(bf_sublayer)  , intent(inout) :: added_sublayer

          call this%border_interface%update_integration_borders(added_sublayer)

        end subroutine update_integration_borders          

      end module bf_interface_class
