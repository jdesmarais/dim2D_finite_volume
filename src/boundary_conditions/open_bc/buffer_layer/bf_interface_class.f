      !> @file
      !> bf_interface_icr object augmented with procedures to remove
      !> the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_icr object augmented with procedures to remove
      !> the buffer layers
      !
      !> @date
      ! 26_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_class

        use bf_interface_dcr_class, only :
     $       bf_interface_dcr

        use bf_restart_module, only :
     $       get_restart_alignment

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nf90_operators_module, only :
     $       nf90_close_file

        use nf90_operators_read_module, only :
     $       nf90_open_file_for_reading,
     $       nf90_get_varid,
     $       nf90_get_var_model_nopt,
     $       nf90_read_borders

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: bf_interface
        

        !>@class bf_interface
        !> bf_interface_icr object augmented with procedures to remove
        !> the buffer layers
        !
        !>@param remove_inactivated_bf_layers
        !> remove the buffer layers that do not have activated grid-points
        !> at the interface with the interior domain
        !
        !>@param adapt_domain_extension
        !> adapt the configuration and extents of the domain extension
        !---------------------------------------------------------------
        type, extends(bf_interface_dcr) :: bf_interface

           contains

           procedure, pass :: restart

        end type bf_interface

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> re-initialize the buffer layers
        !
        !> @date
        !> 29_03_2015 - initial version - J.L. Desmarais
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
                bf_index_format = floor(log10(real(i))+1.0)

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
                bf_sublayer_restart => this%allocate_sublayer(
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

             end do
          end do          

        end subroutine restart

      end module bf_interface_class
