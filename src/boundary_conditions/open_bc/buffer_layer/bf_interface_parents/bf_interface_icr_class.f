      !> @file
      !> module encapsulating the bf_interface_icr object.
      !> bf_interface_grdpts_id_update augmented with procedures
      !> detecting how the domain extension should be increased
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_grdpts_id_update augmented with procedures
      !> detecting how the domain extension should be increased
      !
      !> @date
      ! 24_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_icr_class

        use bf_increase_coords_module, only :
     $       get_mainlayer_coord

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_layer_bc_sections_overlap_module, only :
     $       determine_edge_grdpts_computed,
     $       determine_corner_or_anti_corner_grdpts_computed

        use bf_layer_errors_module, only :
     $       error_bc_section_type

        use bf_sublayer_class, only :
     $       bf_sublayer

        use icr_interface_class, only :
     $       icr_interface

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,     
     $       
     $       dct_icr_distance,
     $       
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       cptnot_type,
     $       
     $       BF_SUCCESS

        use parameters_constant, only :
     $       N,S,E,W,
     $       no_interior_mainlayer,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type

        use parameters_input, only :
     $       nx,ny,ne,
     $       adapt_N_choice,
     $       adapt_S_choice,
     $       adapt_E_choice,
     $       adapt_W_choice,
     $       obc_crenel_removal_ac

        use parameters_constant, only :
     $       interior,
     $       fixed_domain_choice,
     $       adapt_domain_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only : 
     $       pmodel_eq

        implicit none

        private
        public :: bf_interface_icr

        
        !> @class bf_interface_icr
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !> @param analyze_and_update_boundary
        !> analyze and update the boundary identified by its
        !> mainlayer_id
        !
        !> @param adapt_domain_extension
        !> adapt the configuration and extents of the domain
        !> extension by increasing the buffer layers where
        !> bc_interior_pt are activated
        !  
        !> @param get_interior_bcs_mainlayer_id
        !> determine the cardinal coordinate (N,S,E,W) identifying
        !> to which interiro boundary a bc_sction is belonging
        !
        !> @param analyze_interior_bc_sections
        !> analyze the bc_interior_pt of an interior bc_section
        !> and stage for update the ones triggered by nodes 
        !> activated in the interior domain
        !
        !> @param analyze_bf_layer_bc_sections
        !> analyze the bc_interior_pt of a buffer layer bc_section
        !> and stage for update the ones triggered by nodes 
        !> activated
        !  
        !> @param stage
        !> stage a bc_interior_pt for update
        !
        !> @param analyze_bc_section_edge_y
        !> analyze a bc_section of type edge_y (N_edge_type or
        !> S_edge_type)
        !
        !> @param analyze_bc_section_edge_x
        !> analyze a bc_section of type edge_x (E_edge_type or
        !> W_edge_type)
        !
        !> @param analyze_bc_section_xy
        !> analyze a bc_section of type square
        !> (NE_corner_type, NW_crner_type, SE_corner_type,
        !>  SW_corner_type, NE_edge_type, NW_edge_type,
        !>  SW_edge_type, SE_edge_type)
        !
        !> @param analyze_bc_section_square_bounds
        !> determine the grid-points triggering the activation of
        !> the bc_interior_pt in a square-like bc_section
        !> (NE_corner_type, NW_crner_type, SE_corner_type,
        !>  SW_corner_type, NE_edge_type, NW_edge_type,
        !>  SW_edge_type, SE_edge_type)
        !
        !> @param analyze_bc_section_square_xy
        !> check the gridpoints triggering the activation of the
        !> bc_interior_pt in a square-like bc_section
        !
        !> @param analyze_bc_section
        !> analyze the bc_interior_pt of an interior or buffer
        !> layer bc_section and stage for update the ones
        !> triggered by the activation of upstream grid-points
        !
        !> @param is_interior_node_activated
        !> determine whether a grid-point belonging to the interior
        !> domain is activated
        !
        !> @param is_bf_layer_node_activated
        !> determine whether a grid-point belonging to a buffer
        !> layer is activated
        !
        !> @param is_node_activated
        !> determine whether a grid-point is activated
        !
        !> @param update_entire_edge
        !> check whether the entire edge should be updated
        !> to prevent boundary crenel
        !-------------------------------------------------------------
        type, extends(bf_interface_coords) :: bf_interface_icr

          contains

          ! adapt the extent of the domain extension
          procedure,   pass :: analyze_and_update_boundary
          procedure,   pass :: adapt_domain_extension


          ! procedure to know which bc_interior_pt are 
          ! activated in the interior_bc_sections and
          ! the buffer layers
          procedure, nopass :: get_interior_bcs_mainlayer_id
          procedure,   pass :: analyze_interior_bc_sections
          procedure,   pass :: analyze_bf_layer_bc_sections
          

          ! procedures to know which bc_interior_pt are
          ! activated in a bc_section
          procedure,   pass :: stage
          procedure,   pass :: update_entire_edge
          procedure,   pass :: analyze_bc_section_edge_y
          procedure,   pass :: analyze_bc_section_edge_x
          procedure,   pass :: analyze_bc_section_xy
          procedure, nopass :: analyze_bc_section_square_bounds
          procedure,   pass :: analyze_bc_section_square_xy
          procedure,   pass :: analyze_bc_section


          ! procedures to know whether a grid-point is activated
          procedure, nopass :: is_interior_node_activated
          procedure,   pass :: is_bf_layer_node_activated
          procedure,   pass :: is_node_activated

        end type bf_interface_icr

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> adapt the extent of the domain extension by updating
        !> the bc_interior_pt at the edge of teh computational
        !> domain
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !>@param interior_nodes1
        !> nodes of the interior domain at t
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> timestep
        !
        !>@param interior_bc_sections
        !> boundary sections for the interior domain
        !--------------------------------------------------------------
        subroutine adapt_domain_extension(
     $       this,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes0,
     $       interior_nodes1,
     $       p_model,
     $       t,
     $       dt,
     $       interior_bc_sections)

          implicit none

          class(bf_interface_icr)                    , intent(inout) :: this
          real(rkind), dimension(nx)                 , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)                 , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes1
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: interior_bc_sections

          
          type(icr_interface) :: icr_interface_used
          integer             :: mainlayer_id


          call icr_interface_used%ini()


          ! check whether the South boundary should be updated
          if(adapt_S_choice.eq.adapt_domain_choice) then
             
             mainlayer_id = S

             call analyze_and_update_boundary(
     $            this,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            p_model,
     $            t,
     $            dt,
     $            interior_bc_sections,
     $            mainlayer_id)

          end if


          ! check whether the West boundary should be updated
          if(adapt_W_choice.eq.adapt_domain_choice) then
             
             mainlayer_id = W

             call analyze_and_update_boundary(
     $            this,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            p_model,
     $            t,
     $            dt,
     $            interior_bc_sections,
     $            mainlayer_id)

          end if


          ! check whether the East boundary should be updated
          if(adapt_E_choice.eq.adapt_domain_choice) then
             
             mainlayer_id = E

             call analyze_and_update_boundary(
     $            this,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            p_model,
     $            t,
     $            dt,
     $            interior_bc_sections,
     $            mainlayer_id)

          end if


          ! check whether the North boundary should be updated
          if(adapt_N_choice.eq.adapt_domain_choice) then
             
             mainlayer_id = N

             call analyze_and_update_boundary(
     $            this,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes0,
     $            interior_nodes1,
     $            p_model,
     $            t,
     $            dt,
     $            interior_bc_sections,
     $            mainlayer_id)

          end if


          ! finalize the update of the boundary by committing the 
          ! paths of the icr_interface that were not processed
          call icr_interface_used%finalize_domain_increase(
     $         this,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)


          ! uniformize the size of the buffer layers at the interface
          ! between the main layers
          call this%uniformize_mainlayer_interfaces(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1)

        end subroutine adapt_domain_extension


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze and update a boundary (N,S,E or W)
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes0
        !> nodes of the interior domain at t-dt
        !
        !>@param interior_nodes1
        !> nodes of the interior domain at t
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param dt
        !> time step
        !
        !>@param interior_bc_sections
        !> boundary sections for the interior domain
        !--------------------------------------------------------------
        subroutine analyze_and_update_boundary(
     $     this,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     p_model,
     $     t,
     $     dt,
     $     interior_bc_sections,
     $     mainlayer_id)

          implicit none

          class(bf_interface_icr)                    , intent(inout) :: this
          type(icr_interface)                        , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)                 , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)                 , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           , intent(in)    :: interior_nodes1
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)                                , intent(in)    :: t
          real(rkind)                                , intent(in)    :: dt
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: interior_bc_sections
          integer                                    , intent(in)    :: mainlayer_id

          
          !analyze the interior_bc_sections belonging to the
          !boundary investigated
          call analyze_interior_bc_sections(
     $         this,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1,
     $         p_model,
     $         mainlayer_id,
     $         interior_bc_sections)

          !analyze the buffer layer bc_sections belonging to
          !the boundary investigated
          call analyze_bf_layer_bc_sections(
     $         this,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1,
     $         p_model,
     $         mainlayer_id)

          !apply the update operations on the buffer layers
          !of the boundary investigated
          call icr_interface_used%commit(
     $         mainlayer_id,
     $         this,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

        end subroutine analyze_and_update_boundary


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer computed
        !--------------------------------------------------------------
        subroutine analyze_bf_layer_bc_sections(
     $       this,
     $       icr_interface_used,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       mainlayer_id)

          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer                         , intent(in)    :: mainlayer_id


          logical                                     :: interior_domain
          integer                                     :: nb_sublayers
          integer                                     :: k,m
          type(bf_sublayer), pointer                  :: bf_sublayer_ptr
          integer(ikind), dimension(2)                :: match_table
          integer(ikind), dimension(:,:), allocatable :: bf_layer_bc_sections

          integer, dimension(5) :: modified_bc_section
          integer               :: bcs_mainlayer_id

          interior_domain = .false.


          ! determine the number of buffer layers in the main layer
          ! to estimate the size of the loop over the buffer layers
          !------------------------------------------------------------
          nb_sublayers = this%mainlayer_pointers(mainlayer_id)%get_nb_sublayers()

          if(nb_sublayers.gt.0) then


             ! analyze each buffer layer by extracting its bc_sections
             ! and then analyzing the bc_interior_pt in each bc_section
             !------------------------------------------------------------
             !k: loop index for the buffer layers
             !m: loop index for the bc_section of the buffer layer
             !------------------------------------------------------------
             bf_sublayer_ptr => this%mainlayer_pointers(mainlayer_id)%get_head_sublayer()

             do k=1, nb_sublayers

                match_table = bf_sublayer_ptr%get_general_to_local_coord_tab()

                call bf_sublayer_ptr%get_dct_sections(bf_layer_bc_sections)

                if(allocated(bf_layer_bc_sections)) then

                   do m=1, size(bf_layer_bc_sections,2)

                      !get the interior mainlayer_id to which the bc_section
                      !is belonging
                      modified_bc_section(2) = bf_layer_bc_sections(2,m)+match_table(1)
                      modified_bc_section(3) = bf_layer_bc_sections(3,m)+match_table(2)
                      bcs_mainlayer_id = get_interior_bcs_mainlayer_id(modified_bc_section)

                      !verify that this interior mainlayer can be adapted
                      !otherwise refuse to analyze this boundary section
                      if(  (bcs_mainlayer_id.eq.no_interior_mainlayer).or.(.not.(
     $                     ((bcs_mainlayer_id.eq.N).and.(adapt_N_choice.eq.fixed_domain_choice)).or.
     $                     ((bcs_mainlayer_id.eq.S).and.(adapt_S_choice.eq.fixed_domain_choice)).or.
     $                     ((bcs_mainlayer_id.eq.E).and.(adapt_E_choice.eq.fixed_domain_choice)).or.
     $                     ((bcs_mainlayer_id.eq.W).and.(adapt_W_choice.eq.fixed_domain_choice))))) then

                         call analyze_bc_section(
     $                        this,
     $                        icr_interface_used,
     $                        bf_sublayer_ptr,
     $                        match_table,
     $                        interior_x_map,
     $                        interior_y_map,
     $                        interior_nodes,
     $                        p_model,
     $                        interior_domain,
     $                        bf_layer_bc_sections(:,m))

                      end if

                   end do

                   deallocate(bf_layer_bc_sections)

                end if

                bf_sublayer_ptr => bf_sublayer_ptr%get_next()

             end do

          end if

        end subroutine analyze_bf_layer_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer computed
        !
        !>@param interior_bc_sections
        !> boundary sections for the interior domain
        !--------------------------------------------------------------
        subroutine analyze_interior_bc_sections(
     $       this,
     $       icr_interface_used,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       mainlayer_id,
     $       interior_bc_sections)

          implicit none

          class(bf_interface_icr)                      , intent(in)    :: this
          type(icr_interface)                          , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)                   , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)                   , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          type(pmodel_eq)                              , intent(in)    :: p_model
          integer                                      , intent(in)    :: mainlayer_id
          integer(ikind), dimension(:,:)  , allocatable, intent(in)    :: interior_bc_sections


          integer                      :: k
          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: match_table
          logical                      :: interior_domain
          integer                      :: bcs_mainlayer_id

          interior_domain = .true.


          !> the interior_bc_sections are analyzed and only
          !> the ones matching the mainlayer_id requested
          !> have their bc_interior_pt analyzed
          if(allocated(interior_bc_sections)) then

             do k=1, size(interior_bc_sections,2)

                bcs_mainlayer_id = get_interior_bcs_mainlayer_id(
     $               interior_bc_sections(:,k))

                if(bcs_mainlayer_id.eq.mainlayer_id) then

                   call analyze_bc_section(
     $                  this,
     $                  icr_interface_used,
     $                  bf_sublayer_ptr,
     $                  match_table,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes,
     $                  p_model,
     $                  interior_domain,
     $                  interior_bc_sections(:,k))

                end if                   

             end do

          end if

        end subroutine analyze_interior_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the cardinal coordinate of the main layer 
        !> to which the interior bc_section is belonging
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section
        !
        !>@return mainlayer_id
        !> cardinal coordinate of the main layer to which the
        !> bc_section is belonging
        !--------------------------------------------------------------
        function get_interior_bcs_mainlayer_id(
     $       bc_section)
     $       result(mainlayer_id)

          implicit none

          integer(ikind), dimension(5), intent(in) :: bc_section
          integer                                  :: mainlayer_id


          if(bc_section(3).eq.(align_S-1)) then
             mainlayer_id = S

          else
             if(bc_section(3).eq.align_N) then
                mainlayer_id = N

             else
                
                if(bc_section(2).eq.(align_W-1)) then
                   mainlayer_id = W

                else
                   if(bc_section(2).eq.align_E) then
                      mainlayer_id = E

                   else
                      mainlayer_id = no_interior_mainlayer

                   end if
                end if
             end if
          end if

        end function get_interior_bcs_mainlayer_id


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_bc_section(
     $       this,
     $       icr_interface_used,
     $       bf_sublayer_ptr,
     $       match_table,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       interior_domain,
     $       bc_section)

          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr
          integer(ikind), dimension(2)    , intent(in)    :: match_table
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                         , intent(in)    :: interior_domain
          integer(ikind), dimension(5)    , intent(in)    :: bc_section


          select case(bc_section(1))

            case(N_edge_type,S_edge_type)

               call analyze_bc_section_edge_y(
     $              this,
     $              icr_interface_used,
     $              bf_sublayer_ptr,
     $              match_table,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              p_model,
     $              interior_domain,
     $              bc_section)
               
            case(E_edge_type,W_edge_type)
               
               call analyze_bc_section_edge_x(
     $              this,
     $              icr_interface_used,
     $              bf_sublayer_ptr,
     $              match_table,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              p_model,
     $              interior_domain,
     $              bc_section)

            case(NE_corner_type, NW_corner_type,
     $           SE_corner_type, SW_corner_type,
     $           NE_edge_type  , NW_edge_type,
     $           SE_edge_type  , SW_edge_type)

               call analyze_bc_section_square_xy(
     $              this,
     $              icr_interface_used,
     $              bf_sublayer_ptr,
     $              match_table,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              p_model,
     $              interior_domain,
     $              bc_section)

            case default

               call error_bc_section_type(
     $              'bf_interface_icr_class',
     $              'analyze_bc_section',
     $              bc_section(1))

          end select
               
        end subroutine analyze_bc_section


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_bc_section_square_xy(
     $       this,
     $       icr_interface_used,
     $       bf_sublayer_ptr,
     $       match_table,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       interior_domain,
     $       bc_section)

          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr
          integer(ikind), dimension(2)    , intent(in)    :: match_table
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                         , intent(in)    :: interior_domain
          integer(ikind), dimension(5)    , intent(in)    :: bc_section


          integer(ikind), dimension(2,2,2) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: activated_grdpts
          integer                          :: nb_activated_grdpts
          logical                          :: no_activation


          !> determine the edges for the analyze of the grid-points
          !> and the bc_interior_pt activated
          call analyze_bc_section_square_bounds(
     $         bc_section,
     $         analyzed_grdpts_bounds,
     $         activated_grdpts,
     $         nb_activated_grdpts,
     $         no_activation)
          

          !> check whether the bc_interior_pt are indeed activated
          call analyze_bc_section_xy(
     $         this,
     $         icr_interface_used,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         interior_domain,
     $         analyzed_grdpts_bounds,
     $         activated_grdpts,
     $         nb_activated_grdpts,
     $         no_activation)

        end subroutine analyze_bc_section_square_xy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param analyzed_grdpts_bounds
        !> bounds for the loops checking the grid-points triggering
        !> the update of the bc_interior_pt
        !
        !>@param activated_grdpts
        !> bc_interior_pt triggered by the activation of the
        !> grid-points
        !
        !>@param nb_activated_grdpts
        !> total number of bc_interior_pt triggered by the activation
        !> of the grid-points
        !
        !>@param no_activation
        !> determines whether the grid-points should be checked
        !> for the activation of the bc_interior_pt
        !--------------------------------------------------------------
        subroutine analyze_bc_section_square_bounds(
     $     bc_section,
     $     analyzed_grdpts_bounds,
     $     activated_grdpts,
     $     nb_activated_grdpts,
     $     no_activation)

          implicit none

          integer(ikind), dimension(5)    , intent(in)  :: bc_section
          integer(ikind), dimension(2,2,2), intent(out) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(out) :: activated_grdpts
          integer                         , intent(out) :: nb_activated_grdpts
          logical                         , intent(out) :: no_activation


          integer(ikind)               :: i_min
          integer(ikind)               :: j_min
          integer       , dimension(4) :: compute_point
          integer(ikind), dimension(2) :: pt


          no_activation = .false.


          !> determine which grid-points are computed in the 
          !> square (corner or anti-corner)
          call determine_corner_or_anti_corner_grdpts_computed(
     $         bc_section(4),
     $         bc_section(5),
     $         compute_point)

          
          i_min = bc_section(2)
          j_min = bc_section(3)


          !> determine the bounds for the activation and the
          !> bc_interior_pt activated
          select case(bc_section(1))
            case(NE_corner_type)

               if(compute_point(1).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min,j_min]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 
     $                 pt(1)-dct_icr_distance-1,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance/),
     $                 (/2,2,2/))

               end if


            case(NW_corner_type)

               if(compute_point(2).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min+1,j_min]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+dct_icr_distance+1,
     $                 pt(2)-dct_icr_distance/),
     $                 (/2,2,2/))

               end if


            case(SE_corner_type)

               if(compute_point(3).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min,j_min+1]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance-1,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1/),
     $                 (/2,2,2/))

               end if

            case(SW_corner_type)

               if(compute_point(4).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+dct_icr_distance+1,
     $                 pt(2)+dct_icr_distance,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1/),
     $                 (/2,2,2/))

               end if

            case(NE_edge_type)
               
               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt = [i_min,j_min]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+2,
     $                 pt(2)-dct_icr_distance,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance+1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+2/),
     $                 (/2,2,2/))

               end if

            case(NW_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min+1,j_min]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-2,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance+1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+2/),
     $                 (/2,2,2/))

               end if

            case(SE_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min,j_min+1]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-2,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance-1,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+2,
     $                 pt(2)+dct_icr_distance/),
     $                 (/2,2,2/))

               end if

            case(SW_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min+1,j_min+1]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-2,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance-1,
     $                 
     $                 pt(1)-2,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance/),
     $                 (/2,2,2/))

               end if
               
          end select

        end subroutine analyze_bc_section_square_bounds


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@param analyzed_grdpts_bounds
        !> bounds for the grid-points to be analyzed
        !
        !>@param activated_grdpts
        !> general coordinates of the activated bc_interior_pt if
        !> the grid-points analyzed are activated
        !
        !>@param nb_activated_grdpts
        !> number of bc_interior_pt activated
        !
        !>@param no_activation
        !> determines whether the grid-points should be checked
        !> for the activation of the bc_interior_pt
        !--------------------------------------------------------------
        subroutine analyze_bc_section_xy(
     $     this,
     $     icr_interface_used,
     $     bf_sublayer_ptr,
     $     match_table,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     interior_domain,
     $     analyzed_grdpts_bounds,
     $     activated_grdpts,
     $     nb_activated_grdpts,
     $     no_activation)

          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr
          integer(ikind), dimension(2)    , intent(in)    :: match_table
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                         , intent(in)    :: interior_domain
          integer(ikind), dimension(2,2,2), intent(in)    :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(in)    :: activated_grdpts
          integer                         , intent(in)    :: nb_activated_grdpts
          logical                         , intent(in)    :: no_activation

          logical        :: node_activated
          integer(ikind) :: i,j
          integer        :: k


          if(.not.no_activation) then

             node_activated = .false.
   
             !there are two loops of grid-points that should be analyzed
             !and trigger the activation of the bc_interior_pt
             do k=1,2
                do j=analyzed_grdpts_bounds(2,1,k),analyzed_grdpts_bounds(2,2,k)
                   do i=analyzed_grdpts_bounds(1,1,k),analyzed_grdpts_bounds(1,2,k)
   
                      node_activated = this%is_node_activated(
     $                     [i,j],
     $                     bf_sublayer_ptr,
     $                     match_table,
     $                     interior_x_map,
     $                     interior_y_map,
     $                     interior_nodes,
     $                     p_model,
     $                     interior_domain)
   
                      if(node_activated) then
                         exit
                      end if
   
                   end do
   
                   if(node_activated) then
                      exit
                   end if
   
                end do
   
                if(node_activated) then
                   exit
                end if
   
             end do
   
             !if only one node is activated in the grid-points that are
             !analyzed, the activated bc_interior_pt are staged
             if(node_activated) then

                do k=1,nb_activated_grdpts
                   
                   call this%stage(
     $                  icr_interface_used,
     $                  activated_grdpts(:,k),
     $                  match_table,
     $                  interior_domain)
                   
                end do
                
             end if

          end if

        end subroutine analyze_bc_section_xy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of an E_edge or W_edge bc_section
        !> and stage the activated bc_interior_pt for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_bc_section_edge_x(
     $     this,
     $     icr_interface_used,
     $     bf_sublayer_ptr,
     $     match_table,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     interior_domain,
     $     bc_section)

          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr
          integer(ikind), dimension(2)    , intent(in)    :: match_table
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                         , intent(in)    :: interior_domain
          integer(ikind), dimension(5)    , intent(in)    :: bc_section

          logical       , dimension(2) :: compute_edge
          integer                      :: i_analyzed
          integer                      :: dir_analyzed
          logical                      :: analyze
          integer(ikind)               :: loc_i_analyzed
          integer(ikind), dimension(2) :: loc_central_coords
          logical                      :: node_activated
          integer(ikind)               :: last_j_added
          integer(ikind)               :: j

          logical :: update_edge
          

          call determine_edge_grdpts_computed(
     $         bc_section(5),
     $         compute_edge)
          

          !determine :
          ! - i_analyzed: the position of the bc_interior_pt on the edge
          ! - dir_analyzed: left or right direction
          ! - analyze : whether the integration overlap the bc_interior_pt
          select case(bc_section(1))

            case(E_edge_type)

               i_analyzed   =  bc_section(2)
               dir_analyzed = -1

               analyze = compute_edge(1)

            case(W_edge_type)

               i_analyzed   =  bc_section(2)+1
               dir_analyzed = +1

               analyze = compute_edge(2)

            case default
               call error_bc_section_type(
     $              'icr_interface_bc_section_module',
     $              'analyze_interior_edge_x',
     $              bc_section(1))

          end select


          ! if the time integration borders do not overlap the
          ! bc_interior_pt in the bc_section, the grid-points
          ! are analyzed
          if(analyze) then

             ! check whether the entire edge is activated to
             ! prevent a crenel to the inside of the computational
             ! domain
             update_edge = update_entire_edge(this,bf_sublayer_ptr,bc_section)

             if(update_edge) then

                do j=bc_section(3), bc_section(4)

                   call this%stage(
     $                  icr_interface_used,
     $                  [i_analyzed,j],
     $                  match_table,
     $                  interior_domain)

                end do

             ! otherwise, analyze the bc_interior_pt individually
             else

                loc_i_analyzed =
     $               i_analyzed +
     $               dir_analyzed*dct_icr_distance

                ! loop over the bc_interior_pt in the x-direction
                j=bc_section(3)
                last_j_added=j-3
                do j=bc_section(3), bc_section(4)
                
                   loc_central_coords = [loc_i_analyzed,j]
                
                   !check whether the node from which the grid-point
                   !depends is activated
                   node_activated = this%is_node_activated(
     $                  loc_central_coords,
     $                  bf_sublayer_ptr,
     $                  match_table,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes,
     $                  p_model,
     $                  interior_domain)
                
                   ! if the node is activated, the current grid-point
                   ! and its nearest neighbors are staged
                   if(node_activated) then
                      
                      if(j.ge.(last_j_added+3)) then
                      
                        ! grid-point [i,j-1] is staged
                         if(j.gt.bc_section(3)) then
                
                            call this%stage(
     $                           icr_interface_used,
     $                           [i_analyzed,j-1],
     $                           match_table,
     $                           interior_domain)
                            
                         end if
                      end if
                
                      if(j.ge.(last_j_added+2)) then
                
                         ! grid-point [i,j] is staged
                         call this%stage(
     $                           icr_interface_used,
     $                           [i_analyzed,j],
     $                           match_table,
     $                           interior_domain)
                
                      end if
                
                      last_j_added = j
                
                      ! grid-point [i,j+1] is staged
                      if(j.lt.bc_section(4)) then
                
                         call this%stage(
     $                        icr_interface_used,
     $                        [i_analyzed,j+1],
     $                        match_table,
     $                        interior_domain)
                
                      end if
                         
                   end if
                
                end do

             end if

          end if

        end subroutine analyze_bc_section_edge_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of an N_edge or S_edge bc_section
        !> and stage the activated bc_interior_pt for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_bc_section_edge_y(
     $     this,
     $     icr_interface_used,
     $     bf_sublayer_ptr,
     $     match_table,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     interior_domain,
     $     bc_section)
        
          implicit none

          class(bf_interface_icr)         , intent(in)    :: this
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_sublayer), pointer      , intent(in)    :: bf_sublayer_ptr
          integer(ikind), dimension(2)    , intent(in)    :: match_table
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                         , intent(in)    :: interior_domain
          integer(ikind), dimension(5)    , intent(in)    :: bc_section

          logical       , dimension(2) :: compute_edge
          integer                      :: j_analyzed
          integer                      :: dir_analyzed
          logical                      :: analyze
          integer(ikind)               :: loc_j_analyzed
          integer(ikind), dimension(2) :: loc_central_coords
          logical                      :: node_activated
          integer(ikind)               :: last_i_added
          integer(ikind)               :: i

          logical :: update_edge
          
          
          call determine_edge_grdpts_computed(
     $         bc_section(5),
     $         compute_edge)


          !determine :
          ! - j_analyzed: the position of the bc_interior_pt on the edge
          ! - dir_analyzed: top or bottom direction
          ! - analyze : whether the integration overlap the bc_interior_pt
          select case(bc_section(1))

            case(N_edge_type)

               j_analyzed   =  bc_section(3)
               dir_analyzed = -1
               analyze      = compute_edge(1)

            case(S_edge_type)

               j_analyzed   =  bc_section(3)+1
               dir_analyzed = +1
               analyze      = compute_edge(2)

            case default
               call error_bc_section_type(
     $              'icr_interface_bc_section_module',
     $              'analyze_interior_edge_y',
     $              bc_section(1))

          end select


          ! if the time integration borders do not overlap the
          ! bc_interior_pt in the bc_section, the grid-points
          ! are analyzed
          if(analyze) then

             ! check whether the entire edge is activated to
             ! prevent a crenel to the inside of the computational
             ! domain
             update_edge = update_entire_edge(this,bf_sublayer_ptr,bc_section)

             if(update_edge) then

                do i=bc_section(2), bc_section(4)

                   call this%stage(
     $                  icr_interface_used,
     $                  [i,j_analyzed],
     $                  match_table,
     $                  interior_domain)

                end do

             ! otherwise, analyze the bc_interior_pt individually
             else             

                loc_j_analyzed =
     $               j_analyzed +
     $               dir_analyzed*dct_icr_distance
                
                ! loop over the bc_interior_pt in the x-direction
                i=bc_section(2)
                last_i_added=i-3
                do i=bc_section(2),bc_section(4)
                
                   loc_central_coords = [i,loc_j_analyzed]
                
                   !check whether the node from which the grid-point
                   !depends is activated
                   node_activated = this%is_node_activated(
     $                  loc_central_coords,
     $                  bf_sublayer_ptr,
     $                  match_table,
     $                  interior_x_map,
     $                  interior_y_map,
     $                  interior_nodes,
     $                  p_model,
     $                  interior_domain)
                
                   ! if the node is activated, the current grid-point
                   ! and its nearest neighbors are staged
                   if(node_activated) then
                      
                      if(i.ge.(last_i_added+3)) then
                      
                        ! grid-point [i-1,j] is staged
                         if(i.gt.bc_section(2)) then
                
                            call this%stage(
     $                           icr_interface_used,
     $                           [i-1,j_analyzed],
     $                           match_table,
     $                           interior_domain)
                            
                         end if
                      end if
                
                      if(i.ge.(last_i_added+2)) then
                
                        ! grid-point [i,j] is staged
                         call this%stage(
     $                        icr_interface_used,
     $                        [i,j_analyzed],
     $                        match_table,
     $                        interior_domain)
                
                      end if
                
                      last_i_added = i
                
                      ! grid-point [i+1,j] is staged
                      if(i.lt.bc_section(4)) then
                
                         call this%stage(
     $                        icr_interface_used,
     $                        [i+1,j_analyzed],
     $                        match_table,
     $                        interior_domain)
                
                      end if
                         
                   end if
                
                end do

             end if

          end if

        end subroutine analyze_bc_section_edge_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> stage a grid-point for the update of the domain extension
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param loc_coords
        !> coordinates of the grid-point staged
        !
        !>@param match_table
        !> in case the grid-point is staged from a buffer layer, the
        !> local coordinaates should be turned into general coordinates
        !
        !>@param interior_domain
        !> indicates whether this is grid-point from the interior or
        !> from a buffer layer
        !--------------------------------------------------------------
        subroutine stage(
     $     this,
     $     icr_interface_used,
     $     local_coords,
     $     match_table,
     $     interior_domain)

          implicit none

          class(bf_interface_icr)     , intent(in)    :: this
          type(icr_interface)         , intent(inout) :: icr_interface_used
          integer(ikind), dimension(2), intent(in)    :: local_coords
          integer(ikind), dimension(2), intent(in)    :: match_table
          logical                     , intent(in)    :: interior_domain

          
          if(interior_domain) then

             call icr_interface_used%stage(
     $            local_coords,
     $            this)

          else

             call icr_interface_used%stage(
     $            [local_coords(1)+match_table(1),
     $             local_coords(2)+match_table(2)],
     $            this)

          end if

        end subroutine stage

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether a grid-point belonging to a
        !> buffer layer is activated
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param loc_coords
        !> coordinates of the grid-point investigated in the
        !> reference frame of the buffer layer
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param interior_domain
        !> indicates whether the grid-point belongs to an interior
        !> boundary section or a buffer layer
        !
        !>@return node_activated
        !> determines whether the grid-point is activated
        !--------------------------------------------------------------
        function is_node_activated(
     $       this,
     $       loc_coords,
     $       bf_sublayer_ptr,
     $       match_table,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       interior_domain)
     $       result(node_activated)

          implicit none

          class(bf_interface_icr)           , intent(in) :: this
          integer(ikind), dimension(2)      , intent(in) :: loc_coords
          type(bf_sublayer), pointer        , intent(in) :: bf_sublayer_ptr
          integer(ikind), dimension(2)      , intent(in) :: match_table
          real(rkind)  , dimension(nx)      , intent(in) :: interior_x_map
          real(rkind)  , dimension(ny)      , intent(in) :: interior_y_map
          real(rkind)  , dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                   , intent(in) :: p_model
          logical                           , intent(in) :: interior_domain
          logical                                        :: node_activated


          if(interior_domain) then

             node_activated = is_interior_node_activated(
     $            loc_coords,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model)
             
          else

             node_activated = is_bf_layer_node_activated(
     $            this,
     $            loc_coords,
     $            bf_sublayer_ptr,
     $            match_table,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model)

          end if


        end function is_node_activated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine whether a grid-point belonging to a
        !> buffer layer is activated
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param loc_coords
        !> coordinates of the grid-point investigated in the
        !> reference frame of the buffer layer
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer to which the bc_section
        !> belongs
        !
        !>@param match_table
        !> array matching the local coodinates of the buffer layer
        !> into coordinates expressed in the general reference
        !> frame
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@return node_activated
        !> determines whether the grid-point is activated
        !--------------------------------------------------------------
        function is_bf_layer_node_activated(
     $       this,
     $       loc_coords,
     $       bf_sublayer_ptr,
     $       match_table,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model)
     $       result(node_activated)

          implicit none

          class(bf_interface_icr)           , intent(in) :: this
          integer(ikind), dimension(2)      , intent(in) :: loc_coords
          type(bf_sublayer), pointer        , intent(in) :: bf_sublayer_ptr
          integer(ikind), dimension(2)      , intent(in) :: match_table
          real(rkind)  , dimension(nx)      , intent(in) :: interior_x_map
          real(rkind)  , dimension(ny)      , intent(in) :: interior_y_map
          real(rkind)  , dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                   , intent(in) :: p_model
          logical                                        :: node_activated


          logical                      :: ierror
          integer                      :: mainlayer_id
          integer(ikind), dimension(2) :: gen_coords
          type(bf_sublayer), pointer   :: bf_sublayer2_ptr
          integer(ikind), dimension(2) :: match_table2


          !> ask the buffer layer to determine whether the 
          !> grid-point is activated
          node_activated = bf_sublayer_ptr%is_node_activated(
     $         loc_coords,
     $         p_model,
     $         ierror)


          !> it is possible that the buffer layer did not
          !> have all the grid-points needed to decide whether
          !> the node was activated or not
          !> in this case, we determine to which buffer layer or
          !> interior domain, the grid-point is belonging to
          if(ierror.neqv.BF_SUCCESS) then
             
             gen_coords(1) = loc_coords(1) + match_table(1)
             gen_coords(2) = loc_coords(2) + match_table(2)             

             mainlayer_id = get_mainlayer_coord(gen_coords)

             !> if the node and its neighbors are inside the
             !> interior domain, we can simply compute whether
             !> the node is activated
             if(mainlayer_id.eq.interior) then

                node_activated = is_interior_node_activated(
     $               gen_coords,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model)
                
             !> if the node belongs to a main buffer layer, we
             !> need to find the corresponding buffer layer
             else

                bf_sublayer2_ptr => this%get_bf_layer_from_gen_coords(
     $               gen_coords,
     $               mainlayer_id=mainlayer_id)

                if(associated(bf_sublayer2_ptr)) then

                   match_table2 = bf_sublayer2_ptr%get_general_to_local_coord_tab()

                   node_activated = bf_sublayer2_ptr%is_node_activated(
     $                  [gen_coords(1)-match_table2(1),
     $                   gen_coords(2)-match_table2(2)],
     $                  p_model,
     $                  ierror)

                   !remark: even if the node and its nearest grid-points are
                   !not inside the buffer layer, the node is set to desactivated
                   !by default

                else

                   node_activated = .false.

                end if
                
             end if

          end if

        end function is_bf_layer_node_activated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the interior grid-point leads to undermined
        !> open boundary conditions
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param loc_coords
        !> general coordinates of the grid-point checked
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !> @return node_activated
        !> whether the grid-point is activated
        !--------------------------------------------------------------
        function is_interior_node_activated(
     $     loc_coords,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model)
     $     result(node_activated)

          implicit none

          integer(ikind) , dimension(2)       , intent(in) :: loc_coords
          real(rkind)    , dimension(nx)      , intent(in) :: interior_x_map
          real(rkind)    , dimension(ny)      , intent(in) :: interior_y_map
          real(rkind)    , dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                     , intent(in) :: p_model
          logical                                          :: node_activated

          integer(ikind) :: i_min
          integer(ikind) :: i_max
          integer(ikind) :: j_min
          integer(ikind) :: j_max

          i_min = loc_coords(1)-1
          i_max = loc_coords(1)+1
          j_min = loc_coords(2)-1
          j_max = loc_coords(2)+1

          node_activated = p_model%are_openbc_undermined(
     $         interior_x_map(i_min:i_max),
     $         interior_y_map(j_min:j_max),
     $         interior_nodes(i_min:i_max,j_min:j_max,:))

        end function is_interior_node_activated


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether a crenel to the inside of the computational
        !> should be removed by activating all the bc_interior_pt of
        !> the edge-like bc_section
        !
        !> @date
        !> 23_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface_grdpts_id_update augmented with procedures
        !> detecting how the domain extension should be increased
        !
        !>@param bf_sublayer_ptr
        !> pointer to the buffer layer where the edge is checked
        !
        !>@param bc_section
        !> edge-like bc_section
        !
        !> @return update_entire_edge
        !> check whether the entire edge should be activated to remove
        !> a crenel to the inside of the computational domain
        !--------------------------------------------------------------
        function update_entire_edge(this,bf_sublayer_ptr,bc_section)

          implicit none

          class(bf_interface_icr)     , intent(in) :: this
          type(bf_sublayer), pointer  , intent(in) :: bf_sublayer_ptr
          integer(ikind), dimension(5), intent(in) :: bc_section
          logical                                  :: update_entire_edge

          ! check if the removal of crenel
          ! to the inside of the computational
          ! domain is activated
          if(obc_crenel_removal_ac) then

             ! check whether there is a crenel to
             ! the inside of the computational domain
             if(associated(bf_sublayer_ptr)) then
                update_entire_edge = this%mainlayer_interfaces%analyze_bc_section_edge(
     $               bc_section,
     $               bf_sublayer_ptr)
             else
                update_entire_edge = .false.
             end if
             
          else
             update_entire_edge = .false.
          end if          

        end function update_entire_edge

      end module bf_interface_icr_class
