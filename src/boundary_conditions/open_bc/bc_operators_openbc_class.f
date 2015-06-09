      !> @file
      !> bc_operators_abstract augmented with interfaces to compute
      !> the time derivatives at the egdes, corners and anti-corners
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> bc_operators_abstract augmented with interfaces to compute
      !> the time derivatives at the egdes, corners and anti-corners
      !
      !> @date
      !> 22_10_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_openbc_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use bf_layer_bc_checks_module, only :
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W

        use bf_layer_bc_sections_overlap_module, only :
     $       determine_corner_or_anti_corner_grdpts_computed

        use bf_layer_errors_module, only :
     $       error_bc_section_type

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table

        use interface_primary, only :
     $       gradient_proc
        
        use parameters_bf_layer, only :
     $       cptnot_type

        use parameters_constant, only :
     $       bc_flux_and_node_choice,
     $       bc_timedev_choice,
     $       bc_fluxes_choice,
     $       N,S,E,W,
     $       left,right,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SW_edge_type,
     $       SE_edge_type,
     $       NW_edge_type,
     $       NE_edge_type

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       bc_NE_type_choice,
     $       bc_NW_type_choice,
     $       bc_SE_type_choice,
     $       bc_SW_type_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $     sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $     sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $     sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $     sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $     sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $     sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $     sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $     sd_operators_y_oneside_R0

        implicit none

        private
        public :: bc_operators_openbc


        !> @class bc_operators_openbc
        !> abstract class encapsulating interfaces to compute the 
        !> time derivatives at the egdes, corners and anti-corners
        !
        !>@param check_x_flux_interactions_btw_bcs
        !> check whether the computation of the x-fluxes overlap between
        !> the N/S and E/W layers
        !
        !>@param check_y_flux_interactions_btw_bcs
        !> check whether the computation of the y-fluxes overlap between
        !> the N/S and E/W layers
        !
        !>@param apply_bc_on_timedev_nopt
        !> compute the time derivatives based on the boundary conditions
        !> using the bc_section
        !
        !>@param compute_timedev_corner
        !> compute the time derivatives for a corner
        !---------------------------------------------------------------
        type, extends(bc_operators_default), abstract :: bc_operators_openbc

           contains

           procedure                , nopass         :: check_x_flux_interactions_btw_bcs
           procedure                , nopass         :: check_y_flux_interactions_btw_bcs

           procedure                , pass           :: apply_bc_on_timedev_nopt

           procedure(tdev_N_edge)   , pass, deferred :: apply_bc_on_timedev_N_edge
           procedure(tdev_S_edge)   , pass, deferred :: apply_bc_on_timedev_S_edge
           procedure(tdev_E_edge)   , pass, deferred :: apply_bc_on_timedev_E_edge
           procedure(tdev_W_edge)   , pass, deferred :: apply_bc_on_timedev_W_edge
           procedure(tdev_xy_corner), pass, deferred :: apply_bc_on_timedev_xy_corner

           procedure                , pass           :: compute_timedev_corner
           procedure(tdev_xy_edge)  , pass, deferred :: compute_timedev_anti_corner

        end type bc_operators_openbc


        abstract interface

           subroutine tdev_N_edge(
     $       this,
     $       t,
     $       bf_alignment,
     $       bf_grdpts_id,
     $       bf_x_map,
     $       bf_y_map,
     $       bf_nodes,
     $       interior_nodes,
     $       s_y_R1, s_y_R0,
     $       p_model,
     $       i_min, i_max, j_min,
     $       overlap_type,
     $       flux_x,
     $       timedev)

            import bc_operators_openbc
            import ikind
            import nx,ny,ne
            import pmodel_eq
            import rkind
            import sd_operators_y_oneside_R1
            import sd_operators_y_oneside_R0
            
            class(bc_operators_openbc)         , intent(in)    :: this
            real(rkind)                        , intent(in)    :: t
            integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
            integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
            real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
            real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
            real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
            real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
            type(sd_operators_y_oneside_R1)    , intent(in)    :: s_y_R1
            type(sd_operators_y_oneside_R0)    , intent(in)    :: s_y_R0
            type(pmodel_eq)                    , intent(in)    :: p_model
            integer(ikind)                     , intent(in)    :: i_min
            integer(ikind)                     , intent(in)    :: i_max
            integer(ikind)                     , intent(in)    :: j_min
            integer                            , intent(in)    :: overlap_type
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
           
          end subroutine tdev_N_edge


          subroutine tdev_S_edge(
     $       this,
     $       t,
     $       bf_alignment,
     $       bf_grdpts_id,
     $       bf_x_map,
     $       bf_y_map,
     $       bf_nodes,
     $       interior_nodes,
     $       s_y_L0, s_y_L1,
     $       p_model,
     $       i_min, i_max, j_min,
     $       overlap_type,
     $       flux_x,
     $       timedev)

            import bc_operators_openbc
            import ikind
            import nx,ny,ne
            import pmodel_eq
            import rkind
            import sd_operators_y_oneside_L0
            import sd_operators_y_oneside_L1          
            
            class(bc_operators_openbc)         , intent(in)    :: this
            real(rkind)                        , intent(in)    :: t
            integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
            integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
            real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
            real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
            real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
            real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
            type(sd_operators_y_oneside_L0)    , intent(in)    :: s_y_L0
            type(sd_operators_y_oneside_L1)    , intent(in)    :: s_y_L1
            type(pmodel_eq)                    , intent(in)    :: p_model
            integer(ikind)                     , intent(in)    :: i_min
            integer(ikind)                     , intent(in)    :: i_max
            integer(ikind)                     , intent(in)    :: j_min
            integer                            , intent(in)    :: overlap_type
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev

          end subroutine tdev_S_edge


          subroutine tdev_E_edge(
     $      this,
     $      t,
     $      bf_alignment,
     $      bf_grdpts_id,
     $      bf_x_map,
     $      bf_y_map,
     $      bf_nodes,
     $      interior_nodes,
     $      s_x_R1, s_x_R0,
     $      p_model,
     $      i_min, j_min, j_max,
     $      overlap_type,
     $      flux_y,
     $      timedev)
           
            import bc_operators_openbc
            import nx,ny,ne
            import pmodel_eq
            import ikind
            import rkind
            import sd_operators_x_oneside_R1
            import sd_operators_x_oneside_R0
           
            class(bc_operators_openbc)         , intent(in)    :: this
            real(rkind)                        , intent(in)    :: t
            integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
            integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
            real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
            real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
            real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
            real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
            type(sd_operators_x_oneside_R1)    , intent(in)    :: s_x_R1
            type(sd_operators_x_oneside_R0)    , intent(in)    :: s_x_R0
            type(pmodel_eq)                    , intent(in)    :: p_model
            integer(ikind)                     , intent(in)    :: i_min
            integer(ikind)                     , intent(in)    :: j_min
            integer(ikind)                     , intent(in)    :: j_max
            integer                            , intent(in)    :: overlap_type
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
           
          end subroutine tdev_E_edge


          subroutine tdev_W_edge(
     $      this,
     $      t,
     $      bf_alignment,
     $      bf_grdpts_id,
     $      bf_x_map,
     $      bf_y_map,
     $      bf_nodes,
     $      interior_nodes,
     $      s_x_L0, s_x_L1,
     $      p_model,
     $      i_min, j_min, j_max,
     $      overlap_type,
     $      flux_y,
     $      timedev)
           
            import bc_operators_openbc
            import nx,ny,ne
            import pmodel_eq
            import ikind
            import rkind
            import sd_operators_x_oneside_L0
            import sd_operators_x_oneside_L1
           
            class(bc_operators_openbc)         , intent(in)    :: this
            real(rkind)                        , intent(in)    :: t
            integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
            integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
            real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
            real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
            real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
            real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
            type(sd_operators_x_oneside_L0)    , intent(in)    :: s_x_L0
            type(sd_operators_x_oneside_L1)    , intent(in)    :: s_x_L1
            type(pmodel_eq)                    , intent(in)    :: p_model
            integer(ikind)                     , intent(in)    :: i_min
            integer(ikind)                     , intent(in)    :: j_min
            integer(ikind)                     , intent(in)    :: j_max
            integer                            , intent(in)    :: overlap_type
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
            real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
           
         end subroutine tdev_W_edge


         function tdev_xy_corner(
     $     this,
     $     t,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     p_model,
     $     i,j,
     $     side_x, side_y)
     $     result(timedev)
         
           import bc_operators_openbc
           import ikind
           import ne
           import pmodel_eq
           import rkind
           
           class(bc_operators_openbc)   , intent(in) :: this
           real(rkind)                  , intent(in) :: t
           real(rkind), dimension(:)    , intent(in) :: bf_x_map
           real(rkind), dimension(:)    , intent(in) :: bf_y_map
           real(rkind), dimension(:,:,:), intent(in) :: bf_nodes
           type(pmodel_eq)              , intent(in) :: p_model
           integer(ikind)               , intent(in) :: i
           integer(ikind)               , intent(in) :: j
           logical                      , intent(in) :: side_x
           logical                      , intent(in) :: side_y
           real(rkind), dimension(ne)                :: timedev
           
         end function tdev_xy_corner


         subroutine tdev_xy_edge(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     s_x_L1,
     $     s_x_R1,
     $     s_y_L1,
     $     s_y_R1,
     $     p_model,
     $     bc_section,
     $     flux_x,
     $     flux_y,
     $     timedev)
           
           import bc_operators_openbc
           import pmodel_eq
           import ikind
           import rkind

           import nx,ny,ne

           import sd_operators_x_oneside_L1
           import sd_operators_x_oneside_R1
           import sd_operators_y_oneside_L1
           import sd_operators_y_oneside_R1
           
           class(bc_operators_openbc)         , intent(in)    :: this
           real(rkind)                        , intent(in)    :: t
           integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
           integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
           real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
           real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
           real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
           real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
           type(sd_operators_x_oneside_L1)    , intent(in)    :: s_x_L1
           type(sd_operators_x_oneside_R1)    , intent(in)    :: s_x_R1
           type(sd_operators_y_oneside_L1)    , intent(in)    :: s_y_L1
           type(sd_operators_y_oneside_R1)    , intent(in)    :: s_y_R1
           type(pmodel_eq)                    , intent(in)    :: p_model
           integer       , dimension(5)       , intent(in)    :: bc_section
           real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
           real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
           real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
           
         end subroutine tdev_xy_edge

        end interface

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the x-fluxes are re-computed when
        !> several boundary conditions interact at a corner
        !
        !> @date
        !> 09_06_2015 - initial version - J.L. Desmarais
        !
        !>@param i_min
        !> min x-index where the x-fluxes are computed
        !
        !>@param i_max
        !> max x-index where the x-fluxes are computed
        !-------------------------------------------------------------
        subroutine check_x_flux_interactions_btw_bcs(
     $       i_min,i_max,
     $       side,
     $       bf_alignment)

          implicit none

          integer(ikind)                          , intent(inout) :: i_min
          integer(ikind)                          , intent(inout) :: i_max
          logical                                 , intent(in)    :: side
          integer(ikind), dimension(2,2), optional, intent(in)    :: bf_alignment


          integer(ikind), dimension(2) :: match_table
          integer(ikind)               :: i_min_g
          integer(ikind)               :: i_max_g


          ! when the x-flxues are computed for an open b.c., this is
          ! either for the North or for the South layers
          ! if the boundary conditions on the East or West layers are
          ! modifying the x-fluxes, then the x-fluxes should not be
          ! recomputed by the curretn boundary conditions, therefore
          ! the min-max indices can be shifted to avoid to re-compute
          ! the fluxes
          if(present(bf_alignment)) then

             match_table = get_bf_layer_match_table(
     $            bf_alignment)

             i_min_g = i_min + match_table(1)
             i_max_g = i_max + match_table(1)


             ! South layer
             if(side.eqv.left) then

                ! if the West b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=bc_size+1
                if((i_min_g.eq.(bc_size+1)).and.(
     $             (bc_SW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SW_type_choice.eq.bc_fluxes_choice))) then
                   i_min=i_min+1
                end if

                ! if the East b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=nx-bc_size+1
                if((i_max_g.eq.(nx-bc_size+1)).and.(
     $             (bc_SE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SE_type_choice.eq.bc_fluxes_choice))) then
                   i_max=i_max-1
                end if

             ! North layer
             else

                ! if the West b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=bc_size+1
                if((i_min_g.eq.(bc_size+1)).and.(
     $             (bc_NW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NW_type_choice.eq.bc_fluxes_choice))) then
                   i_min=i_min+1
                end if

                ! if the East b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=nx-bc_size+1
                if((i_max_g.eq.(nx-bc_size+1)).and.(
     $             (bc_NE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NE_type_choice.eq.bc_fluxes_choice))) then
                   i_max=i_max-1
                end if                

             end if

          else          
          
             ! South layer
             if(side.eqv.left) then

                ! if the West b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=bc_size+1
                if((bc_SW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SW_type_choice.eq.bc_fluxes_choice)) then
                   i_min=i_min+1
                end if
             
                ! if the East b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=nx-bc_size+1
                if((bc_SE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SE_type_choice.eq.bc_fluxes_choice))then
                   i_max=i_max-1
                end if

             ! North layer
             else

                ! if the West b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=bc_size+1
                if((bc_NW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NW_type_choice.eq.bc_fluxes_choice)) then
                   i_min=i_min+1
                end if
             
                ! if the East b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at i=nx-bc_size+1
                if((bc_NE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NE_type_choice.eq.bc_fluxes_choice))then
                   i_max=i_max-1
                end if

             end if

          end if

        end subroutine check_x_flux_interactions_btw_bcs


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the y-fluxes are re-computed when
        !> several boundary conditions interact at a corner
        !
        !> @date
        !> 09_06_2015 - initial version - J.L. Desmarais
        !
        !>@param j_min
        !> min y-index where the x-fluxes are computed
        !
        !>@param j_max
        !> max y-index where the x-fluxes are computed
        !-------------------------------------------------------------
        subroutine check_y_flux_interactions_btw_bcs(
     $     j_min,j_max,
     $     side,
     $     bf_alignment)

          implicit none

          integer(ikind)                          , intent(inout) :: j_min
          integer(ikind)                          , intent(inout) :: j_max
          logical                                 , intent(in)    :: side
          integer(ikind), dimension(2,2), optional, intent(in)    :: bf_alignment


          integer(ikind), dimension(2) :: match_table
          integer(ikind)               :: j_min_g
          integer(ikind)               :: j_max_g


          ! when the y-flxues are computed for an open b.c., this is
          ! either for the West or for the East layers
          ! if the boundary conditions on the North or South layers are
          ! modifying the y-fluxes, then the y-fluxes should not be
          ! recomputed by the curretn boundary conditions, therefore
          ! the min-max indices can be shifted to avoid to re-compute
          ! the fluxes
          if(present(bf_alignment)) then

             match_table = get_bf_layer_match_table(
     $            bf_alignment)

             j_min_g = j_min + match_table(2)
             j_max_g = j_max + match_table(2)

             
             ! West layer
             if(side.eqv.left) then

                ! if the South b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=bc_size+1
                if((j_min_g.eq.(bc_size+1)).and.(
     $             (bc_SW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SW_type_choice.eq.bc_fluxes_choice))) then
                   j_min=j_min+1
                end if
                
                ! if the North b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=ny-bc_size+1
                if((j_max_g.eq.(ny-bc_size+1)).and.(
     $             (bc_NW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NW_type_choice.eq.bc_fluxes_choice))) then
                   j_max=j_max-1
                end if

             ! East layer
             else

                ! if the South b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=bc_size+1
                if((j_min_g.eq.(bc_size+1)).and.(
     $               (bc_SE_type_choice.eq.bc_flux_and_node_choice).or.
     $               (bc_SE_type_choice.eq.bc_fluxes_choice))) then
                   j_min=j_min+1
                end if
                
                ! if the North b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=ny-bc_size+1
                if((j_max_g.eq.(ny-bc_size+1)).and.(
     $               (bc_NE_type_choice.eq.bc_flux_and_node_choice).or.
     $               (bc_NE_type_choice.eq.bc_fluxes_choice))) then
                   j_max=j_max-1
                end if

             end if

          else

             ! West layer
             if(side.eqv.left) then
                          
                ! if the South b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=bc_size+1
                if((bc_SW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SW_type_choice.eq.bc_fluxes_choice)) then
                   j_min=j_min+1
                end if
               
                ! if the North b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=ny-bc_size+1
                if((bc_NW_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NW_type_choice.eq.bc_fluxes_choice)) then
                   j_max=j_max-1
                end if

             ! East layer
             else

                ! if the South b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=bc_size+1
                if((bc_SE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_SE_type_choice.eq.bc_fluxes_choice)) then
                   j_min=j_min+1
                end if
               
                ! if the North b.c. modifies the fluxes, we should not
                ! re-compute the fluxes at j=ny-bc_size+1
                if((bc_NE_type_choice.eq.bc_flux_and_node_choice).or.
     $             (bc_NE_type_choice.eq.bc_fluxes_choice)) then
                   j_max=j_max-1
                end if

             end if

          end if

        end subroutine check_y_flux_interactions_btw_bcs


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the time-derivatives for a specific
        !> boundary section on a sub-domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> section computed for the boundary conditions
        !> of the sub-domain
        !
        !>@param bf_alignment
        !> relative alignment of the sub-domain compared
        !> to the interior domain
        !
        !>@param bf_grdpts_id
        !> identification of the grid-points in the
        !> sub-domain
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param nodes
        !> governing variables
        !
        !>@param interior_nodes
        !> grid-points of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param timedev
        !> tim derivatives of the sub-domain
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_nopt(
     $     this,
     $     bc_section,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     t,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     p_model,
     $     flux_x, flux_y,
     $     timedev)
        
          implicit none

          class(bc_operators_openbc)         , intent(in)    :: this
          integer       , dimension(5)       , intent(in)    :: bc_section
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
          real(rkind)                        , intent(in)    :: t
          real(rkind)   , dimension(:)       , intent(in)    :: bf_x_map
          real(rkind)   , dimension(:)       , intent(in)    :: bf_y_map
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev

          
          !spatial discretisation operators
          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0

          !intermediate variables
          real(rkind)    :: dx,dy
          integer(ikind) :: i_min, i_max
          integer(ikind) :: j_min, j_max
          integer        :: overlap_type
          

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          overlap_type = bc_section(5)


          !identify the type of boundary layer
          select case(bc_section(1))

            case(N_edge_type)

               j_min = bc_section(3)
               i_min = bc_section(2)
               i_max = bc_section(4)

               call this%apply_bc_on_timedev_N_edge(
     $              t,
     $              bf_alignment,
     $              bf_grdpts_id,
     $              bf_x_map,
     $              bf_y_map,
     $              bf_nodes,
     $              interior_nodes,
     $              s_y_R1, s_y_R0,
     $              p_model,
     $              i_min, i_max, j_min,
     $              overlap_type,
     $              flux_x,
     $              timedev)
               

            case(S_edge_type)
             
                j_min = bc_section(3)
                i_min = bc_section(2)
                i_max = bc_section(4)
             
                call this%apply_bc_on_timedev_S_edge(
     $               t,
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_x_map,
     $               bf_y_map,
     $               bf_nodes,
     $               interior_nodes,
     $               s_y_L0, s_y_L1,
     $               p_model,
     $               i_min, i_max, j_min,
     $               overlap_type,
     $               flux_x,
     $               timedev)

             
             case(E_edge_type)
             
                i_min = bc_section(2)
                j_min = bc_section(3)
                j_max = bc_section(4)
             
                call this%apply_bc_on_timedev_E_edge(
     $               t,
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_x_map,
     $               bf_y_map,
     $               bf_nodes,
     $               interior_nodes,
     $               s_x_R1, s_x_R0,
     $               p_model,
     $               i_min, j_min, j_max,
     $               overlap_type,
     $               flux_y,
     $               timedev)

             
             case(W_edge_type)
             
                i_min = bc_section(2)
                j_min = bc_section(3)
                j_max = bc_section(4)
                
                call this%apply_bc_on_timedev_W_edge(
     $               t,
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_x_map,
     $               bf_y_map,
     $               bf_nodes,
     $               interior_nodes,
     $               s_x_L0, s_x_L1,
     $               p_model,
     $               i_min, j_min, j_max,
     $               overlap_type,
     $               flux_y,
     $               timedev)


             ! corner type bc_section
             case(NE_corner_type,
     $            NW_corner_type,
     $            SE_corner_type,
     $            SW_corner_type)

               call this%compute_timedev_corner(
     $            t,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            bc_section,
     $            timedev)


             ! anti-corner type bc_section
             case(NE_edge_type,
     $            NW_edge_type,
     $            SE_edge_type,
     $            SW_edge_type)

               call this%compute_timedev_anti_corner(
     $            t,
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            interior_nodes,
     $            s_x_L1,
     $            s_x_R1,
     $            s_y_L1,
     $            s_y_R1,
     $            p_model,
     $            bc_section,
     $            flux_x,
     $            flux_y,
     $            timedev)


             case default
                call error_bc_section_type(
     $               'bc_operators_openbc_class',
     $               'apply_bc_on_timedev_nopt',
     $               bc_section(1))

           end select
                     
        end subroutine apply_bc_on_timedev_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an xy edge: NE_edge, NW_edge, SE_edge, SW_edge
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary operator
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param bf_x_map
        !> coordinates along the x-direction on the sub-domain
        !
        !>@param bf_y_map
        !> coordinates along the y-direction on the sub-domain
        !
        !>@param bf_nodes
        !> object encapsulating the main variables on the sub-domain
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param bc_section
        !> type of corner + properties on the location
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine compute_timedev_corner(
     $     this,
     $     t,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     p_model,
     $     bc_section,
     $     timedev)
        
          implicit none
        
          class(bc_operators_openbc)     , intent(in)    :: this
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:)      , intent(in)    :: bf_x_map
          real(rkind), dimension(:)      , intent(in)    :: bf_y_map
          real(rkind), dimension(:,:,:)  , intent(in)    :: bf_nodes
          type(pmodel_eq)                , intent(in)    :: p_model
          integer    , dimension(5)      , intent(in)    :: bc_section
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev
          
          
          integer(ikind)        :: i_min, j_min
          logical               :: side_x, side_y
          integer, dimension(4) :: compute_point


          i_min=bc_section(2)
          j_min=bc_section(3)

          select case(bc_section(1))

            case(SW_corner_type)

               side_x = left
               side_y = left


            case(SE_corner_type)
               
               side_x = right
               side_y = left

                     
            case(NW_corner_type)
               
               side_x = left
               side_y = right


            case(NE_corner_type)
               
               side_x = right
               side_y = right


            case default
               call error_bc_section_type(
     $              'bc_operators_openbc_class',
     $              'compute_timedev_corner',
     $              bc_section(1))

          end select

          !computation of the corner pts
          call determine_corner_or_anti_corner_grdpts_computed(
     $         bc_section(4),
     $         bc_section(5),
     $         compute_point)
          
          call compute_timedev_corner_pts(
     $         this,
     $         t,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes,
     $         p_model,
     $         i_min, j_min,
     $         side_x, side_y,
     $         compute_point,
     $         timedev)

        end subroutine compute_timedev_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives resulting from the
        !> application of the boundary condition + choose
        !> which grid points of the corner are computed
        !
        !> @date
        !> 04_03_2015 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param x_map
        !> coordinates along the x-direction
        !
        !>@param y_map
        !> coordinates along the y-direction
        !
        !>@param side_x
        !> direction in which the waves are incoming in the
        !> x-direction
        !
        !>@param side_y
        !> direction in which the waves are incoming in the
        !> y-direction
        !
        !>@param i_min
        !> index identifying the SW border of the corner in the
        !> x-direction
        !
        !>@param j_min
        !> index identifying the SW border of the corner in the
        !> y-direction
        !
        !>@param compute_point
        !> grid-points computed in the corner
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine compute_timedev_corner_pts(
     $     this,
     $     t,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     p_model,
     $     i_min, j_min,
     $     side_x, side_y,
     $     compute_point,
     $     timedev)

          implicit none

          class(bc_operators_openbc)     , intent(in)    :: this
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:)      , intent(in)    :: bf_x_map
          real(rkind), dimension(:)      , intent(in)    :: bf_y_map
          real(rkind), dimension(:,:,:)  , intent(in)    :: bf_nodes
          type(pmodel_eq)                , intent(in)    :: p_model
          integer(ikind)                 , intent(in)    :: i_min
          integer(ikind)                 , intent(in)    :: j_min
          logical                        , intent(in)    :: side_x
          logical                        , intent(in)    :: side_y
          integer    , dimension(4)      , intent(in)    :: compute_point
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev


          integer(ikind) :: i,j


          if(compute_point(1).ne.cptnot_type) then

             i=i_min
             j=j_min

             timedev(i,j,:) = this%apply_bc_on_timedev_xy_corner(
     $            t,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            i,j,
     $            side_x, side_y)

          end if

          if(compute_point(2).ne.cptnot_type) then

             i=i_min+1
             j=j_min

             timedev(i,j,:) = this%apply_bc_on_timedev_xy_corner(
     $            t,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            i,j,
     $            side_x, side_y)

          end if

          if(compute_point(3).ne.cptnot_type) then

             i=i_min
             j=j_min+1

             timedev(i,j,:) = this%apply_bc_on_timedev_xy_corner(
     $            t,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            i,j,
     $            side_x, side_y)

          end if

          if(compute_point(4).ne.cptnot_type) then

             i=i_min+1
             j=j_min+1

             timedev(i,j,:) = this%apply_bc_on_timedev_xy_corner(
     $            t,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            p_model,
     $            i,j,
     $            side_x, side_y)

          end if

        end subroutine compute_timedev_corner_pts

      end module bc_operators_openbc_class 
