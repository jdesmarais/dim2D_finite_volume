      !> @file
      !> bc_operators_openbc augmented with subroutines computing
      !> the fluxes and the time derivatives at the egdes
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> bc_operators_openbc augmented with subroutines computing
      !> the fluxes and the time derivatives at the egdes
      !
      !> @date
      !> 22_10_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_openbc_normal_class
      
        use bc_operators_openbc_class, only :
     $       bc_operators_openbc

        use bf_layer_bc_fluxes_module, only :
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       extract_grdpts_to_compute_bc_fluxes        

        use bf_layer_bc_sections_overlap_module, only :
     $       determine_edge_grdpts_computed

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table

        use interface_primary, only :
     $       gradient_proc

        use parameters_constant, only :
     $       N,S,E,W,
     $       left,right

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       debug_initialize_nodes,
     $       debug_real

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only : 
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0

        implicit none

        private
        public  :: bc_operators_openbc_normal


        !> @class bc_operators_openbc_normal
        !> abstract class encapsulating subroutine to compute
        !> the fluxes at the edges of the computational domain
        !> for open boundary conditions
        !
        !>@param compute_fluxes_x_for_bc_x_edge
        !> compute the fluxes at an x-like boundary edge
        !
        !>@param compute_fluxes_y_for_bc_y_edge
        !> compute the fluxes at an y-like boundary edge
        !
        !>@param apply_bc_on_timedev_x_edge
        !> compute the time derivatives for an x-edge (E or W)
        !
        !>@param apply_bc_on_timedev_y_edge
        !> compute the time derivatives for an y-edge (N or S)
        !---------------------------------------------------------------
        type, abstract, extends(bc_operators_openbc) :: bc_operators_openbc_normal

          contains

          procedure, nopass :: compute_fluxes_y_for_bc_x_edge
          procedure, nopass :: compute_fluxes_x_for_bc_y_edge

          procedure,   pass :: apply_bc_on_timedev_N_edge
          procedure,   pass :: apply_bc_on_timedev_S_edge
          procedure,   pass :: apply_bc_on_timedev_E_edge
          procedure,   pass :: apply_bc_on_timedev_W_edge

          procedure(tdev_x_edge_l), pass, deferred :: apply_bc_on_timedev_x_edge
          procedure(tdev_y_edge_l), pass, deferred :: apply_bc_on_timedev_y_edge

        end type bc_operators_openbc_normal


        abstract interface

           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> compute the time derivatives at (i,j) resulting
           !> of the application of the boundary condition on
           !> and x edge: W_edge or E_edge
           !
           !> @date
           !> 21_10_2014 - initial version - J.L. Desmarais
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
           !>@param dx
           !> grid size along the x-axis
           !
           !>@param dy
           !> grid size along the y-axis
           !
           !>@param i
           !> grid point index along the x-axis
           !
           !>@param j
           !> grid point index along the y-axis
           !
           !>@param flux_y
           !> fluxes along the y-direction
           !
           !>@param side_x
           !> edge side to determine the boundary normal vector
           !
           !>@param gradient_x
           !> procedure to compute the gradient along the x-direction
           !> at (i,j)
           !
           !>@param timedev
           !> time derivatives of the grid points
           !--------------------------------------------------------------
           function tdev_x_edge_l(
     $       this,
     $       t,
     $       bf_x_map,
     $       bf_y_map,
     $       bf_nodes,
     $       p_model,
     $       gradient_x,
     $       i,j,
     $       flux_y,
     $       side_x)
     $       result(timedev)
           
             import bc_operators_openbc_normal
             import gradient_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_openbc_normal), intent(in) :: this
             real(rkind)                      , intent(in) :: t
             real(rkind), dimension(:)        , intent(in) :: bf_x_map
             real(rkind), dimension(:)        , intent(in) :: bf_y_map
             real(rkind), dimension(:,:,:)    , intent(in) :: bf_nodes
             type(pmodel_eq)                  , intent(in) :: p_model
             procedure(gradient_proc)                      :: gradient_x
             integer(ikind)                   , intent(in) :: i
             integer(ikind)                   , intent(in) :: j
             real(rkind), dimension(:,:,:)    , intent(in) :: flux_y
             logical                          , intent(in) :: side_x
             real(rkind), dimension(ne)                    :: timedev
           
           end function tdev_x_edge_l


           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> compute the time derivatives at (i,j) resulting
           !> of the application of the boundary condition on
           !> an y edge: N_edge or S_edge
           !
           !> @date
           !> 21_10_2014 - initial version - J.L. Desmarais
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
           !>@param dx
           !> grid size along the x-axis
           !
           !>@param dy
           !> grid size along the y-axis
           !
           !>@param i
           !> grid point index along the x-axis
           !
           !>@param j
           !> grid point index along the y-axis
           !
           !>@param flux_x
           !> fluxes along the y-direction
           !
           !>@param side_y
           !> edge side to determine the boundary normal vector
           !
           !>@param gradient_y
           !> procedure to compute the gradient along the y-direction
           !> at (i,j)
           !
           !>@param timedev
           !> time derivatives of the grid points
           !--------------------------------------------------------------
           function tdev_y_edge_l(
     $        this,
     $        t,
     $        bf_x_map,
     $        bf_y_map,
     $        bf_nodes,
     $        p_model,
     $        gradient_y,
     $        i,j,
     $        flux_x,
     $        side_y)
     $        result(timedev)
           
             import bc_operators_openbc_normal
             import gradient_proc
             import ne
             import pmodel_eq
             import ikind
             import rkind
           
             class(bc_operators_openbc_normal), intent(in) :: this
             real(rkind)                      , intent(in) :: t
             real(rkind), dimension(:)        , intent(in) :: bf_x_map
             real(rkind), dimension(:)        , intent(in) :: bf_y_map
             real(rkind), dimension(:,:,:)    , intent(in) :: bf_nodes
             type(pmodel_eq)                  , intent(in) :: p_model
             procedure(gradient_proc)                      :: gradient_y
             integer(ikind)                   , intent(in) :: i
             integer(ikind)                   , intent(in) :: j
             real(rkind), dimension(:,:,:)    , intent(in) :: flux_x
             logical                          , intent(in) :: side_y
             real(rkind), dimension(ne)                    :: timedev
           
           end function tdev_y_edge_l
      
        end interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the x-direction so that
        !> the time derivatives for an edge in the y-direction
        !> can be computed
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param compute_edge
        !> determine which grid points are computed
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_y_for_bc_x_edge(
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_op_edge1,
     $     sd_op_edge2,
     $     p_model,
     $     i,j_min, j_max,
     $     compute_edge,
     $     flux_y)
        
          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_op_edge1
          class(sd_operators)                , intent(in)    :: sd_op_edge2
          type(pmodel_eq)                    , intent(in)    :: p_model          
          integer(ikind)                     , intent(in)    :: i
          integer(ikind)                     , intent(in)    :: j_min
          integer(ikind)                     , intent(in)    :: j_max
          logical       , dimension(2)       , intent(in)    :: compute_edge
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y

          integer(ikind) :: j
          integer(ikind) :: j_lim_min
          integer(ikind) :: j_lim_max
          integer(ikind) :: j_mid_min
          integer(ikind) :: j_mid_max
          integer(ikind) :: size_y


          size_y    = size(bf_nodes,2)
          j_lim_min = bc_size
          j_lim_max = size_y-bc_size+2
          j_mid_min = max(j_min,j_lim_min+1)
          j_mid_max = min(j_max,j_lim_max-1)


          ! for the computation of the fluxes requiring to gather
          ! grid-points from the interior domain (left side)
          !============================================================
          if(j_min.le.j_lim_min) then
             
             if(compute_edge(1)) then
                call compute_flux_y_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge1,
     $               p_model,
     $               i  ,j_min,min(j_lim_min,j_max),
     $               flux_y)
             end if

             if(compute_edge(2)) then
                call compute_flux_y_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge2,
     $               p_model,
     $               i+1,j_min,min(j_lim_min,j_max),
     $               flux_y)
             end if

          end if


          ! for the computation of the fluxes that do not require
          ! to gather grid-points from the interior domain
          !============================================================
          if(compute_edge(1)) then
             
             ! compute both edges
             if(compute_edge(2)) then

                do j=j_mid_min,j_mid_max
                   
                   flux_y(i  ,j,:) = p_model%compute_flux_y_oneside(
     $                  bf_nodes,dx,dy,
     $                  i  ,j,
     $                  sd_op_edge1)
                   
                   flux_y(i+1,j,:) = p_model%compute_flux_y_oneside(
     $                  bf_nodes,dx,dy,
     $                  i+1,j,
     $                  sd_op_edge2)
                
                end do
                
             else

                ! compute only first edge
                do j=j_mid_min,j_mid_max
                   
                   flux_y(i  ,j,:) = p_model%compute_flux_y_oneside(
     $                  bf_nodes,dx,dy,
     $                  i  ,j,
     $                  sd_op_edge1)
                
                end do

             end if

          else

             if(compute_edge(2)) then

                ! compute only second edge
                do j=j_mid_min,j_mid_max
                   
                   flux_y(i+1,j,:) = p_model%compute_flux_y_oneside(
     $                  bf_nodes,dx,dy,
     $                  i+1,j,
     $                  sd_op_edge2)
                   
                end do

             end if

          end if


          ! for the computation of the fluxes requiring to gather
          ! grid-points from the interior domain (right side)
          !============================================================
          if(j_max.ge.j_lim_max) then
             
             if(compute_edge(1)) then
                call compute_flux_y_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge1,
     $               p_model,
     $               i  ,max(j_lim_max,j_min),j_max,
     $               flux_y)
             end if

             if(compute_edge(2)) then
                call compute_flux_y_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge2,
     $               p_model,
     $               i+1,max(j_lim_max,j_min),j_max,
     $               flux_y)
             end if

          end if

        end subroutine compute_fluxes_y_for_bc_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the x-direction so that
        !> the time derivatives for an edge in the y-direction
        !> can be computed
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param compute_edge
        !> determine which grid points are computed
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_x_for_bc_y_edge(
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_op_edge1,
     $     sd_op_edge2,
     $     p_model,
     $     i_min, i_max, j,
     $     compute_edge,
     $     flux_x)
        
          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          integer       , dimension(:,:)     , intent(in)    :: bf_grdpts_id
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_op_edge1
          class(sd_operators)                , intent(in)    :: sd_op_edge2
          type(pmodel_eq)                    , intent(in)    :: p_model          
          integer(ikind)                     , intent(in)    :: i_min
          integer(ikind)                     , intent(in)    :: i_max
          integer(ikind)                     , intent(in)    :: j
          logical       , dimension(2)       , intent(in)    :: compute_edge
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x

          integer(ikind) :: i
          integer(ikind) :: i_lim_min
          integer(ikind) :: i_lim_max
          integer(ikind) :: i_mid_min
          integer(ikind) :: i_mid_max
          integer(ikind) :: size_x


          size_x    = size(bf_nodes,1)
          i_lim_min = bc_size
          i_lim_max = size_x-bc_size+2
          i_mid_min = max(i_min,i_lim_min+1)
          i_mid_max = min(i_max,i_lim_max-1)


          if(compute_edge(1)) then

             ! fluxes that can only be computed by combining
             ! the current nodes with the interior_nodes
             !============================================================
             if(i_min.le.i_lim_min) then

                call compute_flux_x_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge1,
     $               p_model,
     $               i_min,min(i_lim_min,i_max),j,
     $               flux_x)

             end if


             ! fluxes that can be computed using directly 
             ! nodes
             do i=i_mid_min, i_mid_max

                flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $               bf_nodes,dx,dy,
     $               i,j,
     $               sd_op_edge1)
                
             end do


             ! fluxes that can only be computed by combining
             ! the current nodes with the interior_nodes
             if(i_max.ge.i_lim_max) then

                call compute_flux_x_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge1,
     $               p_model,
     $               max(i_lim_max,i_min),i_max,j,
     $               flux_x)

             end if

          end if


          if(compute_edge(2)) then


             ! fluxes that can only be computed by combining
             ! the current nodes with the interior_nodes
             if(i_min.le.bc_size) then

                call compute_flux_x_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge2,
     $               p_model,
     $               i_min,bc_size,j+1,
     $               flux_x)
     $                    
             end if


             ! fluxes that can be computed using directly 
             ! nodes
             do i=max(i_min,bc_size+1), min(i_max,size_x-bc_size+1)

                flux_x(i,j+1,:) = p_model%compute_flux_x_oneside(
     $               bf_nodes,dx,dy,
     $               i,j+1,
     $               sd_op_edge2)
                
             end do


             ! fluxes that can only be computed by combining
             ! the current nodes with the interior_nodes
             if(i_max.ge.(size_x-bc_size+2)) then

                call compute_flux_x_by_combining_grdpts(
     $               bf_alignment,
     $               bf_grdpts_id,
     $               bf_nodes,
     $               interior_nodes,
     $               dx,dy,
     $               sd_op_edge2,
     $               p_model,
     $               size_x-bc_size+2,i_max,j+1,
     $               flux_x)
     $                    
             end if

          end if

        end subroutine compute_fluxes_x_for_bc_y_edge
          

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for a
        !> North edge bc_section
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param s_y_L0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (-y)-direction
        !
        !>@param s_y_L1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (-y)-direction
        !
        !>@param s_y_R1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (+y)-direction
        !
        !>@param s_y_R0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (+y)-direction
        !
        !>@param dx
        !> space step in the x-direction
        !
        !>@param dy
        !> space step in the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j_min
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param overlap_type
        !> determine which grid points are computed
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_N_edge(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     s_y_R1, s_y_R0,
     $     p_model,
     $     i_min, i_max, j_min,
     $     overlap_type,
     $     flux_x,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal)  , intent(in)    :: this
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

          logical, dimension(2) :: compute_edge
          logical               :: side_y
          integer(ikind)        :: i,j
          real(rkind)           :: dx,dy


          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          
          call determine_edge_grdpts_computed(overlap_type,compute_edge)


          if((compute_edge(1).or.compute_edge(2)).and.
     $       ((i_max-i_min+1).gt.0)) then

             !compute the fluxes at the edges
             call compute_fluxes_x_for_bc_y_edge(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            dx,dy,
     $            s_y_R1, s_y_R0,
     $            p_model,
     $            i_min, i_max+1, j_min,
     $            compute_edge,
     $            flux_x)

             !compute the time derivatives
             side_y = right
          
             !1st section: j=j_min
             if(compute_edge(1)) then

                j=j_min

                do i=i_min,i_max
                
                   timedev(i,j,:) = this%apply_bc_on_timedev_y_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_y_y_oneside_R0,
     $                  i,j,
     $                  flux_x,
     $                  side_y)
                
                end do
             end if
             
             !2nd section: j=j_min+1
             if(compute_edge(2)) then
                j=j_min+1
                do i=i_min,i_max
                   
                   timedev(i,j,:) = this%apply_bc_on_timedev_y_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_y_y_oneside_R0,
     $                  i,j,
     $                  flux_x,
     $                  side_y)
                   
                end do
             end if

          end if

        end subroutine apply_bc_on_timedev_N_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for a
        !> South edge bc_section
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param s_y_L0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (-y)-direction
        !
        !>@param s_y_L1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (-y)-direction
        !
        !>@param s_y_R1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (+y)-direction
        !
        !>@param s_y_R0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (+y)-direction
        !
        !>@param dx
        !> space step in the x-direction
        !
        !>@param dy
        !> space step in the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j_min
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param overlap_type
        !> determine which grid points are computed
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_S_edge(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     s_y_L0, s_y_L1,
     $     p_model,
     $     i_min, i_max, j_min,
     $     overlap_type,
     $     flux_x,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal)  , intent(in)    :: this
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

          logical, dimension(2) :: compute_edge
          logical               :: side_y
          integer(ikind)        :: i,j
          real(rkind)           :: dx,dy


          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          call determine_edge_grdpts_computed(overlap_type,compute_edge)


          if((compute_edge(1).or.compute_edge(2)).and.
     $       ((i_max-i_min+1).gt.0)) then

             !compute the fluxes at the edges
             call compute_fluxes_x_for_bc_y_edge(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            dx,dy,
     $            s_y_L0, s_y_L1,
     $            p_model,
     $            i_min, i_max+1, j_min,
     $            compute_edge,
     $            flux_x)


             !compute the time derivatives
             side_y = left

             !1st section: j=j_min
             if(compute_edge(1)) then
                j=j_min
                do i=i_min,i_max
                   
                   timedev(i,j,:) = this%apply_bc_on_timedev_y_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_y_y_oneside_L0,
     $                  i,j,
     $                  flux_x,
     $                  side_y)

                end do
             end if

             !2nd section: j=j_min+1
             if(compute_edge(2)) then
                j=j_min+1
                do i=i_min,i_max

                   timedev(i,j,:) = this%apply_bc_on_timedev_y_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_y_y_oneside_L0,
     $                  i,j,
     $                  flux_x,
     $                  side_y)

                end do
             end if

          end if

        end subroutine apply_bc_on_timedev_S_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for a
        !> East edge bc_section
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param s_x_L0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (-x)-direction
        !
        !>@param s_x_L1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (-x)-direction
        !
        !>@param s_x_R1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (+x)-direction
        !
        !>@param s_x_R0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (+x)-direction
        !
        !>@param dx
        !> space step in the x-direction
        !
        !>@param dy
        !> space step in the y-direction
        !
        !>@param j_min
        !> index min along the y-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param j_max
        !> index max along the y-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param i_min
        !> index along the x-direction positioning the
        !> the edge boundary section
        !
        !>@param overlap_type
        !> determine which grid points are computed
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_E_edge(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     s_x_R1, s_x_R0,
     $     p_model,
     $     i_min, j_min, j_max,
     $     overlap_type,
     $     flux_y,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal)  , intent(in)    :: this
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

          logical, dimension(2) :: compute_edge
          logical               :: side_x
          integer(ikind)        :: i,j
          real(rkind)           :: dx,dy


          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          call determine_edge_grdpts_computed(overlap_type,compute_edge)


          if((compute_edge(1).or.compute_edge(2)).and.
     $       ((j_max-j_min+1).gt.0)) then

             !compute the fluxes at the edges
             call compute_fluxes_y_for_bc_x_edge(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            dx,dy,
     $            s_x_R1, s_x_R0,
     $            p_model,
     $            i_min, j_min, j_max+1,
     $            compute_edge,
     $            flux_y)

             !compute the time derivatives
             side_x = right
          
             
             ! compute the time derivatives
             !============================================================
             if(compute_edge(1)) then

                !both sides: i_min and i_min+1
                if(compute_edge(2)) then
                   do j=j_min,j_max
             
                      i=i_min
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_R0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                      i=i_min+1
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_R0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                   end do

                else

                   !one side: i_min
                   do j=j_min,j_max
             
                      i=i_min
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_R0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                   end do

                end if

             else

                !one side: i_min+1
                do j=j_min,j_max
             
                   i=i_min+1
                   timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_x_x_oneside_R0,
     $                  i,j,
     $                  flux_y,
     $                  side_x)
                   
                end do

             end if

          end if

        end subroutine apply_bc_on_timedev_E_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the time derivatives for a
        !> East edge bc_section
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param s_x_L0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (-x)-direction
        !
        !>@param s_x_L1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (-x)-direction
        !
        !>@param s_x_R1
        !> space operators needed to compute the fluxes
        !> when only one grid point is available in the (+x)-direction
        !
        !>@param s_x_R0
        !> space operators needed to compute the fluxes
        !> when no grid point is available in the (+x)-direction
        !
        !>@param dx
        !> space step in the x-direction
        !
        !>@param dy
        !> space step in the y-direction
        !
        !>@param j_min
        !> index min along the y-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param j_max
        !> index max along the y-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param i_min
        !> index along the x-direction positioning the
        !> the edge boundary section
        !
        !>@param overlap_type
        !> determine which grid points are computed
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_W_edge(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     s_x_L0, s_x_L1,
     $     p_model,
     $     i_min, j_min, j_max,
     $     overlap_type,
     $     flux_y,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal)  , intent(in)    :: this
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

          logical, dimension(2) :: compute_edge
          logical               :: side_x
          integer(ikind)        :: i,j
          real(rkind)           :: dx,dy


          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          call determine_edge_grdpts_computed(overlap_type,compute_edge)


          if((compute_edge(1).or.compute_edge(2)).and.
     $       ((j_max-j_min+1).gt.0)) then

             !compute the fluxes at the edges
             call compute_fluxes_y_for_bc_x_edge(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            dx,dy,
     $            s_x_L0, s_x_L1,
     $            p_model,
     $            i_min, j_min, j_max+1,
     $            compute_edge,
     $            flux_y)

             !compute the time derivatives
             side_x = left
          
             
             ! compute the time derivatives
             !============================================================
             if(compute_edge(1)) then

                !both sides: i_min and i_min+1
                if(compute_edge(2)) then
                   do j=j_min,j_max
             
                      i=i_min
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_L0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                      i=i_min+1
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_L0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                   end do

                else

                   !one side: i_min
                   do j=j_min,j_max
             
                      i=i_min
                      timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                     t,
     $                     bf_x_map,
     $                     bf_y_map,
     $                     bf_nodes,
     $                     p_model,
     $                     gradient_x_x_oneside_L0,
     $                     i,j,
     $                     flux_y,
     $                     side_x)
                   
                   end do

                end if

             else

                !one side: i_min+1
                do j=j_min,j_max
             
                   i=i_min+1
                   timedev(i,j,:) = this%apply_bc_on_timedev_x_edge(
     $                  t,
     $                  bf_x_map,
     $                  bf_y_map,
     $                  bf_nodes,
     $                  p_model,
     $                  gradient_x_x_oneside_L0,
     $                  i,j,
     $                  flux_y,
     $                  side_x)
                   
                end do

             end if

          end if

        end subroutine apply_bc_on_timedev_W_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the x-fluxes at the edge of the buffer layer
        !> by combining grid-points from the interior_domain
        !> and the buffer layer
        !
        !> @date
        !> 31_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_grdpts_id
        !> configuration of the grid-points in the buffer layer
        !
        !>@param bf_nodes
        !> nodes for the buffer layer
        !
        !>@param interior_nodes
        !> nodes in the interior domain
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model
        !
        !>@param i_min
        !> lower border of the flux_x computed along the x-direction
        !
        !>@param i_max
        !> upper border of the flux_x computed along the x-direction
        !
        !>@param j
        !> y-coordinate where the flux_x are computed
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_flux_x_by_combining_grdpts(
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i_min,i_max,j,
     $     flux_x)

          integer(ikind)     , dimension(2,2)      , intent(in)    :: bf_alignment
          integer            , dimension(:,:)      , intent(in)    :: bf_grdpts_id
          real(rkind)        , dimension(:,:,:)    , intent(in)    :: bf_nodes
          real(rkind)        , dimension(nx,ny,ne) , intent(in)    :: interior_nodes
          real(rkind)                              , intent(in)    :: dx
          real(rkind)                              , intent(in)    :: dy
          class(sd_operators)                      , intent(in)    :: sd_used
          type(pmodel_eq)                          , intent(in)    :: p_model
          integer(ikind)                           , intent(in)    :: i_min
          integer(ikind)                           , intent(in)    :: i_max
          integer(ikind)                           , intent(in)    :: j
          real(rkind)        , dimension(:,:,:)    , intent(inout) :: flux_x

          
          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind)   , dimension(:,:,:), allocatable :: tmp_nodes
          integer(ikind), dimension(2)                  :: match_table
          integer(ikind), dimension(2,2)                :: gen_coords

          integer(ikind) :: i_s,j_s
          integer        :: k_s

          integer(ikind) :: i


          ! determine the coordinates of the grid-points to be
          ! extracted to be able to compute the fluxes
          !============================================================
          !border_coords: general coordinates identifying the SW and NE
          !               corners of the grid-points extracted
          !
          !cpt_coords   : coordinates in the temporary array of the first
          !               grid-point computed for the fluxes
          !
          !REMARK: we use i_min,i_max"-1" b/c the grdpts_needed_for_flux_x
          !        compute the grdpts needed to compute the fluxes at i
          !        for the computation of the time derivatives, thus
          !        for the fluxes at i-1/2 and i+1/2 while here only the
          !        fluxes at i-1/2 are needed
          !============================================================
          grdpts_needed = are_grdpts_needed_for_flux_x(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i_min,i_max-1,j,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)
          

          ! extract the grid-points corresponding to border_coords
          !============================================================
          !grdpts_needed: identify whether the grid-points are needed to
          !               compute the fluxes, as i_min.le.bc_size or
          !               i_max.ge.(size_x-bc_size+1), it should always
          !               be the case
          !
          !tmp_nodes    : temporary array gathering the grid-points needed
          !               to compute the fluxes
          !
          !match_table  : array needed to turn coordinates expressed in the
          !               local frame of the buffer layer into coordinates
          !               in the general frame (interior_domain)
          !============================================================
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! DEBUG: to create instabilities if the grid-points used
             ! to compute the fluxes are not correctlty initialized
             if(debug_initialize_nodes) then
                tmp_nodes = reshape((/
     $               (((debug_real,
     $               i_s=1,size(tmp_nodes,1)),
     $               j_s=1,size(tmp_nodes,2)),
     $               k_s=1,ne)/),
     $               (/size(tmp_nodes,1),size(tmp_nodes,2),ne/))
             end if

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)


             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call extract_grdpts_to_compute_bc_fluxes(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the x-fluxes
             do i=i_min, i_max

                flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $               tmp_nodes,dx,dy,
     $               cpt_coords(1)+(i-i_min),cpt_coords(2),
     $               sd_used)

             end do

             deallocate(tmp_nodes)

          else

             print '(''bc_operators_openbc_normal'')'
             print '(''compute_flux_x_by_combining_grdpts'')'
             print '(''grdpts_needed=.false.'')'
             print '(''this is a problem since:'')'
             print '(''i_min.le.bc_size or'')'
             print '(''i_max.ge.size_x'')'
             print '(''i_min: '',I2)', i_min
             print '(''i_max: '',I2)', i_max
             stop ''

          end if

        end subroutine compute_flux_x_by_combining_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the y-fluxes at the edge of the buffer layer
        !> by combining grid-points from the interior_domain
        !> and the buffer layer
        !
        !> @date
        !> 31_03_2015 - initial version - J.L. Desmarais
        !
        !>@param bf_alignment
        !> relative position of the buffer layer compared to the
        !> interior domain
        !
        !>@param bf_grdpts_id
        !> configuration of the grid-points in the buffer layer
        !
        !>@param bf_nodes
        !> nodes for the buffer layer
        !
        !>@param interior_nodes
        !> nodes in the interior domain
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model        
        !
        !>@param i
        !> x-coordinate where the flux_y are computedx
        !
        !>@param j_min
        !> lower border of the flux_y computed along the y-direction
        !
        !>@param j_max
        !> upper border of the flux_y computed along the y-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine compute_flux_y_by_combining_grdpts(
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i,j_min,j_max,
     $     flux_y)

          integer(ikind)     , dimension(2,2)      , intent(in)    :: bf_alignment
          integer            , dimension(:,:)      , intent(in)    :: bf_grdpts_id
          real(rkind)        , dimension(:,:,:)    , intent(in)    :: bf_nodes
          real(rkind)        , dimension(nx,ny,ne) , intent(in)    :: interior_nodes
          real(rkind)                              , intent(in)    :: dx
          real(rkind)                              , intent(in)    :: dy
          class(sd_operators)                      , intent(in)    :: sd_used
          type(pmodel_eq)                          , intent(in)    :: p_model
          integer(ikind)                           , intent(in)    :: i
          integer(ikind)                           , intent(in)    :: j_min
          integer(ikind)                           , intent(in)    :: j_max
          real(rkind)        , dimension(:,:,:)    , intent(inout) :: flux_y

          
          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind)   , dimension(:,:,:), allocatable :: tmp_nodes
          integer(ikind), dimension(2)                  :: match_table
          integer(ikind), dimension(2,2)                :: gen_coords

          integer(ikind) :: i_s,j_s
          integer        :: k_s

          integer(ikind) :: j


          ! determine the coordinates of the grid-points to be
          ! extracted to be able to compute the fluxes
          !============================================================
          !border_coords: general coordinates identifying the SW and NE
          !               corners of the grid-points extracted
          !
          !cpt_coords   : coordinates in the temporary array of the first
          !               grid-point computed for the fluxes
          !
          !REMARK: we use i_min,i_max"-1" b/c the grdpts_needed_for_flux_x
          !        compute the grdpts needed to compute the fluxes at i
          !        for the computation of the time derivatives, thus
          !        for the fluxes at i-1/2 and i+1/2 while here only the
          !        fluxes at i-1/2 are needed
          !============================================================
          grdpts_needed = are_grdpts_needed_for_flux_y(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i,j_min,j_max-1,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)
          

          ! extract the grid-points corresponding to border_coords
          !============================================================
          !grdpts_needed: identify whether the grid-points are needed to
          !               compute the fluxes, as i_min.le.bc_size or
          !               i_max.ge.(size_x-bc_size+1), it should always
          !               be the case
          !
          !tmp_nodes    : temporary array gathering the grid-points needed
          !               to compute the fluxes
          !
          !match_table  : array needed to turn coordinates expressed in the
          !               local frame of the buffer layer into coordinates
          !               in the general frame (interior_domain)
          !============================================================
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! DEBUG: to create instabilities if the grid-points used
             ! to compute the fluxes are not correctlty initialized
             if(debug_initialize_nodes) then
                tmp_nodes = reshape((/
     $               (((debug_real,
     $               i_s=1,size(tmp_nodes,1)),
     $               j_s=1,size(tmp_nodes,2)),
     $               k_s=1,ne)/),
     $               (/size(tmp_nodes,1),size(tmp_nodes,2),ne/))
             end if

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)
             

             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call extract_grdpts_to_compute_bc_fluxes(
     $            bf_alignment,
     $            bf_grdpts_id,
     $            bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the y-fluxes
             do j=j_min, j_max

                flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $               tmp_nodes,dx,dy,
     $               cpt_coords(1),cpt_coords(2)+(j-j_min),
     $               sd_used)

             end do

             deallocate(tmp_nodes)

          else

             print '(''bc_operators_openbc_normal'')'
             print '(''compute_flux_y_by_combining_grdpts'')'
             print '(''grdpts_needed=.false.'')'
             print '(''this is a problem since:'')'
             print '(''j_min.le.bc_size or'')'
             print '(''j_max.ge.size_y'')'
             print '(''j_min: '',I2)', j_min
             print '(''j_max: '',I2)', j_max
             stop ''

          end if

        end subroutine compute_flux_y_by_combining_grdpts

      end module bc_operators_openbc_normal_class
