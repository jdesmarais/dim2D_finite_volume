      !> @file
      !> abstract class encapsulating subroutine to compute
      !> the fluxes at the edges of the computational domain
      !> for open boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutine to compute
      !> the fluxes at the edges of the computational domain
      !> for open boundary conditions
      !
      !> @date
      !> 22_10_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_openbc_class

        use bc_operators_default_class, only :
     $     bc_operators_default

        use bc_operators_nopt_module, only :
     $     compute_edge_N,
     $     compute_edge_S,
     $     compute_edge_E,
     $     compute_edge_W  

        use bf_layer_bc_procedure_module, only :
     $     N_edge_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     SW_corner_type,
     $     SE_corner_type,
     $     NW_corner_type,
     $     NE_corner_type,
     $     SW_edge_type,
     $     SE_edge_type,
     $     NW_edge_type,
     $     NE_edge_type

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use parameters_constant, only :
     $     bc_timedev_choice,
     $     N,S,E,W,
     $     left,right

        use parameters_input, only :
     $     nx, bc_size, ne

        use parameters_kind, only :
     $     ikind,
     $     rkind

        use pmodel_eq_class, only :
     $     pmodel_eq

        use sd_operators_fd_module, only :
     $     gradient_x_x_oneside_L0,
     $     gradient_x_x_oneside_L1,
     $     gradient_x_x_oneside_R1,
     $     gradient_x_x_oneside_R0,
     $     gradient_y_y_oneside_L0,
     $     gradient_y_y_oneside_L1,
     $     gradient_y_y_oneside_R1,
     $     gradient_y_y_oneside_R0

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
        !> abstract class encapsulating subroutine to compute
        !> the fluxes at the edges of the computational domain
        !> for open boundary conditions
        !
        !>@param compute_fluxes_for_bc_x_edge
        !> compute the fluxes at an x-like boundary edge
        !
        !>@param compute_fluxes_for_bc_y_edge
        !> compute the fluxes at an y-like boundary edge
        !---------------------------------------------------------------
        type, extends(bc_operators_default), abstract :: bc_operators_openbc

           contains

           procedure                , pass           :: apply_bc_on_timedev_nopt
           procedure                , pass           :: compute_fluxes_for_bc_x_edge
           procedure                , pass           :: compute_fluxes_for_bc_y_edge
           procedure(tdev_x_edge)   , pass, deferred :: apply_bc_on_timedev_x_edge
           procedure(tdev_y_edge)   , pass, deferred :: apply_bc_on_timedev_y_edge
           procedure(tdev_xy_corner), pass, deferred :: apply_bc_on_timedev_xy_corner

        end type bc_operators_openbc


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
           function tdev_x_edge(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_y,
     $        side_x,
     $        gradient_x)
     $        result(timedev)
           
             import bc_operators_openbc
             import gradient_x_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_openbc)   , intent(in) :: this
             type(pmodel_eq)              , intent(in) :: p_model
             real(rkind)                  , intent(in) :: t
             real(rkind), dimension(:,:,:), intent(in) :: nodes
             real(rkind), dimension(:)    , intent(in) :: x_map
             real(rkind), dimension(:)    , intent(in) :: y_map
             integer(ikind)               , intent(in) :: i
             integer(ikind)               , intent(in) :: j
             real(rkind), dimension(:,:,:), intent(in) :: flux_y
             logical                      , intent(in) :: side_x
             procedure(gradient_x_proc)                :: gradient_x
             real(rkind), dimension(ne)                :: timedev
           
           end function tdev_x_edge


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
           function tdev_y_edge(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_x, side_y, gradient_y)
     $        result(timedev)
           
             import bc_operators_openbc
             import gradient_y_proc
             import ne
             import pmodel_eq
             import ikind
             import rkind
           
             class(bc_operators_openbc)   , intent(in) :: this
             type(pmodel_eq)              , intent(in) :: p_model
             real(rkind)                  , intent(in) :: t
             real(rkind), dimension(:,:,:), intent(in) :: nodes
             real(rkind), dimension(:)    , intent(in) :: x_map
             real(rkind), dimension(:)    , intent(in) :: y_map
             integer(ikind)               , intent(in) :: i
             integer(ikind)               , intent(in) :: j
             real(rkind), dimension(:,:,:), intent(in) :: flux_x
             logical                      , intent(in) :: side_y
             procedure(gradient_y_proc)                :: gradient_y
             real(rkind), dimension(ne)                :: timedev
           
           end function tdev_y_edge
           
           
           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> compute the time derivatives at (i,j) resulting
           !> of the application of the boundary condition on
           !> a corner: SE_corner, SW_corner, NE_corner, NW_corner
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
           function tdev_xy_corner(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        side_x, side_y,
     $        gradient_x, gradient_y)
     $        result(timedev)
           
             import bc_operators_openbc
             import gradient_x_proc
             import gradient_y_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_openbc)   , intent(in) :: this
             type(pmodel_eq)              , intent(in) :: p_model
             real(rkind)                  , intent(in) :: t
             real(rkind), dimension(:,:,:), intent(in) :: nodes
             real(rkind), dimension(:)    , intent(in) :: x_map
             real(rkind), dimension(:)    , intent(in) :: y_map
             integer(ikind)               , intent(in) :: i
             integer(ikind)               , intent(in) :: j
             logical                      , intent(in) :: side_x
             logical                      , intent(in) :: side_y
             procedure(gradient_x_proc)                :: gradient_x
             procedure(gradient_y_proc)                :: gradient_y
             real(rkind), dimension(ne)                :: timedev
           
           end function tdev_xy_corner

        end interface

        contains

        !apply the boundary conditions on the time derivatives
        subroutine apply_bc_on_timedev_nopt(
     $       this,
     $       p_model,
     $       t, nodes, x_map, y_map,
     $       flux_x, flux_y,
     $       timedev,
     $       bc_sections)
        
          implicit none

          class(bc_operators_openbc)              , intent(in)    :: this
          type(pmodel_eq)                         , intent(in)    :: p_model
          real(rkind)                             , intent(in)    :: t
          real(rkind), dimension(:,:,:)           , intent(in)    :: nodes
          real(rkind), dimension(:)               , intent(in)    :: x_map
          real(rkind), dimension(:)               , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)           , intent(inout) :: flux_x
          real(rkind), dimension(:,:,:)           , intent(inout) :: flux_y
          real(rkind), dimension(:,:,:)           , intent(inout) :: timedev
          integer    , dimension(:,:), allocatable, intent(in)    :: bc_sections

          
          !spatial discretisation operators
          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0
          

          !if there are effectively boundary layers
          !in the buffer layer computed, the time
          !derivatives corresponding to the boundary
          !grid points are computed
          if(allocated(bc_sections)) then


             !compute the fluxes at the grid points
             !corresponding to edge-type boundary
             !layers
             call compute_fluxes_bc_sections(
     $            this,
     $            p_model,
     $            nodes,x_map,y_map,
     $            s_x_L0,s_x_L1,
     $            s_x_R1,s_x_R0,
     $            s_y_L0,s_y_L1,
     $            s_y_R1,s_y_R0,
     $            flux_x,flux_y,
     $            bc_sections)


             !compute the time derivatives
             !corresponding to the buffer layers
             call compute_timedev_bc_sections(
     $            this,
     $            p_model,
     $            t,nodes,x_map,y_map,
     $            flux_x, flux_y,
     $            timedev,
     $            bc_sections)

          end if          

        end subroutine apply_bc_on_timedev_nopt


        !compute the edge fluxes needed for the
        !computation of the time derivatives in
        !the boundary layers
        subroutine compute_fluxes_bc_sections(
     $     this,
     $     p_model,
     $     nodes,x_map,y_map,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     flux_x, flux_y,
     $     bc_sections)

          implicit none
          
          class(bc_operators_openbc)     , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y
          integer    , dimension(:,:)    , intent(in)    :: bc_sections

          real(rkind)    :: dx,dy
          integer(ikind) :: i_min, i_max
          integer(ikind) :: j_min, j_max
          integer(ikind) :: i,j
          integer        :: k
          logical        :: compute_edge

          
          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)


          !go through the boundary layers
          !if the boundary actually needs the computation
          !of the fluxes in the direction of the edge, the
          !fluxes are computed
          do k=1, size(bc_sections,2)

             !identify the type of boundary layer
             select case(bc_sections(1,k))

               case(N_edge_type)

                  !do not compute the edge fluxes only if
                  !(y.ge.bc_y_max).and.
                  !(bc_N_type_choice.ne.bc_timedev_choice)

                  j = bc_sections(3,k)

                  compute_edge = compute_edge_N(j,y_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     i_min = bc_sections(2,k)
                     i_max = bc_sections(4,k)+1

                     call this%compute_fluxes_for_bc_y_edge(
     $                    p_model,
     $                    nodes,
     $                    s_y_L0, s_y_L1,
     $                    s_y_R1, s_y_R0,
     $                    dx, dy,
     $                    i_min, i_max, j,
     $                    N,
     $                    flux_x)

                  end if
                     
               case(S_edge_type)

                  !do not compute the edge fluxes only if
                  !(y.le.bc_y_min).and.
                  !(bc_S_type_choice.ne.bc_timedev_choice)

                  j = bc_sections(3,k)

                  compute_edge = compute_edge_S(j,y_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     i_min = bc_sections(2,k)
                     i_max = bc_sections(4,k)+1

                     call this%compute_fluxes_for_bc_y_edge(
     $                    p_model,
     $                    nodes,
     $                    s_y_L0, s_y_L1,
     $                    s_y_R1, s_y_R0,
     $                    dx, dy,
     $                    i_min, i_max, j,
     $                    S,
     $                    flux_x)

                  end if

               case(E_edge_type)

                  !do not compute the edge fluxes only if
                  !(x.ge.bc_x_max).and.
                  !(bc_E_type_choice.ne.bc_timedev_choice)

                  i = bc_sections(2,k)

                  compute_edge = compute_edge_E(i,x_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     j_min = bc_sections(3,k)
                     j_max = bc_sections(4,k)+1

                     call this%compute_fluxes_for_bc_x_edge(
     $                    p_model,
     $                    nodes,
     $                    s_x_L0, s_x_L1,
     $                    s_x_R1, s_x_R0,
     $                    dx, dy,
     $                    j_min, j_max, i,
     $                    E,
     $                    flux_y)

                  end if

               case(W_edge_type)

                  !do not compute the edge fluxes only if
                  !(x.le.bc_x_min).and.
                  !(bc_W_type_choice.ne.bc_timedev_choice)

                  i = bc_sections(2,k)

                  compute_edge = compute_edge_W(i,x_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the fluxes
                  if(compute_edge) then
                  
                     j_min = bc_sections(3,k)
                     j_max = bc_sections(4,k)+1

                     call this%compute_fluxes_for_bc_x_edge(
     $                    p_model,
     $                    nodes,
     $                    s_x_L0, s_x_L1,
     $                    s_x_R1, s_x_R0,
     $                    dx, dy,
     $                    j_min, j_max, i,
     $                    W,
     $                    flux_y)

                  end if

             end select

          end do

        end subroutine compute_fluxes_bc_sections


        !compute the time derivatives using the boundary
        !conditions
        subroutine compute_timedev_bc_sections(
     $     this,
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     flux_x, flux_y,
     $     timedev,
     $     bc_sections)

          implicit none

          class(bc_operators_openbc)     , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind)                    , intent(in)    :: t
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          real(rkind), dimension(:)      , intent(in)    :: x_map
          real(rkind), dimension(:)      , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)  , intent(in)    :: flux_x
          real(rkind), dimension(:,:,:)  , intent(in)    :: flux_y
          real(rkind), dimension(:,:,:)  , intent(inout) :: timedev
          integer    , dimension(:,:)    , intent(in)    :: bc_sections

          integer(ikind) :: i,j
          integer(ikind) :: i_min, i_max
          integer(ikind) :: j_min, j_max
          integer        :: k
          logical        :: side_x, side_y
          logical        :: compute_edge


          !go through the boundary layers
          do k=1, size(bc_sections,2)

             !identify the type of boundary layer
             select case(bc_sections(1,k))

               case(N_edge_type)

                  !do not compute the time derivatives if
                  !(y.ge.bc_y_max).and.
                  !(bc_N_type_choice.ne.bc_timedev_choice)

                  j_min = bc_sections(3,k)

                  compute_edge = compute_edge_N(j_min,y_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the time derivatives
                  if(compute_edge) then
                  
                     i_min  = bc_sections(2,k)
                     i_max  = bc_sections(4,k)
                     side_y = right

                     j=j_min
                     do i=i_min,i_max

                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_y_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_x,
     $                       side_y,
     $                       gradient_y_y_oneside_R1)

                     end do

                     j=j_min+1
                     do i=i_min,i_max

                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_y_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_x,
     $                       side_y,
     $                       gradient_y_y_oneside_R0)

                     end do

                  end if
                     
               case(S_edge_type)

                  !do not compute the time derivatives if
                  !(y.le.bc_y_min).and.
                  !(bc_S_type_choice.ne.bc_timedev_choice)

                  j_min = bc_sections(3,k)

                  compute_edge = compute_edge_S(j_min,y_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the time derivatives
                  if(compute_edge) then
                  
                     i_min  = bc_sections(2,k)
                     i_max  = bc_sections(4,k)
                     side_y = left

                     j=j_min
                     do i=i_min,i_max

                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_y_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_x,
     $                       side_y,
     $                       gradient_y_y_oneside_L0)

                     end do

                     j=j_min+1
                     do i=i_min,i_max

                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_y_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_x,
     $                       side_y,
     $                       gradient_y_y_oneside_L1)

                     end do

                  end if

               case(E_edge_type)

                  !do not compute the time derivatives if
                  !(x.ge.bc_x_max).and.
                  !(bc_E_type_choice.ne.bc_timedev_choice)

                  i_min = bc_sections(2,k)

                  compute_edge = compute_edge_E(i_min,x_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the time derivatives
                  if(compute_edge) then
                  
                     j_min  = bc_sections(3,k)
                     j_max  = bc_sections(4,k)
                     side_x = right

                     do j=j_min, j_max

                        i=i_min
                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_x_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_y,
     $                       side_x,
     $                       gradient_x_x_oneside_R1)

                        i=i_min+1
                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_x_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_y,
     $                       side_x,
     $                       gradient_x_x_oneside_R0)

                     end do

                  end if

               case(W_edge_type)

                  !do not compute the time derivatives if
                  !(x.le.bc_x_min).and.
                  !(bc_W_type_choice.ne.bc_timedev_choice)

                  i_min = bc_sections(2,k)

                  compute_edge = compute_edge_W(i_min,x_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the time derivatives
                  if(compute_edge) then
                  
                     j_min  = bc_sections(3,k)
                     j_max  = bc_sections(4,k)
                     side_x = left

                     do j=j_min,j_max

                        i=i_min
                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_x_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_y,
     $                       side_x,
     $                       gradient_x_x_oneside_L0)

                        i=i_min+1
                        timedev(i,j,:) = 
     $                       this%apply_bc_on_timedev_x_edge(
     $                       p_model,
     $                       t,nodes,
     $                       x_map,y_map,i,j,
     $                       flux_y,
     $                       side_x,
     $                       gradient_x_x_oneside_L1)

                     end do

                  end if

               case(SW_corner_type,SW_edge_type)

                  i_min=bc_sections(2,k)
                  j_min=bc_sections(3,k)

                  compute_edge =
     $                 compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $                 compute_edge_W(i_min,x_map,bc_timedev_choice)

                  !determine the extent of the edge from the
                  !bc_section and compute the time derivatives
                  if(compute_edge) then

                     side_x = left
                     side_y = left

                     i=i_min
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

                     i=i_min+1
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L1,
     $                    gradient_y_y_oneside_L0)

                     i=i_min
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L1)

                     i=i_min+1
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L1,
     $                    gradient_y_y_oneside_L1)

                  end if

               case(SE_corner_type,SE_edge_type)

                  i_min=bc_sections(2,k)
                  j_min=bc_sections(3,k)

                  compute_edge =
     $                 compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $                 compute_edge_E(i_min,x_map,bc_timedev_choice)

                  if(compute_edge) then

                     side_x = right
                     side_y = left

                     i=i_min
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R1,
     $                    gradient_y_y_oneside_L0)
                     
                     i=i_min+1
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)
                     
                     i=i_min
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R1,
     $                    gradient_y_y_oneside_L1)
                     
                     i=i_min+1
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L1)

                  end if
                     
               case(NW_corner_type,NW_edge_type)

                  i_min=bc_sections(2,k)
                  j_min=bc_sections(3,k)

                  compute_edge =
     $                 compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $                 compute_edge_W(i_min,x_map,bc_timedev_choice)
                  
                  if(compute_edge) then

                     side_x = left
                     side_y = right
                     
                     i=i_min
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R1)

                     i=i_min+1
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L1,
     $                    gradient_y_y_oneside_R1)

                     i=i_min
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

                     i=i_min+1
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_L1,
     $                    gradient_y_y_oneside_R0)

                  end if

               case(NE_corner_type,NE_edge_type)

                  i_min=bc_sections(2,k)
                  j_min=bc_sections(3,k)

                  compute_edge =
     $                 compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $                 compute_edge_E(i_min,x_map,bc_timedev_choice)
                  
                  if(compute_edge) then

                     side_x = right
                     side_y = right

                     i=i_min
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R1,
     $                    gradient_y_y_oneside_R1)

                     i=i_min+1
                     j=j_min
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R1)

                     i=i_min
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R1,
     $                    gradient_y_y_oneside_R0)

                     i=i_min+1
                     j=j_min+1
                     timedev(i,j,:) =
     $                    this%apply_bc_on_timedev_xy_corner(
     $                    p_model,
     $                    t,nodes,x_map,y_map,i,j,
     $                    side_x, side_y,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

                  end if

             end select

          end do

        end subroutine compute_timedev_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the y-direction so that
        !> the time derivatives for an edge in the x-direction
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
        !>@param j_min
        !> index min along the y-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param j_max
        !> index max along the y-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param i
        !> index along the x-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_x_edge(
     $       this,
     $       p_model,
     $       nodes,
     $       s_x_L0, s_x_L1,
     $       s_x_R1, s_x_R0,
     $       dx, dy,
     $       j_min, j_max, i,
     $       edge_card_coord,
     $       flux_y)
        
          implicit none            
        
          class(bc_operators_openbc)     , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          type(sd_operators_x_oneside_L0), intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1), intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1), intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0), intent(in)    :: s_x_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: j_min
          integer(ikind)                 , intent(in)    :: j_max
          integer(ikind)                 , intent(in)    :: i
          integer                        , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_y

          integer(ikind)        :: i_f
          integer(ikind)        :: j
          integer, dimension(4) :: bc_s

          bc_s = this%bc_type


          select case(edge_card_coord)
            case(W)
               do j=j_min,j_max

                  flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i,j,
     $                 s_x_L0)

                  flux_y(i+1,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i+1,j,
     $                 s_x_L1)

               end do
               
            case(E)
               do j=j_min,j_max

                  flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i,j,
     $                 s_x_R1)

                  flux_y(i+1,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i+1,j,
     $                 s_x_R0)

               end do

            case(E+W)
               do j=j_min,j_max

                  i_f=1
                  flux_y(i_f,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i_f,j,
     $                 s_x_L0)

                  i_f=bc_size
                  flux_y(i_f,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i_f,j,
     $                 s_x_L1)

                  i_f=nx-1
                  flux_y(i_f,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i_f,j,
     $                 s_x_R1)

                  i_f=nx
                  flux_y(i_f,j,:) = p_model%compute_flux_y_oneside(
     $                 nodes,dx,dy,
     $                 i_f,j,
     $                 s_x_R0)

               end do


            case default
               print '(''bc_operators_openbc_class'')'
               print '(''compute_fluxes_for_bc_x_edge'')'
               stop 'case not recognized'
          end select

                  
        end subroutine compute_fluxes_for_bc_x_edge


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
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_y_edge(
     $     this,
     $     p_model,
     $     nodes,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     dx, dy,
     $     i_min, i_max, j,
     $     edge_card_coord,
     $     flux_x)
        
          implicit none
        
          class(bc_operators_openbc)     , intent(in)    :: this
          type(pmodel_eq)                , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)  , intent(in)    :: nodes
          type(sd_operators_y_oneside_L0), intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1), intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1), intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0), intent(in)    :: s_y_R0
          real(rkind)                    , intent(in)    :: dx
          real(rkind)                    , intent(in)    :: dy
          integer(ikind)                 , intent(in)    :: i_min
          integer(ikind)                 , intent(in)    :: i_max
          integer(ikind)                 , intent(in)    :: j
          integer                        , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)  , intent(inout) :: flux_x


          integer(ikind)        :: i
          integer, dimension(4) :: bc_s

          bc_s = this%bc_type


          select case(edge_card_coord)
            case(S)
               do i=i_min, i_max
                  flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $                 nodes,dx,dy,
     $                 i,j,
     $                 s_y_L0)
               end do

               do i=i_min, i_max
                  flux_x(i,j+1,:) = p_model%compute_flux_x_oneside(
     $                 nodes,dx,dy,
     $                 i,j+1,
     $                 s_y_L1)
               end do

            case(N)
               do i=i_min,i_max
                  flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $                 nodes,dx,dy,
     $                 i,j,
     $                 s_y_R1)
               end do
          
               do i=i_min, i_max
                  flux_x(i,j+1,:) = p_model%compute_flux_x_oneside(
     $                 nodes,dx,dy,
     $                 i,j+1,
     $                 s_y_R0)
               end do

            case default
               print '(''bc_operators_openbc_class'')'
               print '(''compute_fluxes_for_bc_y_edge'')'
               stop 'case not recognized'
          end select
        
        end subroutine compute_fluxes_for_bc_y_edge

      end module bc_operators_openbc_class 
