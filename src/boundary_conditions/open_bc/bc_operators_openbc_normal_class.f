      module bc_operators_openbc_normal_class
      
        use bc_operators_openbc_class, only :
     $       bc_operators_openbc

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use parameters_constant, only :
     $     N,S,E,W,
     $     left,right

        use parameters_input, only :
     $       nx,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only : 
     $       pmodel_eq

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
        public  :: bc_operators_openbc_normal


        type, abstract, extends(bc_operators_openbc) :: bc_operators_openbc_normal

          contains

          procedure, nopass :: compute_fluxes_for_bc_x_edge
          procedure, nopass :: compute_fluxes_for_bc_y_edge

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
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_y,
     $        side_x,
     $        gradient_x)
     $        result(timedev)
           
             import bc_operators_openbc_normal
             import gradient_x_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_openbc_normal), intent(in) :: this
             type(pmodel_eq)                  , intent(in) :: p_model
             real(rkind)                      , intent(in) :: t
             real(rkind), dimension(:,:,:)    , intent(in) :: nodes
             real(rkind), dimension(:)        , intent(in) :: x_map
             real(rkind), dimension(:)        , intent(in) :: y_map
             integer(ikind)                   , intent(in) :: i
             integer(ikind)                   , intent(in) :: j
             real(rkind), dimension(:,:,:)    , intent(in) :: flux_y
             logical                          , intent(in) :: side_x
             procedure(gradient_x_proc)                    :: gradient_x
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
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_x, side_y, gradient_y)
     $        result(timedev)
           
             import bc_operators_openbc_normal
             import gradient_y_proc
             import ne
             import pmodel_eq
             import ikind
             import rkind
           
             class(bc_operators_openbc_normal), intent(in) :: this
             type(pmodel_eq)                  , intent(in) :: p_model
             real(rkind)                      , intent(in) :: t
             real(rkind), dimension(:,:,:)    , intent(in) :: nodes
             real(rkind), dimension(:)        , intent(in) :: x_map
             real(rkind), dimension(:)        , intent(in) :: y_map
             integer(ikind)                   , intent(in) :: i
             integer(ikind)                   , intent(in) :: j
             real(rkind), dimension(:,:,:)    , intent(in) :: flux_x
             logical                          , intent(in) :: side_y
             procedure(gradient_y_proc)                    :: gradient_y
             real(rkind), dimension(ne)                    :: timedev
           
           end function tdev_y_edge_l
      
        end interface


        contains


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
     $       p_model,
     $       nodes,
     $       s_x_L0, s_x_L1,
     $       s_x_R1, s_x_R0,
     $       dx, dy,
     $       j_min, j_max, i,
     $       edge_card_coord,
     $       flux_y)
        
          implicit none            
        
          type(pmodel_eq)                  , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)    , intent(in)    :: nodes
          type(sd_operators_x_oneside_L0)  , intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1)  , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)  , intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0)  , intent(in)    :: s_x_R0
          real(rkind)                      , intent(in)    :: dx
          real(rkind)                      , intent(in)    :: dy
          integer(ikind)                   , intent(in)    :: j_min
          integer(ikind)                   , intent(in)    :: j_max
          integer(ikind)                   , intent(in)    :: i
          integer                          , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)    , intent(inout) :: flux_y

          integer(ikind)        :: i_f
          integer(ikind)        :: j


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
     $     p_model,
     $     nodes,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     dx, dy,
     $     i_min, i_max, j,
     $     edge_card_coord,
     $     flux_x)
        
          implicit none
        
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



        subroutine apply_bc_on_timedev_N_edge(
     $       this,
     $       p_model,
     $       t,nodes,
     $       x_map, y_map,
     $       flux_x,
     $       s_y_L0, s_y_L1,
     $       s_y_R1, s_y_R0,
     $       dx,dy,
     $       i_min, i_max, j_min,
     $       timedev)

          implicit none

          class(bc_operators_openbc_normal), intent(in)    :: this
          type(pmodel_eq)                  , intent(in)    :: p_model
          real(rkind)                      , intent(in)    :: t
          real(rkind), dimension(:,:,:)    , intent(in)    :: nodes
          real(rkind), dimension(:)        , intent(in)    :: x_map
          real(rkind), dimension(:)        , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)    , intent(inout) :: flux_x
          type(sd_operators_y_oneside_L0)  , intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1)  , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)  , intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0)  , intent(in)    :: s_y_R0
          real(rkind)                      , intent(in)    :: dx
          real(rkind)                      , intent(in)    :: dy
          integer(ikind)                   , intent(in)    :: i_min
          integer(ikind)                   , intent(in)    :: i_max
          integer(ikind)                   , intent(in)    :: j_min
          real(rkind), dimension(:,:,:)    , intent(inout) :: timedev

          logical        :: side_y
          integer(ikind) :: i,j


          !compute the fluxes at the edges
          call compute_fluxes_for_bc_y_edge(
     $         p_model,
     $         nodes,
     $         s_y_L0, s_y_L1,
     $         s_y_R1, s_y_R0,
     $         dx, dy,
     $         i_min, i_max+1, j_min,
     $         N,
     $         flux_x)


          !compute the time derivatives
          side_y = right
          
          j=j_min
          do i=i_min,i_max
             
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_y_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_x,
     $            side_y,
     $            gradient_y_y_oneside_R1)
             
          end do
          
          j=j_min+1
          do i=i_min,i_max
             
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_y_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_x,
     $            side_y,
     $            gradient_y_y_oneside_R0)

          end do

        end subroutine apply_bc_on_timedev_N_edge


        subroutine apply_bc_on_timedev_S_edge(
     $     this,
     $     p_model,
     $     t,nodes,
     $     x_map, y_map,
     $     flux_x,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     dx,dy,
     $     i_min, i_max, j_min,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal), intent(in)    :: this
          type(pmodel_eq)                  , intent(in)    :: p_model
          real(rkind)                      , intent(in)    :: t
          real(rkind), dimension(:,:,:)    , intent(in)    :: nodes
          real(rkind), dimension(:)        , intent(in)    :: x_map
          real(rkind), dimension(:)        , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)    , intent(inout) :: flux_x
          type(sd_operators_y_oneside_L0)  , intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1)  , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)  , intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0)  , intent(in)    :: s_y_R0
          real(rkind)                      , intent(in)    :: dx
          real(rkind)                      , intent(in)    :: dy
          integer(ikind)                   , intent(in)    :: i_min
          integer(ikind)                   , intent(in)    :: i_max
          integer(ikind)                   , intent(in)    :: j_min
          real(rkind), dimension(:,:,:)    , intent(inout) :: timedev

          logical        :: side_y
          integer(ikind) :: i,j


          !compute the fluxes at the edges
          call compute_fluxes_for_bc_y_edge(
     $         p_model,
     $         nodes,
     $         s_y_L0, s_y_L1,
     $         s_y_R1, s_y_R0,
     $         dx, dy,
     $         i_min, i_max+1, j_min,
     $         S,
     $         flux_x)


          !compute the time derivatives
          side_y = left

          j=j_min
          do i=i_min,i_max

             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_y_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_x,
     $            side_y,
     $            gradient_y_y_oneside_L0)

          end do

          j=j_min+1
          do i=i_min,i_max

             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_y_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_x,
     $            side_y,
     $            gradient_y_y_oneside_L1)

          end do

        end subroutine apply_bc_on_timedev_S_edge


        subroutine apply_bc_on_timedev_E_edge(
     $     this,
     $     p_model,
     $     t,nodes,
     $     x_map, y_map,
     $     flux_y,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     dx,dy,
     $     j_min, j_max, i_min,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal), intent(in)    :: this
          type(pmodel_eq)                  , intent(in)    :: p_model
          real(rkind)                      , intent(in)    :: t
          real(rkind), dimension(:,:,:)    , intent(in)    :: nodes
          real(rkind), dimension(:)        , intent(in)    :: x_map
          real(rkind), dimension(:)        , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)    , intent(inout) :: flux_y
          type(sd_operators_x_oneside_L0)  , intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1)  , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)  , intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0)  , intent(in)    :: s_x_R0
          real(rkind)                      , intent(in)    :: dx
          real(rkind)                      , intent(in)    :: dy
          integer(ikind)                   , intent(in)    :: j_min
          integer(ikind)                   , intent(in)    :: j_max
          integer(ikind)                   , intent(in)    :: i_min
          real(rkind), dimension(:,:,:)    , intent(inout) :: timedev

          logical        :: side_x
          integer(ikind) :: i,j


          !compute the fluxes at the edges
          call compute_fluxes_for_bc_x_edge(
     $         p_model,
     $         nodes,
     $         s_x_L0, s_x_L1,
     $         s_x_R1, s_x_R0,
     $         dx, dy,
     $         j_min, j_max+1, i_min,
     $         E,
     $         flux_y)

          !compute the time derivatives
          side_x = right
          
          do j=j_min, j_max
             
             i=i_min
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_x_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_y,
     $            side_x,
     $            gradient_x_x_oneside_R1)
             
             i=i_min+1
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_x_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_y,
     $            side_x,
     $            gradient_x_x_oneside_R0)
             
          end do

        end subroutine apply_bc_on_timedev_E_edge


        subroutine apply_bc_on_timedev_W_edge(
     $     this,
     $     p_model,
     $     t,nodes,
     $     x_map, y_map,
     $     flux_y,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     dx,dy,
     $     j_min, j_max, i_min,
     $     timedev)

          implicit none

          class(bc_operators_openbc_normal), intent(in)    :: this
          type(pmodel_eq)                  , intent(in)    :: p_model
          real(rkind)                      , intent(in)    :: t
          real(rkind), dimension(:,:,:)    , intent(in)    :: nodes
          real(rkind), dimension(:)        , intent(in)    :: x_map
          real(rkind), dimension(:)        , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)    , intent(inout) :: flux_y
          type(sd_operators_x_oneside_L0)  , intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1)  , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)  , intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0)  , intent(in)    :: s_x_R0
          real(rkind)                      , intent(in)    :: dx
          real(rkind)                      , intent(in)    :: dy
          integer(ikind)                   , intent(in)    :: j_min
          integer(ikind)                   , intent(in)    :: j_max
          integer(ikind)                   , intent(in)    :: i_min
          real(rkind), dimension(:,:,:)    , intent(inout) :: timedev

          logical        :: side_x
          integer(ikind) :: i,j


          !compute the fluxes
          call compute_fluxes_for_bc_x_edge(
     $         p_model,
     $         nodes,
     $         s_x_L0, s_x_L1,
     $         s_x_R1, s_x_R0,
     $         dx, dy,
     $         j_min, j_max+1, i_min,
     $         W,
     $         flux_y)

          !compute the time derivatives
          side_x = left

          do j=j_min,j_max

             i=i_min
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_x_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_y,
     $            side_x,
     $            gradient_x_x_oneside_L0)

             i=i_min+1
             timedev(i,j,:) = 
     $            this%apply_bc_on_timedev_x_edge(
     $            p_model,
     $            t,nodes,
     $            x_map,y_map,i,j,
     $            flux_y,
     $            side_x,
     $            gradient_x_x_oneside_L1)

          end do

        end subroutine apply_bc_on_timedev_W_edge
          
      end module bc_operators_openbc_normal_class
