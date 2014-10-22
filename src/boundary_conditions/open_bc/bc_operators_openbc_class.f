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

        use parameters_constant, only :
     $     N,S,E,W

        use parameters_input, only :
     $     nx, bc_size

        use parameters_kind, only :
     $     ikind,
     $     rkind

        use pmodel_eq_class, only :
     $     pmodel_eq

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

           procedure, pass :: compute_fluxes_for_bc_x_edge
           procedure, pass :: compute_fluxes_for_bc_y_edge

        end type bc_operators_openbc

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

          integer(ikind) :: i_f
          integer(ikind) :: j
          integer        :: bc_s

          bc_s = this%bcx_type


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


          integer(ikind) :: i
          integer :: bc_s

          bc_s = this%bcy_type


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
