      !> @file
      !> module encapsulating the functions shared by the open boundary
      !> conditions procedures
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the functions shared by the open boundary
      !> conditions procedures
      !
      !> @date
      ! 01_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module openbc_operators_module

        use parameters_input  , only : nx,ny,ne,bc_size
        use parameters_kind   , only : ikind, rkind
        use pmodel_eq_class   , only : pmodel_eq

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
        

        private
        public :: compute_fluxes_at_the_edges_2ndorder

        
        contains


        subroutine compute_fluxes_at_the_edges_2ndorder(
     $       nodes, dx, dy,
     $       s_x_L0, s_x_L1, s_x_R1, s_x_R0,
     $       s_y_L0, s_y_L1, s_y_R1, s_y_R0,
     $       p_model,
     $       flux_x, flux_y)

          implicit none

          real(rkind)        , dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                               , intent(in)    :: dx
          real(rkind)                               , intent(in)    :: dy
          type(sd_operators_x_oneside_L0)           , intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1)           , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)           , intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0)           , intent(in)    :: s_x_R0
          type(sd_operators_y_oneside_L0)           , intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1)           , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)           , intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0)           , intent(in)    :: s_y_R0
          type(pmodel_eq)                           , intent(in)    :: p_model
          real(rkind)        , dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind)        , dimension(nx,ny+1,ne), intent(inout) :: flux_y


          integer(ikind) :: i,j


          !south layer of y-fluxes
          j=1
          do i=bc_size+1, nx-bc_size
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_L0)
          end do

          j=bc_size
          do i=bc_size+1, nx-bc_size
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_L1)
          end do

          
          !west layer of x-fluxes
          do j=bc_size+1, ny-bc_size
             i=1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_L0)

             i=bc_size
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_L1)
          end do


          !east layer of x-fluxes
          do j=bc_size+1, ny-bc_size
             i=nx
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_R1)

             i=nx+1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_R0)
          end do


          !north layer of y-fluxes
          j=ny
          do i=bc_size+1, nx-bc_size
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_R1)
          end do

          j=ny+1
          do i=bc_size+1, nx-bc_size
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_R0)
          end do

        end subroutine compute_fluxes_at_the_edges_2ndorder

      end module openbc_operators_module
