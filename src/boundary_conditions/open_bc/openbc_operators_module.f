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

        use parameters_input, only :
     $     nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

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
        public :: 
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right,
     $       inflow_left,
     $       inflow_right,
     $       compute_fluxes_at_the_edges_2ndorder,
     $       add_body_forces


        interface

           function incoming_proc(eigenvalue) result(is_incoming)

             import rkind

             real(rkind), intent(in) :: eigenvalue
             logical                 :: is_incoming

           end function incoming_proc

        end interface

        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the type of the characteristic wave reaching
        !> the edge of the computational domain at the east or south
        !> boundary : incoming or not
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param eigenvalues
        !> eigenvalues evaluated at the edge
        !
        !>@return incoming
        !> check whether the characteristic wave is entering the
        !> computational domain
        !--------------------------------------------------------------
        function incoming_left(eigenvalue) result(is_incoming)

          implicit none

          real(rkind), intent(in) :: eigenvalue
          logical                 :: is_incoming

          is_incoming = eigenvalue.ge.0
          
        end function incoming_left


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the type of the characteristic wave reaching
        !> the edge of the computational domain at the west or north
        !> boundary: incoming or not
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param eigenvalues
        !> eigenvalues evaluated at the edge
        !
        !>@return incoming
        !> check whether the characteristic wave is entering the
        !> computational domain
        !--------------------------------------------------------------
        function incoming_right(eigenvalue) result(is_incoming)

          implicit none

          real(rkind), intent(in) :: eigenvalue
          logical                 :: is_incoming

          is_incoming = eigenvalue.le.0
          
        end function incoming_right



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the type of the characteristic wave reaching
        !> the edge of the computational domain at the east or south
        !> boundary : incoming or not
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param eigenvalues
        !> eigenvalues evaluated at the edge
        !
        !>@return incoming
        !> check whether the characteristic wave is entering the
        !> computational domain
        !--------------------------------------------------------------
        function inflow_left(eigenvalue) result(is_incoming)

          implicit none

          real(rkind), intent(in) :: eigenvalue
          logical                 :: is_incoming

          is_incoming = eigenvalue.ge.0
          
        end function inflow_left


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the type of the characteristic wave reaching
        !> the edge of the computational domain at the west or north
        !> boundary: incoming or not
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param eigenvalues
        !> eigenvalues evaluated at the edge
        !
        !>@return incoming
        !> check whether the characteristic wave is entering the
        !> computational domain
        !--------------------------------------------------------------
        function inflow_right(eigenvalue) result(is_incoming)

          implicit none

          real(rkind), intent(in) :: eigenvalue
          logical                 :: is_incoming

          is_incoming = eigenvalue.le.0
          
        end function inflow_right


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> computation of the fluxes in the boundary layers
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param s_x_L0
        !> procedure used to compute the gradient along the x-axis
        !> with no point on the left
        !
        !>@param s_x_L1
        !> procedure used to compute the gradient along the x-axis
        !> with only one point on the left
        !
        !>@param s_x_R1
        !> procedure used to compute the gradient along the x-axis
        !> with only one point on the right
        !
        !>@param s_x_R0
        !> procedure used to compute the gradient along the x-axis
        !> with no point on the right
        !
        !>@param s_y_L0
        !> procedure used to compute the gradient along the y-axis
        !> with no point on the left
        !
        !>@param s_y_L1
        !> procedure used to compute the gradient along the y-axis
        !> with only one point on the left
        !
        !>@param s_y_R1
        !> procedure used to compute the gradient along the y-axis
        !> with only one point on the right
        !
        !>@param s_y_R0
        !> procedure used to compute the gradient along the y-axis
        !> with no point on the right
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param flux_x
        !> fluxes in the x-direction
        !
        !>@param flux_y
        !> fluxes in the y-direction
        !--------------------------------------------------------------
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


          !south layer of x-fluxes
          j=1
          do i=bc_size+1, nx-bc_size+1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_L0)
          end do

          j=bc_size
          do i=bc_size+1, nx-bc_size+1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_L1)
          end do

          
          !west and east layer of y-fluxes
          do j=bc_size+1, ny-bc_size+1
             i=1
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_L0)

             i=bc_size
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_L1)
          
             i=nx-1
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_R1)

             i=nx
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_x_R0)
          end do


          !north layer of x-fluxes
          j=ny-1
          do i=bc_size+1, nx-bc_size+1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_R1)
          end do

          j=ny
          do i=bc_size+1, nx-bc_size+1
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            nodes,dx,dy,
     $            i,j,
     $            s_y_R0)
          end do


        end subroutine compute_fluxes_at_the_edges_2ndorder


        function add_body_forces(p_model, t,x,y, nodes) result(timedev)

          implicit none

          type(pmodel_eq)           , intent(in) :: p_model
          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(ne)             :: timedev
        
          integer :: k

          do k=1, ne
             timedev(k) = p_model%compute_body_forces(
     $            t,x,y,nodes,k)
          end do
          
        end function add_body_forces

      end module openbc_operators_module
