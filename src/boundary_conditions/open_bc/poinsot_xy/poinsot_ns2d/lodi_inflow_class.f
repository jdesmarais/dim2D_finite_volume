      !> @file
      !> class encapsulating the subroutines computing the LODI
      !> amplitudes in the x and y directions for fixed
      !> \f$ (u_{\text{in}},v_{\text{in}},T_{\text{in}}) \f$
      !> as computed by Poinsot and Lele in "Boundary conditions for
      !> direct simulations of compressible viscous flows",
      !> J. Comput. Phys., Vol 101, No. 1, pp 104 - 129, 1992
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the subroutines computing the LODI
      !> amplitudes in the x and y directions for fixed
      !> \f$ (u_{\text{in}},v_{\text{in}},T_{\text{in}}) \f$ using
      !> Poinsot and Lele's subsonic inflow b.c.
      !
      !> @date
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_inflow_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc
        
        use lodi_class, only :
     $       lodi

        use ns2d_parameters, only :
     $       gamma

        use ns2d_prim_module, only :
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       speed_of_sound

        use parameters_constant, only :
     $       left

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_inflow


        !>@class lodi_inflow_abstract
        !> class encapsulating the subroutines computing the LODI
        !> amplitudes in the x and y directions for fixed
        !> \f$ (u_{\text{in}},v_{\text{in}},T_{\text{in}}) \f$ using
        !> Poinsot and Lele's subsonic inflow b.c.
        !
        !>@param ini
        !> initialization of the functions describing the inlet flow
        !> (ex: u_in, v_in, ...)
        !
        !>@param compute_x_lodi
        !> compute the LODI amplitudes in the x-direction
        !
        !>@param compute_y_lodi
        !> compute the LODI amplitudes in the y-direction
        !
        !>@param compute_duindt
        !> compute the time derivative of the velocity along the x-axis
        !> of the inlet flow
        !
        !>@param compute_dvindt
        !> compute the time derivative of the velocity along the y-axis
        !> of the inlet flow
        !
        !>@param compute_dTindt
        !> compute the time derivative of the temperature of the inlet
        !> flow
        !---------------------------------------------------------------
        type, extends(lodi) :: lodi_inflow

          character(len=19) :: title

          contains

          procedure, pass   :: ini
          procedure, pass   :: compute_x_lodi
          procedure, pass   :: compute_y_lodi

          procedure, nopass :: compute_duindt
          procedure, nopass :: compute_dvindt
          procedure, nopass :: compute_dTindt          

        end type lodi_inflow


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the attributes of the object lodi_inflow
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !---------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(lodi_inflow), intent(inout) :: this

          this%title = 'inflow Poinsot b.c.'

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x-direction enforcing
        !> (u_\text{in}, v_\text{in}, T_\text{in})
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object with the inlet conditions
        !
        !>@param p_model
        !> physical model needed to compute the eigenvalues
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> map for the coordinates along the x-direction
        !
        !>@param y_map
        !> map for the coordinates along the y-direction
        !
        !>@param i
        !> index identifying the grid point along the x-direction
        !
        !>@param j
        !> index identifying the grid point along the y-direction
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return lodi
        !> LODI vector
        !---------------------------------------------------------------
        function compute_x_lodi(
     $     this, p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     side,
     $     gradient)
     $     result(lodi)

          implicit none
          
          class(lodi_inflow)           , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_x_proc)                :: gradient
          real(rkind), dimension(ne)                :: lodi

          real(rkind)                :: c
          real(rkind)                :: dx
          real(rkind)                :: velocity_x_grad
          real(rkind)                :: pressure_grad
          real(rkind)                :: duindt
          real(rkind)                :: dvindt
          real(rkind)                :: dTindt
          real(rkind), dimension(ne) :: eigenvalues


          c  = speed_of_sound(nodes(i,j,:))
          dx = x_map(2)-x_map(1)

          velocity_x_grad = gradient(nodes,i,j,velocity_x,dx)
          pressure_grad   = gradient(nodes,i,j,pressure,dx)

          duindt          = this%compute_duindt(x_map(i),y_map(j),t)
          dvindt          = this%compute_dvindt(x_map(i),y_map(j),t)
          dTindt          = this%compute_dTindt(x_map(i),y_map(j),t)

          eigenvalues     = p_model%compute_x_eigenvalues(nodes(i,j,:))
       

          !computation of the LODI vector
          if(rkind.eq.8) then
             if(side.eqv.left) then
                lodi(3) = eigenvalues(3)*(pressure_grad - nodes(i,j,1)*c*velocity_x_grad)
                lodi(4) = lodi(3) - 2.0d0*nodes(i,j,1)*c*duindt
             else
                lodi(4) = eigenvalues(4)*(pressure_grad + nodes(i,j,1)*c*velocity_x_grad)
                lodi(3) = lodi(4) + 2.0d0*nodes(i,j,1)*c*duindt
             end if
             lodi(2) = 0.5d0*(gamma-1.0d0)*(lodi(3)+lodi(4)) + nodes(i,j,1)*c**2/temperature(nodes,i,j)*dTindt
             lodi(1) = - dvindt
          else
             if(side.eqv.left) then
                lodi(3) = eigenvalues(3)*(pressure_grad - nodes(i,j,1)*c*velocity_x_grad)
                lodi(4) = lodi(3) - 2.0*nodes(i,j,1)*c*duindt
             else
                lodi(4) = eigenvalues(4)*(pressure_grad + nodes(i,j,1)*c*velocity_x_grad)
                lodi(3) = lodi(4) + 2.0*nodes(i,j,1)*c*duindt
             end if
             lodi(2) = 0.5*(gamma-1.0)*(lodi(3)+lodi(4)) + nodes(i,j,1)*c**2/temperature(nodes,i,j)*dTindt
             lodi(1) = - dvindt
          end if

        end function compute_x_lodi


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x-direction enforcing
        !> (u_\text{in}, v_\text{in}, T_\text{in})
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object with the inlet conditions
        !
        !>@param p_model
        !> physical model needed to compute the eigenvalues
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> map for the coordinates along the x-direction
        !
        !>@param y_map
        !> map for the coordinates along the y-direction
        !
        !>@param i
        !> index identifying the grid point along the x-direction
        !
        !>@param j
        !> index identifying the grid point along the y-direction
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return lodi
        !> LODI vector
        !---------------------------------------------------------------
        function compute_y_lodi(
     $     this, p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     side,
     $     gradient)
     $     result(lodi)

          implicit none
          
          class(lodi_inflow)           , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_y_proc)                :: gradient
          real(rkind), dimension(ne)                :: lodi

          real(rkind)                :: c
          real(rkind)                :: dy
          real(rkind)                :: velocity_y_grad
          real(rkind)                :: pressure_grad
          real(rkind)                :: duindt
          real(rkind)                :: dvindt
          real(rkind)                :: dTindt
          real(rkind), dimension(ne) :: eigenvalues


          c  = speed_of_sound(nodes(i,j,:))
          dy = y_map(2)-y_map(1)

          velocity_y_grad = gradient(nodes,i,j,velocity_y,dy)
          pressure_grad   = gradient(nodes,i,j,pressure,dy)

          duindt          = this%compute_duindt(x_map(i),y_map(j),t)
          dvindt          = this%compute_dvindt(x_map(i),y_map(j),t)
          dTindt          = this%compute_dTindt(x_map(i),y_map(j),t)

          eigenvalues     = p_model%compute_y_eigenvalues(nodes(i,j,:))


          !computation of the LODI vector
          if(rkind.eq.8) then
             if(side.eqv.left) then
                lodi(3) = eigenvalues(3)*(pressure_grad - nodes(i,j,1)*c*velocity_y_grad)
                lodi(4) = lodi(3) - 2.0d0*nodes(i,j,1)*c*dvindt
             else
                lodi(4) = eigenvalues(4)*(pressure_grad + nodes(i,j,1)*c*velocity_y_grad)
                lodi(3) = lodi(4) + 2.0d0*nodes(i,j,1)*c*dvindt
             end if
             lodi(2) = 0.5d0*(gamma-1.0d0)*(lodi(3)+lodi(4)) + nodes(i,j,1)*c**2/temperature(nodes,i,j)*dTindt
             lodi(1) = - duindt
          else
             if(side.eqv.left) then
                lodi(3) = eigenvalues(3)*(pressure_grad - nodes(i,j,1)*c*velocity_y_grad)
                lodi(4) = lodi(3) - 2.0*nodes(i,j,1)*c*dvindt
             else
                lodi(4) = eigenvalues(4)*(pressure_grad + nodes(i,j,1)*c*velocity_y_grad)
                lodi(3) = lodi(4) + 2.0*nodes(i,j,1)*c*dvindt
             end if
             lodi(2) = 0.5*(gamma-1.0)*(lodi(3)+lodi(4)) + nodes(i,j,1)*c**2/temperature(nodes,i,j)*dTindt
             lodi(1) = - duindt
          end if

        end function compute_y_lodi


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivative of the velocity along the
        !> x-direction for the inlet flow
        !> \f$ \frac{d u_{\text{in}}}{dt} \f$
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@param t
        !> time
        !---------------------------------------------------------------
        function compute_duindt(x,y,t) result(duindt)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind), intent(in) :: t
          real(rkind)             :: duindt

          real(rkind) :: x_s,y_s,t_s

          duindt = 0.0d0

          !to prevent unused arg. warning
          x_s = x
          y_s = y
          t_s = t

        end function compute_duindt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivative of the velocity along the
        !> x-direction for the inlet flow
        !> \f$ \frac{d v_{\text{in}}}{dt} \f$
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@param t
        !> time
        !---------------------------------------------------------------
        function compute_dvindt(x,y,t) result(dvindt)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind), intent(in) :: t
          real(rkind)             :: dvindt

          real(rkind) :: x_s,y_s,t_s

          dvindt = 0.0d0

          !to prevent unused arg. warning
          x_s = x
          y_s = y
          t_s = t

        end function compute_dvindt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivative of the temperature for
        !> the inlet flow
        !> \f$ \frac{d T_{\text{in}}}{dt} \f$
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@param t
        !> time
        !---------------------------------------------------------------
        function compute_dTindt(x,y,t) result(dTindt)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind), intent(in) :: t
          real(rkind)             :: dTindt

          real(rkind) :: x_s,y_s,t_s

          dTindt = 0.0d0

          !to prevent unused arg. warning
          x_s = x
          y_s = y
          t_s = t

        end function compute_dTindt

      end module lodi_inflow_class
