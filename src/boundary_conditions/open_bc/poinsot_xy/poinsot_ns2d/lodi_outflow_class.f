      !> @file
      !> class encapsulating the subroutines computing the LODI
      !> amplitudes in the x and y directions for fixed pressure
      !> \f$ P_\infty \f$ as computed by Poinsot and Lele in
      !> "Boundary conditions for direct simulations of compressible
      !> viscous flows", J. Comput. Phys., Vol 101, No. 1, pp 104 - 129,
      !> 1992
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the subroutines computing the LODI
      !> amplitudes in the x and y directions for fixed pressure
      !> \f$ P_\infty \f$ using Poinsot and Lele's subsonic outflow b.c.
      !
      !> @date
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_outflow_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc
        
        use lodi_ns2d_class, only :
     $       lodi_ns2d

        use ns2d_parameters, only :
     $       gamma, mach_infty

        use ns2d_prim_module, only :
     $       mass_density,
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       speed_of_sound

        use parameters_constant, only :
     $       left

        use parameters_input, only :
     $       ne,
     $       sigma_P

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_outflow


        !>@class lodi_inflow_abstract
        !> class encapsulating the subroutines computing the LODI
        !> amplitudes in the x and y directions for fixed pressure
        !> \f$ P_\infty \f$ using Poinsot and Lele's subsonic outflow b.c.
        !
        !>@param ini
        !> initialization of the functions describing the outlet flow
        !> (ex: P_\infty, ...)
        !
        !>@param compute_x_lodi
        !> compute the LODI amplitudes in the x-direction
        !
        !>@param compute_y_lodi
        !> compute the LODI amplitudes in the y-direction
        !
        !>@param compute_Pout
        !> compute the outflow pressure
        !---------------------------------------------------------------
        type, extends(lodi_ns2d) :: lodi_outflow

          character(len=20) :: title
          real(rkind)       :: relaxation_P

          contains

          procedure, pass   :: ini
          procedure, pass   :: compute_x_lodi
          procedure, pass   :: compute_y_lodi

          procedure, nopass :: compute_Pout
          procedure,   pass :: get_relaxation_P
          
        end type lodi_outflow


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

          class(lodi_outflow), intent(inout) :: this

          this%title = 'outflow Poinsot b.c.'

c$$$          this%relaxation_P = 0.1249875

          if(rkind.eq.8) then
             this%relaxation_P = sigma_P*(1.0d0-(mach_infty)**2)*1.0d0/2.0d0
          else
             this%relaxation_P = sigma_P*(1.0-(mach_infty)**2)*1.0/2.0
          end if

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x-direction enforcing
        !> non-reflecting outflow b.c. at constant pressure
        !> \f$ P_\infty \f$
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
          
          class(lodi_outflow)          , intent(in) :: this
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
          real(rkind)                :: mass_grad
          real(rkind)                :: velocity_x_grad
          real(rkind)                :: velocity_y_grad
          real(rkind)                :: pressure_grad
          real(rkind)                :: P
          real(rkind)                :: P_out
          real(rkind), dimension(ne) :: eigenvalues


          c  = speed_of_sound(nodes(i,j,:))
          dx = x_map(2)-x_map(1)

          mass_grad       = gradient(nodes,i,j,mass_density,dx)
          velocity_x_grad = gradient(nodes,i,j,velocity_x,dx)
          velocity_y_grad = gradient(nodes,i,j,velocity_y,dx)
          pressure_grad   = gradient(nodes,i,j,pressure,dx)

          P               = pressure(nodes,i,j)
          P_out           = this%compute_Pout(x_map(i),y_map(j),t)

          eigenvalues     = p_model%compute_x_eigenvalues(nodes(i,j,:))


          !computation of the LODI vector
          lodi(1) = eigenvalues(1)*velocity_y_grad
          lodi(2) = eigenvalues(2)*(c**2*mass_grad-pressure_grad)
          if(side.eqv.left) then
             lodi(3) = eigenvalues(3)*(pressure_grad-nodes(i,j,1)*c*velocity_x_grad)
             lodi(4) = this%relaxation_P*(P-P_out)
          else
             lodi(3) = this%relaxation_P*(P-P_out)
             lodi(4) = eigenvalues(4)*(pressure_grad+nodes(i,j,1)*c*velocity_x_grad)
          end if

        end function compute_x_lodi


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the y-direction enforcing
        !> non-reflecting outflow b.c. at constant pressure
        !> \f$ P_\infty \f$
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
          
          class(lodi_outflow)          , intent(in) :: this
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
          real(rkind)                :: mass_grad
          real(rkind)                :: velocity_x_grad
          real(rkind)                :: velocity_y_grad
          real(rkind)                :: pressure_grad
          real(rkind)                :: P
          real(rkind)                :: P_out
          real(rkind), dimension(ne) :: eigenvalues


          c  = speed_of_sound(nodes(i,j,:))
          dy = y_map(2)-y_map(1)

          mass_grad       = gradient(nodes,i,j,mass_density,dy)
          velocity_x_grad = gradient(nodes,i,j,velocity_x,dy)
          velocity_y_grad = gradient(nodes,i,j,velocity_y,dy)
          pressure_grad   = gradient(nodes,i,j,pressure,dy)

          P               = pressure(nodes,i,j)
          P_out           = this%compute_Pout(x_map(i),y_map(j),t)

          eigenvalues     = p_model%compute_y_eigenvalues(nodes(i,j,:))
       

          !computation of the LODI vector
          lodi(1) = eigenvalues(1)*velocity_x_grad
          lodi(2) = eigenvalues(2)*(c**2*mass_grad-pressure_grad)

          if(side.eqv.left) then
             lodi(3) = eigenvalues(3)*(pressure_grad-nodes(i,j,1)*c*velocity_y_grad)
             lodi(4) = this%relaxation_P*(P-P_out)
          else
             lodi(3) = this%relaxation_P*(P-P_out)
             lodi(4) = eigenvalues(4)*(pressure_grad+nodes(i,j,1)*c*velocity_y_grad)
          end if

        end function compute_y_lodi


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the pressure in the far field for the outlet flow
        !> \f$ P_\infty(x,y,t) \f$
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
        function compute_Pout(x,y,t) result(P_out)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind), intent(in) :: t
          real(rkind)             :: P_out

          real(rkind) :: x_s,y_s,t_s

          !P_out = 1.0d0
          P_out = 1.0d0/(gamma*mach_infty**2)

          !to prevent unused arg. warning
          x_s = x
          y_s = y
          t_s = t

        end function compute_Pout

        function get_relaxation_P(this)

          implicit none

          class(lodi_outflow), intent(in) :: this
          real(rkind)                     :: get_relaxation_P

          get_relaxation_P = this%relaxation_P

        end function get_relaxation_P

      end module lodi_outflow_class
