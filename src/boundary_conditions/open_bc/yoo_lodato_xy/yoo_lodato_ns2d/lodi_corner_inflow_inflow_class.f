      !> @file
      !> class implementing the subroutines for the computation
      !> of the LODI amplitudes in the x and y directions for the
      !> inflow/inflow type corners according to the procedure
      !> designed by Yoo et al.
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class implementing the subroutines for the computation
      !> of the LODI amplitudes in the x and y directions for the
      !> inflow/inflow type corners according to the procedure
      !> designed by Yoo et al.
      !
      !> @date
      ! 08_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_corner_inflow_inflow_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use lodi_component_module, only :
     $       get_incoming_acoustic_component,
     $       get_other_acoustic_component,
     $       get_sign_acoustic_component

        use lodi_corner_ns2d_class, only :
     $       lodi_corner_ns2d

        use lodi_relaxation_coeff_module, only :
     $       get_relaxation_normal_velocity,
     $       get_relaxation_trans_velocity,
     $       get_relaxation_temperature,
     $       get_relaxation_lodiT,
     $       get_local_mach

        use ns2d_prim_module, only :
     $       pressure,
     $       velocity_x,
     $       velocity_y,
     $       temperature,
     $       speed_of_sound

        use parameters_constant, only :
     $       inflow_type

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_corner_inflow_inflow,
     $            get_lodi_A_inflow_inflow


        !>@class lodi_corner_inflow_inflow
        !> class implementing the subroutines for the computation
        !> of the LODI amplitudes in the x and y directions for the
        !> inflow/inflow type corners according to the procedure
        !> designed by Yoo et al.
        !
        !>@param ini
        !> initialization of the functions describing the inlet
        !> or outlet flow (ex: u_in, v_in, P_out...)
        !
        !>@param compute_x_and_y_lodi
        !> compute the LODI amplitudes in the x-and y-directions
        !---------------------------------------------------------------
        type, extends(lodi_corner_ns2d) :: lodi_corner_inflow_inflow

          character(len=22) :: title

          contains

          procedure,   pass :: ini
          procedure, nopass :: compute_x_and_y_lodi

        end type lodi_corner_inflow_inflow


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the attributes of the object lodi_edge_inflow
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
        !---------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(lodi_corner_inflow_inflow), intent(inout) :: this

          this%title = 'inflow/inflow Yoo b.c.'

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x- and y-direction
        !> enforcing \f$ (u_\text{in}, v_\text{in}, T_\text{in}) \f$
        !
        !> @date
        !> 08_09_2014 - initial version - J.L. Desmarais
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
        !>@param side_x
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain along
        !> the x-direction
        !
        !>@param side_y
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain along
        !> the y-direction
        !
        !>@param flow_x
        !> boolean designating the flow type at the x-boundary
        !
        !>@param flow_y
        !> boolean designating the flow type at the y-boundary
        !
        !>@param gradient_x
        !> procedure for the gradient computation along the x-direction
        !
        !>@param gradient_y
        !> procedure for the gradient computation along the y-direction
        !
        !>@param lodi_x
        !> LODI vector along the x-direction
        !
        !>@param lodi_y
        !> LODI vector along the y-direction
        !---------------------------------------------------------------
        subroutine compute_x_and_y_lodi(
     $     p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     side_x, side_y,
     $     flow_x, flow_y,
     $     gradient_x, gradient_y,
     $     lodi_x, lodi_y)

          implicit none
          
          type(pmodel_eq)              , intent(in)  :: p_model
          real(rkind)                  , intent(in)  :: t
          real(rkind), dimension(:,:,:), intent(in)  :: nodes
          real(rkind), dimension(:)    , intent(in)  :: x_map
          real(rkind), dimension(:)    , intent(in)  :: y_map
          integer(ikind)               , intent(in)  :: i
          integer(ikind)               , intent(in)  :: j
          logical                      , intent(in)  :: side_x
          logical                      , intent(in)  :: side_y
          logical                      , intent(in)  :: flow_x
          logical                      , intent(in)  :: flow_y
          procedure(gradient_x_proc)                 :: gradient_x
          procedure(gradient_y_proc)                 :: gradient_y
          real(rkind), dimension(ne)   , intent(out) :: lodi_x
          real(rkind), dimension(ne)   , intent(out) :: lodi_y

          real(rkind) :: u
          real(rkind) :: v
          real(rkind) :: temp
          real(rkind) :: c
          
          integer :: ix_in
          integer :: ix_out
          integer :: sign_x_out

          integer :: iy_in
          integer :: iy_out
          integer :: sign_y_out

          real(rkind), dimension(ne) :: eigenvalues
          real(rkind)                :: dx
          real(rkind)                :: dPdx
          real(rkind)                :: dudx
          real(rkind)                :: dy
          real(rkind)                :: dPdy
          real(rkind)                :: dvdy          

          real(rkind) :: u_set
          real(rkind) :: v_set
          real(rkind) :: T_set

          real(rkind) :: L_domain_x
          real(rkind) :: L_domain_y
          real(rkind) :: mach_local
          real(rkind) :: mach_ux_infty
          real(rkind) :: mach_uy_infty

          real(rkind) :: relaxation_u_x
          real(rkind) :: relaxation_v_x
          real(rkind) :: relaxation_T_x
          real(rkind) :: relaxation_u_y
          real(rkind) :: relaxation_v_y
          real(rkind) :: relaxation_T_y
          real(rkind) :: relaxation_lodiT

          real(rkind) :: velocity_x_forcing_x
          real(rkind) :: velocity_y_forcing_x
          real(rkind) :: temperature_forcing_x
          real(rkind) :: velocity_x_forcing_y
          real(rkind) :: velocity_y_forcing_y
          real(rkind) :: temperature_forcing_y

          real(rkind), dimension(6)   :: lodi_forcing
          real(rkind), dimension(6,6) :: lodi_A


          !test whether the boundary is well of type inflow/inflow
          if(.not.(
     $         (flow_x.eqv.inflow_type).and.
     $         (flow_y.eqv.inflow_type))) then
             print '(''lodi_corner_inflow_inflow'')'
             print '(''compute_x_and_y_lodi'')'
             print '(''flow does not have the correct configuration'')'
             print '(''flow_x:'',L2)', flow_x
             print '(''flow_y:'',L2)', flow_y
             stop
          end if


          !primitive variables
          u    = nodes(i,j,2)/nodes(i,j,1)
          v    = nodes(i,j,3)/nodes(i,j,1)
          temp = temperature(nodes,i,j)
          c    = speed_of_sound(nodes(i,j,:))


          !get the indices corresponding to the LODI acoustic components
          ix_in      = get_incoming_acoustic_component(side_x)
          ix_out     = get_other_acoustic_component(ix_in)
          sign_x_out = get_sign_acoustic_component(ix_out)

          iy_in      = get_incoming_acoustic_component(side_y)
          iy_out     = get_other_acoustic_component(iy_in)
          sign_y_out = get_sign_acoustic_component(iy_out)


          !compute the outgoing LODI acoustics components in both directions
          eigenvalues      = p_model%compute_x_eigenvalues(nodes(i,j,:))
          dx               = x_map(2)-x_map(1)
          dPdx             = gradient_x(nodes,i,j,pressure,dx)
          dudx             = gradient_x(nodes,i,j,velocity_x,dx)
          lodi_x(ix_out)   = eigenvalues(ix_out)*(dPdx + sign_x_out*nodes(i,j,1)*c*dudx)

          eigenvalues      = p_model%compute_y_eigenvalues(nodes(i,j,:))
          dy               = y_map(2)-y_map(1)
          dPdy             = gradient_y(nodes,i,j,pressure,dy)
          dvdy             = gradient_y(nodes,i,j,velocity_y,dy)
          lodi_y(iy_out)   = eigenvalues(iy_out)*(dPdy + sign_y_out*nodes(i,j,1)*c*dvdy)


          !get the set values
          u_set = p_model%get_u_in(t,x_map(i),y_map(j))
          v_set = p_model%get_v_in(t,x_map(i),y_map(j))
          T_set = p_model%get_T_in(t,x_map(i),y_map(j))


          !get the domain extensions and the mach numbers
          !in both directions
          L_domain_x       = x_map(size(x_map,1))-x_map(1)
          L_domain_y       = y_map(size(y_map,1))-y_map(1)
          mach_local       = get_local_mach(u,v,c)
          mach_ux_infty    = p_model%get_mach_ux_infty()
          mach_uy_infty    = p_model%get_mach_uy_infty()


          !get the relaxation coefficients
          relaxation_u_x = get_relaxation_normal_velocity(L_domain_x,mach_ux_infty,side_x)
          relaxation_v_x = get_relaxation_trans_velocity(L_domain_x, mach_local)
          relaxation_T_x = get_relaxation_temperature(L_domain_x, mach_local)

          relaxation_u_y = get_relaxation_trans_velocity(L_domain_y, mach_local)
          relaxation_v_y = get_relaxation_normal_velocity(L_domain_y,mach_uy_infty,side_y)
          relaxation_T_y = get_relaxation_temperature(L_domain_y, mach_local)

          relaxation_lodiT = get_relaxation_lodiT(mach_local)


          !get the forcing terms
          velocity_x_forcing_x  = relaxation_u_x*(u-u_set)
          velocity_y_forcing_x  = relaxation_v_x*(v-v_set)
          temperature_forcing_x = relaxation_T_x*(temp-T_set)

          velocity_x_forcing_y  = relaxation_u_y*(u-u_set)
          velocity_y_forcing_y  = relaxation_v_y*(v-v_set)
          temperature_forcing_y = relaxation_T_y*(temp-T_set)


          !compute the forcing LODI vector
          lodi_forcing(1) = velocity_y_forcing_x -
     $                      0.5d0*(1.0d0 - relaxation_lodiT)*
     $                      sign_y_out/(nodes(i,j,1)*c)*lodi_y(iy_out)

          lodi_forcing(2) = temperature_forcing_x

          lodi_forcing(3) = velocity_x_forcing_x -
     $                      0.5d0*(1.0d0 - relaxation_lodiT)*lodi_y(iy_out)

          lodi_forcing(4) = velocity_x_forcing_y -
     $                      0.5d0*(1.0d0 - relaxation_lodiT)*
     $                      sign_x_out/(nodes(i,j,1)*c)*lodi_x(ix_out)

          lodi_forcing(5) = temperature_forcing_y

          lodi_forcing(6) = velocity_y_forcing_y -
     $                      0.5d0*(1.0d0 - relaxation_lodiT)*lodi_x(ix_out)
          

          !get the matrix uncoupling the computation of the
          !LODI vectors
          lodi_A = get_lodi_A_inflow_inflow(
     $         nodes(i,j,1), c,
     $         -sign_x_out, -sign_y_out,
     $         relaxation_lodiT)

          !compute the incoming LODI components
          lodi_forcing = MATMUL(lodi_forcing,lodi_A)


          !transfer the computed LODI components
          !to the lodi x and y vectors
          lodi_x(1)     = lodi_forcing(1)
          lodi_x(2)     = lodi_forcing(2)
          lodi_x(ix_in) = lodi_forcing(3)

          lodi_y(1)     = lodi_forcing(4)
          lodi_y(2)     = lodi_forcing(5)
          lodi_y(iy_in) = lodi_forcing(6)

        end subroutine compute_x_and_y_lodi


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the matrix uncoupling the LODI system for
        !> inflow/inflow corners
        !
        !> @date
        !> 08_09_2014 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !
        !>@param c
        !> speed of sound
        !
        !>@param sign_x
        !> s(ix)
        !>
        !>@param sign_y
        !> s(iy)
        !>
        !>@param tau
        !> relaxation coefficient \f$ \alphda_\tau \f$ for the LODI
        !> transverse terms
        !---------------------------------------------------------------
        function get_lodi_A_inflow_inflow(md,c,sign_x,sign_y,tau)
     $     result(lodi_A)

          implicit none

          real(rkind), intent(in) :: md
          real(rkind), intent(in) :: c
          integer    , intent(in) :: sign_x
          integer    , intent(in) :: sign_y
          real(rkind), intent(in) :: tau

          real(rkind), dimension(6,6) :: lodi_A


          real(rkind) :: a0
          real(rkind) :: a1
          real(rkind) :: a2
          real(rkind) :: a3
          real(rkind) :: a4


          if(rkind.eq.8) then

             a0 = tau - 1.0d0
             a1 = -3.0d0*tau**2+6.0d0*tau+1.0d0
             a2 = - tau**2 + 2.0d0*tau + 1.0d0
             a3 = tau*(1.0d0+tau)*(tau-2.0d0)*(tau-3.0d0)
             a4 = tau*(2.0d0-tau)

             lodi_A(1,1) = a1/a3
             lodi_A(2,1) = 0.0d0
             lodi_A(3,1) = sign_y*a0**2/(md*c*a3)
             lodi_A(4,1) = a0**3*sign_x*sign_y/a3
             lodi_A(5,1) = 0.0d0
             lodi_A(6,1) = sign_y*a0*a2/(md*c*a3)

             lodi_A(1,2) = 0.0d0
             lodi_A(2,2) = 1.0d0/a4
             lodi_A(3,2) = 0.0d0
             lodi_A(4,2) = 0.0d0
             lodi_A(5,2) = a0/a4
             lodi_A(6,2) = 0.0d0

             lodi_A(1,3) = 2.0d0*sign_y*md*c*a0**2/a3
             lodi_A(2,3) = 0.0d0
             lodi_A(3,3) = 2.0d0*a2/a3
             lodi_A(4,3) = 2.0d0*sign_x*md*c*a0*a2/a3
             lodi_A(5,3) = 0.0d0
             lodi_A(6,3) = 2.0d0*a0/a3

             lodi_A(1,4) = a0**3*sign_x*sign_y/a3
             lodi_A(2,4) = 0.0d0
             lodi_A(3,4) = a0*a2*sign_x/(a3*md*c)
             lodi_A(4,4) = a1/a3
             lodi_A(5,4) = 0.0d0
             lodi_A(6,4) = a0**2*sign_x/(md*c*a3)

             lodi_A(1,5) = 0.0d0
             lodi_A(2,5) = a0/a4
             lodi_A(3,5) = 0.0d0
             lodi_A(4,5) = 0.0d0
             lodi_A(5,5) = 1.0d0/a4
             lodi_A(6,5) = 0.0d0

             lodi_A(1,6) = 2.0d0*sign_y*md*c*a0*a2/a3
             lodi_A(2,6) = 0.0d0
             lodi_A(3,6) = 2.0d0*a0/a3
             lodi_A(4,6) = 2.0d0*a0**2*sign_x*md*c/a3
             lodi_A(5,6) = 0.0d0
             lodi_A(6,6) = 2.0d0*a2/a3

          else

             a0 = tau - 1.0
             a1 = -3.0*tau**2+6.0*tau+1.0
             a2 = - tau**2 + 2.0*tau + 1.0
             a3 = tau*(1.0+tau)*(tau-2.0)*(tau-3.0)
             a4 = tau*(2.0-tau)

             lodi_A(1,1) = a1/a3
             lodi_A(2,1) = 0.0
             lodi_A(3,1) = sign_y*a0**2/(md*c*a3)
             lodi_A(4,1) = a0**3*sign_x*sign_y/a3
             lodi_A(5,1) = 0.0
             lodi_A(6,1) = sign_y*a0*a2/(md*c*a3)

             lodi_A(1,2) = 0.0
             lodi_A(2,2) = 1.0/a4
             lodi_A(3,2) = 0.0
             lodi_A(4,2) = 0.0
             lodi_A(5,2) = a0/a4
             lodi_A(6,2) = 0.0

             lodi_A(1,3) = 2.0*sign_y*md*c*a0**2/a3
             lodi_A(2,3) = 0.0
             lodi_A(3,3) = 2.0*a2/a3
             lodi_A(4,3) = 2.0*sign_x*md*c*a0*a2/a3
             lodi_A(5,3) = 0.0
             lodi_A(6,3) = 2*a0/a3

             lodi_A(1,4) = a0**3*sign_x*sign_y/a3
             lodi_A(2,4) = 0.0
             lodi_A(3,4) = a0*a2*sign_x/(a3*md*c)
             lodi_A(4,4) = a1/a3
             lodi_A(5,4) = 0.0
             lodi_A(6,4) = a0**2*sign_x/(md*c*a3)

             lodi_A(1,5) = 0.0
             lodi_A(2,5) = a0/a4
             lodi_A(3,5) = 0.0
             lodi_A(4,5) = 0.0
             lodi_A(5,5) = 1.0/a4
             lodi_A(6,5) = 0.0

             lodi_A(1,6) = 2.0*sign_y*md*c*a0*a2/a3
             lodi_A(2,6) = 0.0
             lodi_A(3,6) = 2.0*a0/a3
             lodi_A(4,6) = 2.0*a0**2*sign_x*md*c/a3
             lodi_A(5,6) = 0.0
             lodi_A(6,6) = 2.0*a2/a3

          end if

        end function get_lodi_A_inflow_inflow

      end module lodi_corner_inflow_inflow_class
