      !> @file
      !> class implementing the subroutines computating the LODI
      !> vectors in the x and y directions for non-reflecting inflow
      !> boundary conditions according to the procedure designed by
      !> Yoo et al. in "Characteristic boundary conditions for direct
      !> simulations of turbulent counterflow flames", Combustion Theory
      !> and Modelling, Vol 9, No. 4, pp 617 - 646, 2012
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the subroutines computing the LODI vectors
      !> for non-reflecting inflow boundary conditions as designed by Yoo
      !> et al.
      !
      !> @date
      ! 04_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_edge_inflow_class

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use lodi_component_module, only :
     $       get_incoming_acoustic_component,
     $       get_other_acoustic_component,
     $       get_sign_acoustic_component

        use lodi_edge_ns2d_class, only :
     $       lodi_edge_ns2d

        use lodi_relaxation_coeff_module, only :
     $       get_local_mach,
     $       get_relaxation_lodiT,
     $       get_relaxation_normal_velocity,
     $       get_relaxation_trans_velocity,
     $       get_relaxation_temperature,
     $       get_relaxation_pressure

        use lodi_transverse_module, only :
     $       get_enhanced_lodi

        use ns2d_prim_module, only :
     $       velocity_x,
     $       velocity_y,
     $       temperature,
     $       pressure,
     $       speed_of_sound

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_edge_inflow

        
        !>@class lodi_edge_abstract
        !> class encapsulating the interfaces for computing the
        !> LODI amplitudes in the x and y directions at the edge
        !> of the computational domain (boundary layers except the
        !> corners)
        !
        !>@param ini
        !> initialization of the functions describing the inlet
        !> or outlet flow (ex: u_in, v_in, P_out...)
        !
        !>@param compute_x_lodi
        !> compute the LODI amplitudes in the x-direction
        !
        !>@param compute_y_lodi
        !> compute the LODI amplitudes in the y-direction
        !
        !>@param compute_x_timedev
        !> compute the contribution to the time derivative of
        !> the LODI amplitudes in the x-direction
        !
        !>@param compute_y_timedev
        !> compute the contribution to the time derivative of
        !> the LODI amplitudes in the y-direction
        !---------------------------------------------------------------
        type, extends(lodi_edge_ns2d) :: lodi_edge_inflow

          character(len=19) :: title

          contains

          procedure,   pass :: ini
          procedure, nopass :: compute_x_lodi
          procedure, nopass :: compute_y_lodi

        end type lodi_edge_inflow

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

          class(lodi_edge_inflow), intent(inout) :: this

          this%title = 'inflow Yoo b.c.'

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x-direction enforcing
        !> \f$ (u_\text{in}, v_\text{in}, T_\text{in}) \f$
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
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
        !>@param transverse_lodi
        !> transverse LODI vector at (i,j)
        !
        !>@param viscous_lodi
        !> viscous LODI vector at (i,j)
        !
        !>@param side
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return lodi
        !> LODI vector
        !---------------------------------------------------------------
        function compute_x_lodi(
     $     p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     transverse_lodi, viscous_lodi,
     $     side,
     $     gradient)
     $     result(lodi)

          implicit none
          
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(ne)   , intent(in) :: transverse_lodi
          real(rkind), dimension(ne)   , intent(in) :: viscous_lodi
          logical                      , intent(in) :: side
          procedure(gradient_x_proc)                :: gradient
          real(rkind), dimension(ne)                :: lodi

          integer :: ix_in    !index of incoming acoustic LODI component
          integer :: ix_out   !index of outgoing acoustic LODI component
          integer :: sign_out !sign for the outgoing acoustic LODI component

          real(rkind) :: u_set
          real(rkind) :: v_set
          real(rkind) :: T_set

          real(rkind) :: c
          real(rkind) :: u
          real(rkind) :: v
          real(rkind) :: temp

          real(rkind)                :: L_domain_x
          real(rkind)                :: mach_local
          real(rkind)                :: mach_ux_infty
          real(rkind)                :: relaxation_u
          real(rkind)                :: relaxation_v
          real(rkind)                :: relaxation_T
          real(rkind)                :: relaxation_lodiT
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind)                :: dx
          real(rkind)                :: dPdx
          real(rkind)                :: dudx
          real(rkind)                :: normal_velocity_forcing
          real(rkind)                :: transv_velocity_forcing
          real(rkind)                :: temperature_forcing
          real(rkind)                :: outgoing_component

          
          !computation of the variables needed for the LODI component
          call get_LODI_inflow_intermediate_variables(
     $         t, nodes(i,j,:), x_map(i), y_map(j), side,
     $         p_model,
     $         ix_in, ix_out, sign_out,
     $         u_set,v_set,T_set,
     $         c,u,v)
          temp = temperature(nodes,i,j)


          !get the variables specific to the x-direction
          L_domain_x       = x_map(size(x_map,1))-x_map(1)
          mach_local       = get_local_mach(u,v,c)
          mach_ux_infty    = p_model%get_mach_ux_infty(side)

          relaxation_u     = get_relaxation_normal_velocity(L_domain_x,mach_ux_infty,side)
          relaxation_v     = get_relaxation_trans_velocity(L_domain_x, mach_local)
          relaxation_T     = get_relaxation_temperature(L_domain_x, mach_local)
          relaxation_lodiT = get_relaxation_lodiT(mach_local)

          eigenvalues      = p_model%compute_x_eigenvalues(nodes(i,j,:))
          dx               = x_map(2)-x_map(1)
          dPdx             = gradient(nodes,i,j,pressure,dx)
          dudx             = gradient(nodes,i,j,velocity_x,dx)

          normal_velocity_forcing = relaxation_u*(u-u_set)
          transv_velocity_forcing = relaxation_v*(v-v_set)
          temperature_forcing     = relaxation_T*(temp-T_set)
          outgoing_component      = eigenvalues(ix_out)*(dPdx + sign_out*nodes(i,j,1)*c*dudx)


          !computation of the LODI components
          lodi = compute_lodi_inflow_components(
     $         ix_in, ix_out,
     $         normal_velocity_forcing,
     $         transv_velocity_forcing,
     $         temperature_forcing,
     $         outgoing_component,
     $         relaxation_lodiT,
     $         transverse_lodi,
     $         viscous_lodi)

        end function compute_x_lodi


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x-direction enforcing
        !> \f$ (u_\text{in}, v_\text{in}, T_\text{in}) \f$
        !
        !> @date
        !> 04_09_2014 - initial version - J.L. Desmarais
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
        !>@param transverse_lodi
        !> transverse LODI vector at (i,j)
        !
        !>@param viscous_lodi
        !> viscous LODI vector at (i,j)
        !
        !>@param side
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return lodi
        !> LODI vector
        !---------------------------------------------------------------
        function compute_y_lodi(
     $     p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     transverse_lodi, viscous_lodi,
     $     side,
     $     gradient)
     $     result(lodi)

          implicit none
          
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(ne)   , intent(in) :: transverse_lodi
          real(rkind), dimension(ne)   , intent(in) :: viscous_lodi
          logical                      , intent(in) :: side
          procedure(gradient_y_proc)                :: gradient
          real(rkind), dimension(ne)                :: lodi

          integer :: iy_in    !index of incoming acoustic LODI component
          integer :: iy_out   !index of outgoing acoustic LODI component
          integer :: sign_out !sign for the outgoing acoustic LODI component

          real(rkind) :: u_set
          real(rkind) :: v_set
          real(rkind) :: T_set

          real(rkind) :: c
          real(rkind) :: u
          real(rkind) :: v
          real(rkind) :: temp

          real(rkind)                :: L_domain_y
          real(rkind)                :: mach_local
          real(rkind)                :: mach_uy_infty
          real(rkind)                :: relaxation_u
          real(rkind)                :: relaxation_v
          real(rkind)                :: relaxation_T
          real(rkind)                :: relaxation_lodiT
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind)                :: dy
          real(rkind)                :: dPdy
          real(rkind)                :: dvdy
          real(rkind)                :: normal_velocity_forcing
          real(rkind)                :: transv_velocity_forcing
          real(rkind)                :: temperature_forcing
          real(rkind)                :: outgoing_component

          
          !computation of the variables needed for the LODI component
          call get_LODI_inflow_intermediate_variables(
     $         t, nodes(i,j,:), x_map(i), y_map(j), side,
     $         p_model,
     $         iy_in, iy_out, sign_out,
     $         u_set,v_set,T_set,
     $         c,u,v)
          temp = temperature(nodes,i,j)


          !get the variables specific to the x-direction
          L_domain_y       = y_map(size(y_map,1))-y_map(1)
          mach_local       = get_local_mach(u,v,c)
          mach_uy_infty    = p_model%get_mach_uy_infty(side)

          relaxation_u     = get_relaxation_trans_velocity(L_domain_y, mach_local)
          relaxation_v     = get_relaxation_normal_velocity(L_domain_y,mach_uy_infty,side)
          relaxation_T     = get_relaxation_temperature(L_domain_y, mach_local)
          relaxation_lodiT = get_relaxation_lodiT(mach_local)

          eigenvalues      = p_model%compute_y_eigenvalues(nodes(i,j,:))
          dy               = y_map(2)-y_map(1)
          dPdy             = gradient(nodes,i,j,pressure,dy)
          dvdy             = gradient(nodes,i,j,velocity_y,dy)

          normal_velocity_forcing = relaxation_v*(v-v_set)
          transv_velocity_forcing = relaxation_u*(u-u_set)
          temperature_forcing     = relaxation_T*(temp-T_set)
          outgoing_component      = eigenvalues(iy_out)*(dPdy + sign_out*nodes(i,j,1)*c*dvdy)


          !computation of the LODI components
          lodi = compute_lodi_inflow_components(
     $         iy_in, iy_out,
     $         normal_velocity_forcing,
     $         transv_velocity_forcing,
     $         temperature_forcing,
     $         outgoing_component,
     $         relaxation_lodiT,
     $         transverse_lodi,
     $         viscous_lodi)

        end function compute_y_lodi


        !get the intermediate variables when computing the LODI
        !inflow components
        subroutine get_LODI_inflow_intermediate_variables(
     $     t, nodes, x, y, side,
     $     p_model,
     $     ix_in, ix_out, sign_out,
     $     u_set,v_set,T_set,
     $     c,u,v)

          implicit none

          real(rkind)               , intent(in) :: t
          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          logical                   , intent(in) :: side
          type(pmodel_eq)           , intent(in) :: p_model

          integer    , intent(out) :: ix_in
          integer    , intent(out) :: ix_out
          integer    , intent(out) :: sign_out
                              
          real(rkind), intent(out) :: u_set
          real(rkind), intent(out) :: v_set
          real(rkind), intent(out) :: T_set
                              
          real(rkind), intent(out) :: c
          real(rkind), intent(out) :: u
          real(rkind), intent(out) :: v


          !get the indices corresponding to the LODI acoustic components
          ix_in    = get_incoming_acoustic_component(side)
          ix_out   = get_other_acoustic_component(ix_in)
          sign_out = get_sign_acoustic_component(ix_out)


          !get the constrained data at the location of the boundary
          u_set = p_model%get_u_in(t,x,y)
          v_set = p_model%get_v_in(t,x,y)
          T_set = p_model%get_T_in(t,x,y)


          !get the primitive variables for the computation
          !of the LODI components
          c = speed_of_sound(nodes)
          u = nodes(2)/nodes(1)
          v = nodes(3)/nodes(1)

        end subroutine get_LODI_inflow_intermediate_variables


        !computation of the LODI inflow components
        function compute_lodi_inflow_components(
     $     ix_in, ix_out,
     $     normal_velocity_forcing,
     $     transv_velocity_forcing,
     $     temperature_forcing,
     $     outgoing_component,
     $     relaxation_lodiT,
     $     transverse_lodi,
     $     viscous_lodi)
     $     result(lodi)

          implicit none
          
          integer                   , intent(in) :: ix_in
          integer                   , intent(in) :: ix_out
          real(rkind)               , intent(in) :: normal_velocity_forcing
          real(rkind)               , intent(in) :: transv_velocity_forcing
          real(rkind)               , intent(in) :: temperature_forcing
          real(rkind)               , intent(in) :: outgoing_component
          real(rkind)               , intent(in) :: relaxation_lodiT
          real(rkind), dimension(ne), intent(in) :: transverse_lodi
          real(rkind), dimension(ne), intent(in) :: viscous_lodi
          real(rkind), dimension(ne)             :: lodi


          !computation of the first two components
          lodi(1) = transv_velocity_forcing + 
     $              get_enhanced_lodi(relaxation_lodiT,     
     $                                transverse_lodi(1),   
     $                                viscous_lodi(1))

          lodi(2) = temperature_forcing + 
     $              get_enhanced_lodi(relaxation_lodiT,  
     $                                transverse_lodi(2),
     $                                viscous_lodi(2))
                                             

          !computation of the acoustic components
          lodi(ix_in)  = normal_velocity_forcing +
     $                   get_enhanced_lodi(relaxation_lodiT,  
     $                                     transverse_lodi(ix_in),
     $                                     viscous_lodi(ix_in))
          lodi(ix_out) = outgoing_component

        end function compute_lodi_inflow_components


      end module lodi_edge_inflow_class
