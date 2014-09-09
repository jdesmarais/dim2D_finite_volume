      !> @file
      !> class implementing the subroutines for the computation
      !> of the LODI amplitudes in the x and y directions for the
      !> inflow/outflow type corners according to the procedure
      !> designed by Yoo et al.
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class implementing the subroutines for the computation
      !> of the LODI amplitudes in the x and y directions for the
      !> inflow/outflow type corners according to the procedure
      !> designed by Yoo et al.
      !
      !> @date
      ! 09_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_corner_inflow_outflow_class

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

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
     $       mass_density,
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       speed_of_sound

        use parameters_constant, only :
     $       inflow_type,
     $       outflow_type

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_corner_inflow_outflow


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
        type, extends(lodi_corner_ns2d) :: lodi_corner_inflow_outflow

          character(len=23) :: title

          contains

          procedure,   pass :: ini
          procedure, nopass :: compute_x_and_y_lodi

        end type lodi_corner_inflow_outflow


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the attributes of the object lodi_edge_inflow
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !---------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(lodi_corner_inflow_outflow), intent(inout) :: this

          this%title = 'inflow/outflow Yoo b.c.'

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the LODI amplitudes in the x- and y-direction
        !> enforcing \f$ (u_{\text{in}}, v_{\text{in}},
        !> T_{\text{in}}) \f$ for the inflow and \f$ P_\text{out}\f$
        !> for the outflow
        !
        !> @date
        !> 09_09_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object with the inlet/outlet conditions
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

          real(rkind) :: u_set
          real(rkind) :: v_set
          real(rkind) :: T_set

          real(rkind) :: mach_local
          real(rkind) :: relaxation_lodiT

          real(rkind), dimension(ne) :: eigenvalues
          real(rkind)                :: dx
          real(rkind)                :: dy
          real(rkind)                :: dmdx
          real(rkind)                :: dudx
          real(rkind)                :: dvdx
          real(rkind)                :: dPdx
          real(rkind)                :: dmdy
          real(rkind)                :: dudy
          real(rkind)                :: dvdy
          real(rkind)                :: dPdy

          real(rkind) :: L_domain_x
          real(rkind) :: L_domain_y          
          real(rkind) :: mach_ux_infty
          real(rkind) :: mach_uy_infty
          real(rkind) :: relaxation_u_x
          real(rkind) :: relaxation_v_x
          real(rkind) :: relaxation_T_x
          real(rkind) :: relaxation_u_y
          real(rkind) :: relaxation_v_y
          real(rkind) :: relaxation_T_y


          !check flow configuration
          if(.not.(
     $         ((flow_x.eqv.inflow_type).and.(flow_y.eqv.outflow_type))
     $         .or.
     $         ((flow_x.eqv.outflow_type).and.(flow_y.eqv.inflow_type)))) then

             print '(''procedure for inflow/outflow corner'')'
             print '(''flow configuration not compatible'')'
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


          !get the set values
          u_set = p_model%get_u_in(t,x_map(i),y_map(j))
          v_set = p_model%get_v_in(t,x_map(i),y_map(j))
          T_set = p_model%get_T_in(t,x_map(i),y_map(j))


          !get the Mach numbers
          mach_local = get_local_mach(u,v,c)
          relaxation_lodiT = get_relaxation_lodiT(mach_local)


          !inflow along x-direction and outflow along y-direction
          if(flow_x.eqv.inflow_type) then

             !compute the LODI components in the y-direction
             eigenvalues      = p_model%compute_y_eigenvalues(nodes(i,j,:))
             dy               = y_map(2)-y_map(1)
             dmdy             = gradient_y(nodes,i,j,mass_density,dy)
             dudy             = gradient_y(nodes,i,j,velocity_x,dy)
             dvdy             = gradient_y(nodes,i,j,velocity_y,dy)
             dPdy             = gradient_y(nodes,i,j,pressure,dy)
             
             lodi_y(1)        = eigenvalues(1)*dudy
             lodi_y(2)        = eigenvalues(2)*(c**2*dmdy-dPdy)
             lodi_y(iy_in)    = 0.0d0
             lodi_y(iy_out)   = eigenvalues(iy_out)*(dPdy + sign_y_out*nodes(i,j,1)*c*dvdy)
             
             
             !compute the outgoing LODI acoustics component in the x-direction
             eigenvalues      = p_model%compute_x_eigenvalues(nodes(i,j,:))
             dx               = x_map(2)-x_map(1)
             dPdx             = gradient_x(nodes,i,j,pressure,dx)
             dudx             = gradient_x(nodes,i,j,velocity_x,dx)
             lodi_x(ix_out)   = eigenvalues(ix_out)*(dPdx + sign_x_out*nodes(i,j,1)*c*dudx)
             
             
             !compute the incoming LODI components in the x-direction
             L_domain_x     = x_map(size(x_map,1))-x_map(1)
             mach_ux_infty  = p_model%get_mach_ux_infty()

             relaxation_u_x = get_relaxation_normal_velocity(L_domain_x,mach_ux_infty,side_x)
             relaxation_v_x = get_relaxation_trans_velocity(L_domain_x, mach_local)
             relaxation_T_x = get_relaxation_temperature(L_domain_x, mach_local)

             lodi_x(1)      = relaxation_v_x*(v-v_set) -
     $                        0.5d0*(1.0d0-relaxation_lodiT)*sign_y_out/
     $                        (nodes(i,j,1)*c)*lodi_y(iy_out)
                            
             lodi_x(2)      = relaxation_T_x*(temp-T_set) -
     $                        (1.0d0-relaxation_lodiT)*lodi_y(2)
                            
             lodi_x(ix_in)  = relaxation_u_x*(u-u_set) -
     $                        (1.0d0-relaxation_lodiT)*(
     $                            0.5d0*lodi_y(iy_out) -
     $                            sign_x_out*nodes(i,j,1)*c*lodi_y(1)
     $                        )


          !inflow along y-direction and outflow along x-direction
          else

             !compute the LODI components in the y-direction
             eigenvalues      = p_model%compute_x_eigenvalues(nodes(i,j,:))
             dx               = x_map(2)-x_map(1)
             dmdx             = gradient_x(nodes,i,j,mass_density,dx)
             dudx             = gradient_x(nodes,i,j,velocity_x,dx)
             dvdx             = gradient_x(nodes,i,j,velocity_y,dx)
             dPdx             = gradient_x(nodes,i,j,pressure,dx)
             
             lodi_x(1)        = eigenvalues(1)*dvdx
             lodi_x(2)        = eigenvalues(2)*(c**2*dmdx-dPdx)
             lodi_x(ix_in)    = 0.0d0
             lodi_x(ix_out)   = eigenvalues(ix_out)*(dPdx + sign_x_out*nodes(i,j,1)*c*dudx)
             
             
             !compute the outgoing LODI acoustics component in the x-direction
             eigenvalues      = p_model%compute_y_eigenvalues(nodes(i,j,:))
             dy               = y_map(2)-y_map(1)
             dPdy             = gradient_y(nodes,i,j,pressure,dy)
             dvdy             = gradient_y(nodes,i,j,velocity_y,dy)
             lodi_y(iy_out)   = eigenvalues(iy_out)*(dPdy + sign_y_out*nodes(i,j,1)*c*dvdy)
             
             
             !compute the incoming LODI components in the x-direction
             L_domain_y     = y_map(size(y_map,1))-y_map(1)
             mach_uy_infty  = p_model%get_mach_uy_infty()

             relaxation_u_y = get_relaxation_trans_velocity(L_domain_y, mach_local)
             relaxation_v_y = get_relaxation_normal_velocity(L_domain_y,mach_uy_infty,side_y)
             relaxation_T_y = get_relaxation_temperature(L_domain_y, mach_local)

             lodi_y(1)      = relaxation_u_y*(u-u_set) -
     $                        0.5d0*(1.0d0-relaxation_lodiT)*sign_x_out/
     $                        (nodes(i,j,1)*c)*lodi_x(ix_out)
                            
             lodi_y(2)      = relaxation_T_y*(temp-T_set) -
     $                        (1.0d0-relaxation_lodiT)*lodi_x(2)
                            
             lodi_y(iy_in)  = relaxation_v_y*(v-v_set) -
     $                        (1.0d0-relaxation_lodiT)*(
     $                            0.5d0*lodi_x(ix_out) -
     $                            sign_y_out*nodes(i,j,1)*c*lodi_x(1)
     $                        )

          end if

        end subroutine compute_x_and_y_lodi

      end module lodi_corner_inflow_outflow_class
