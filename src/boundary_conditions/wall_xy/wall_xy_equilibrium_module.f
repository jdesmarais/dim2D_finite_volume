      !> @file
      !> ghost cells for wall boundary conditions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute the
      !> ghost cells for the wall equilibrium boundary
      !> conditions
      !
      !> @date
      !> 02_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_xy_equilibrium_module

        use check_data_module, only :
     $       is_real_validated

        use dim2d_parameters, only :
     $       cv_r,
     $       Re,Pr,We

        use dim2d_prim_module, only :
     $       mass_density,
     $       velocity_x,
     $       velocity_y,
     $       temperature_eff,
     $       pressure

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid,
     $       get_surface_tension

        use interface_primary, only :
     $       gradient_proc

        use ISO_FORTRAN_ENV, only :
     $       ERROR_UNIT

        use parameters_constant, only :
     $       left,
     $       no_heat_source,
     $       constant_heat_source,
     $       gaussian_heat_source,
     $       sd_interior_type,
     $       sd_L0_type,
     $       sd_L1_type,
     $       sd_R1_type,
     $       sd_R0_type

        use parameters_input, only :
     $       ne,
     $       wall_micro_contact_angle,
     $       wall_heat_source_choice,
     $       wall_maximum_heat_flux,
     $       wall_extra_heat_source_choice,
     $       wall_maximum_extra_heat_flux,
     $       debug_real

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use ridders_method_module, only :
     $       root_fct_abstract,
     $       get_root_ridder_method

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior

        private
        public ::
     $       dmddx,
     $       temperature,
     $       md_average,
     $       temperature_average,
     $       dwallInternalEnergy_dmd,
     $       
     $       wall_x_equilibrium_root_fct,
     $       wall_x_root_fct,
     $       get_wall_x_root_brackets,
     $       get_wall_x_md_ghost_cell,
     $       
     $       get_wall_heat_flux,
     $       get_wall_micro_contact_angle,
     $       
     $       compute_wall_E_ghost_cell,
     $       compute_wall_W_ghost_cell,
     $       compute_wall_N_ghost_cell,
     $       compute_wall_S_ghost_cell,
     $       
     $       flux_x_inviscid_momentum_x,
     $       flux_x_capillarity_momentum_x,
     $       flux_x_viscid_momentum_y,
     $       flux_x_capillarity_momentum_y,
     $       
     $       flux_y_viscid_momentum_x,
     $       flux_y_capillarity_momentum_x,
     $       flux_y_inviscid_momentum_y,
     $       flux_y_capillarity_momentum_y,
     $       
     $       compute_wall_flux_x,
     $       compute_wall_flux_y


        ! sign_normal: -1.0d0 enforces u=0    at boundary
        !               1.0d0 enforces dudn=0 at boundary
        ! sign_trans : -1.0d0 enforces v=0    at boundary
        real(rkind), parameter :: sign_normal = 1.0d0
        real(rkind), parameter :: sign_trans  =-1.0d0


        ! choice of the method to compute the mass
        ! density in the second ghost cell
        integer, parameter :: wall_dmddx_high_order=0
        integer, parameter :: wall_md_reflection   =1
        integer, parameter :: wall_md_copy         =2

        integer, parameter :: wall_md_ghost_cell2_choice = wall_dmddx_high_order


        !> @class wall_x_root_fct
        !> object encapsulating the function computing 
        !> the wall equilibrium mass density ghost value
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !
        !> @param ini
        !> initialization of the attributes determining
        !> the equilibrium function at the wall
        !
        !> @param f
        !> function whose root determine the equilibrium
        !> value of the mass density at the wall
        !------------------------------------------------------------
        type, extends(root_fct_abstract) :: wall_x_root_fct

          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux
          integer                   :: grad_type_x
          integer                   :: grad_type_y

          contains

          procedure, pass :: ini => wall_x_root_ini
          procedure, pass :: f   => wall_x_root_f

        end type wall_x_root_fct


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux at the wall
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@return flux_x
        !> x-flux at the wall
        !--------------------------------------------------------------
        function compute_wall_flux_x(
     $       nodes,
     $       s,
     $       i,j,
     $       t,
     $       x_map,
     $       y_map,
     $       side,
     $       gradient_y_proc)
     $       result(flux_x)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)             , intent(in) :: s
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          real(rkind)                     , intent(in) :: t
          real(rkind)   , dimension(:)    , intent(in) :: x_map
          real(rkind)   , dimension(:)    , intent(in) :: y_map
          logical                         , intent(in) :: side
          procedure(gradient_proc)                     :: gradient_y_proc
          real(rkind)   , dimension(ne)                :: flux_x

          real(rkind) :: dx
          real(rkind) :: dy
          real(rkind) :: wall_extra_heat_flux

          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)

          flux_x(1) = 0.0d0

          flux_x(2) = flux_x_inviscid_momentum_x(t,x_map,y_map,nodes,dx,dy,i,j,side,gradient_y_proc)
     $                -1.0d0/we*flux_x_capillarity_momentum_x(nodes,s,i,j,dx,dy)

          flux_x(3) = -1.0d0/re*flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $                -1.0d0/we*flux_x_capillarity_momentum_y(nodes,s,i,j,dx,dy)

             
          flux_x(4) = -1.0d0/re*(-get_wall_heat_flux(t,x_map(i),y_map(j)))
          if(side.eqv.left) then
             flux_x(4) = - flux_x(4)
          end if

          if(wall_extra_heat_source_choice.ne.no_heat_source) then

             wall_extra_heat_flux = get_wall_extra_heat_flux(t,x_map(i),y_map(j))

             if(side.eqv.left) then
                flux_x(4) = flux_x(4) - wall_extra_heat_flux
             else
                flux_x(4) = flux_x(4) + wall_extra_heat_flux
             end if

          end if

        end function compute_wall_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux at the wall
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@return flux_x
        !> x-flux at the wall
        !--------------------------------------------------------------
        function compute_wall_flux_y(
     $       nodes,
     $       s,
     $       i,j,
     $       t,
     $       x_map,
     $       y_map,
     $       side,
     $       gradient_x_proc)
     $       result(flux_y)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)             , intent(in) :: s
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          real(rkind)                     , intent(in) :: t
          real(rkind)   , dimension(:)    , intent(in) :: x_map
          real(rkind)   , dimension(:)    , intent(in) :: y_map
          logical                         , intent(in) :: side
          procedure(gradient_proc)                     :: gradient_x_proc
          real(rkind)   , dimension(ne)                :: flux_y

          real(rkind) :: dx
          real(rkind) :: dy

          real(rkind) :: wall_extra_heat_flux


          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)


          flux_y(1) = 0.0d0

          flux_y(2) = -1.0d0/re*flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $                -1.0d0/we*flux_y_capillarity_momentum_x(nodes,s,i,j,dx,dy)

          flux_y(3) = flux_y_inviscid_momentum_y(t,x_map,y_map,nodes,dx,dy,i,j,side,gradient_x_proc)
     $                -1.0d0/we*flux_y_capillarity_momentum_y(nodes,s,i,j,dx,dy)

          flux_y(4) = -1.0d0/re*(-get_wall_heat_flux(t,x_map(i),y_map(j)))
          if(side.eqv.left) then
             flux_y(4) = - flux_y(4)
          end if

          if(wall_extra_heat_source_choice.ne.no_heat_source) then

             wall_extra_heat_flux = get_wall_extra_heat_flux(t,x_map(i),y_map(j))

             if(side.eqv.left) then
                flux_y(4) = flux_y(4) - wall_extra_heat_flux
             else
                flux_y(4) = flux_y(4) + wall_extra_heat_flux
             end if

          end if

        end function compute_wall_flux_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the pressure at the wall
        !
        !> @date
        !> 04_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param nodes
        !> governing variables
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i
        !> abscissa in the x-direction
        !
        !>@param j
        !> abscissa in the y-direction
        !
        !>@param side
        !> left/right depending if the pressure is computed
        !> for the W/E layer
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_x_inviscid_momentum_x(
     $     t,x_map,y_map,nodes,
     $     dx,dy,i,j,
     $     side,
     $     gradient_y_proc)
     $     result(var)

          implicit none
          
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_proc)                  :: gradient_y_proc
          real(rkind)                               :: var

          real(rkind) :: Tl
          real(rkind) :: Th
          real(rkind) :: Pl
          real(rkind) :: Ph


          !for the West layer
          if(side.eqv.left) then

             !temperature at (i,j)
             Tl = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $            nodes,i,j,
     $            dx,dy,
     $            gradient_x_interior,
     $            gradient_y_proc)

             !pressure at (i,j)
             Pl = pressure(Tl,nodes(i,j,1))

             !temperature at (i,j-1)
             Th = Tl + dx*Pr*get_wall_heat_flux(t,x_map(i),y_map(j))

             !pressure at (i,j-1)
             Ph = pressure(Th,nodes(i-1,j,1))
             

          !for the East layer
          else
             
             !temperature at (i,j-1)
             Tl = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $            nodes,i-1,j,
     $            dx,dy,
     $            gradient_x_interior,
     $            gradient_y_proc)

             !pressure at (i,j-1)
             Pl = pressure(Tl,nodes(i-1,j,1))

             !temperature at (i,j)
             Th = Tl + dx*Pr*get_wall_heat_flux(t,x_map(i),y_map(j))

             !pressure at (i,j)
             Ph = pressure(Th,nodes(i,j,1))


          end if

          !average pressure at (i,j-1/2)
          var = 0.5d0*(Pl+Ph)

        end function flux_x_inviscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the x-flux at the wall
        !> for the momentum-x component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_x_capillarity_momentum_x(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
                    
          var = 
     $          s%f(nodes,i,j,mass_density)*(
     $            s%d2fdx2(nodes,i,j,mass_density,dx) +
     $            s%d2fdy2(nodes,i,j,mass_density,dy))
     $          +
     $          0.5d0*(
     $            - (s%dfdx(nodes,i,j,mass_density,dx))**2
     $            + (s%dfdy(nodes,i,j,mass_density,dy))**2)

        end function flux_x_capillarity_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the x-flux at the wall
        !> for the momentum-x component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
                    
          var = s%dfdy(nodes,i,j,velocity_x,dy) +
     $          s%dfdx(nodes,i,j,velocity_y,dx)

        end function flux_x_viscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the x-flux at the wall
        !> for the momentum-y component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_x_capillarity_momentum_y(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
                    

          var = - s%dfdx(nodes,i,j,mass_density,dx)*
     $            s%dfdy(nodes,i,j,mass_density,dy)

        end function flux_x_capillarity_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the x-flux at the wall
        !> for the momentum-y component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
                    
          var = s%dgdy(nodes,i,j,velocity_x,dy) +
     $          s%dgdx(nodes,i,j,velocity_y,dx)

        end function flux_y_viscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the y-flux at the wall
        !> for the momentum-x component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_y_capillarity_momentum_x(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
                    

          var = - s%dgdx(nodes,i,j,mass_density,dx)*
     $            s%dgdy(nodes,i,j,mass_density,dy)

        end function flux_y_capillarity_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the y-flux at the wall
        !> for the momentum-x component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_y_inviscid_momentum_y(
     $     t,x_map,y_map,nodes,
     $     dx,dy,i,j,
     $     side,
     $     gradient_x_proc)
     $     result(var)

          implicit none
          
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_proc)                  :: gradient_x_proc
          real(rkind)                               :: var

          real(rkind) :: Tl
          real(rkind) :: Th
          real(rkind) :: Pl
          real(rkind) :: Ph

          !for the South layer
          if(side.eqv.left) then

             !temperature at (i,j)
             Tl = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $            nodes,i,j,
     $            dx,dy,
     $            gradient_x_proc,
     $            gradient_y_interior)

             !pressure at (i,j)
             Pl = pressure(Tl,nodes(i,j,1))

             !temperature at (i,j-1)
             Th = Tl + dy*Pr*get_wall_heat_flux(t,x_map(i),y_map(j))

             !pressure at (i,j-1)
             Ph = pressure(Th,nodes(i,j-1,1))

          !for the North layer
          else
             
             !temperature at (i,j-1)
             Tl = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $            nodes,i,j-1,
     $            dx,dy,
     $            gradient_x_proc,
     $            gradient_y_interior)

             !pressure at (i,j-1)
             Pl = pressure(Tl,nodes(i,j-1,1))

             !temperature at (i,j)
             Th = Tl + dy*Pr*get_wall_heat_flux(t,x_map(i),y_map(j))

             !pressure at (i,j)
             Ph = pressure(Th,nodes(i,j,1))

          end if

          !average pressure at (i,j-1/2)
          var = 0.5d0*(Pl+Ph)

        end function flux_y_inviscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the y-flux at the wall
        !> for the momentum-x component
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@return var
        !> x-flux at the wall
        !--------------------------------------------------------------
        function flux_y_capillarity_momentum_y(nodes,s,i,j,dx,dy)
     $     result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
          
          var = s%g(nodes,i,j,mass_density)*(
     $            s%d2gdx2(nodes,i,j,mass_density,dx) +
     $            s%d2gdy2(nodes,i,j,mass_density,dy)) +
     $         0.5d0*(
     $              (s%dgdx(nodes,i,j,mass_density,dx))**2
     $            - (s%dgdy(nodes,i,j,mass_density,dy))**2)

        end function flux_y_capillarity_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine compute_wall_S_ghost_cell(
     $     i,j,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     T_guess,
     $     md_guess,
     $     grad_type_x)

          implicit none

          integer(ikind)               , intent(in)    :: i
          integer(ikind)               , intent(in)    :: j
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: T_guess
          real(rkind)                  , intent(in)    :: md_guess
          integer                      , intent(in)    :: grad_type_x


          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux


          md_x = [nodes(i,j+2,1), nodes(i,j+1,1)]

          select case(grad_type_x)

            case(sd_interior_type,sd_L1_type,sd_R1_type)
               md_y = [nodes(i-1,j+1,1), nodes(i+1,j+1,1)]

            case(sd_L0_type)
               md_y = [debug_real, nodes(i+1,j+1,1)]

            case(sd_R0_type)
               md_y = [nodes(i-1,j+1,1), debug_real]

            case default
               call error_grad_type(
     $              'wall_xy_equilibrium_module',
     $              'compute_wall_N_ghost_cell',
     $              grad_type_x)
               
          end select

          velocity_x     =-nodes(i,j+1,3)/nodes(i,j+1,1)
          velocity_y     = nodes(i,j+1,2)/nodes(i,j+1,1)
          Ed             = nodes(i,j+1,4)
          
          delta_x        = y_map(2)-y_map(1)
          delta_y        = x_map(2)-x_map(1)
          micro_angle    = get_wall_micro_contact_angle(x_map(i),y_map(j))
          wall_heat_flux = get_wall_heat_flux(t,x_map(i),y_map(j))


          ! compute the mass density of the ghost cell (i,j)
          nodes(i,j,1) = get_wall_x_md_ghost_cell(
     $            T_guess,
     $            md_guess,
     $            md_x,
     $            md_y,
     $            velocity_x,
     $            velocity_y,
     $            Ed,
     $            delta_x,
     $            delta_y,
     $            micro_angle,
     $            wall_heat_flux,
     $            sd_interior_type,
     $            grad_type_x)
          
          ! compute the mass density of the ghost cell (i,j-1)
          nodes(i,j-1,1) = get_wall_md_extra_ghost_cell(
     $         [nodes(i,j+2,1),nodes(i,j+1,1),nodes(i,j,1)])

          ! reflect the velocity_x
          nodes(i,j,2) = sign_trans*nodes(i,j+1,2)/nodes(i,j+1,1)*nodes(i,j,1)

          ! reflect the velocity_y
          nodes(i,j,3) = sign_normal*nodes(i,j+1,3)/nodes(i,j+1,1)*nodes(i,j,1)
          
          
        end subroutine compute_wall_S_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine compute_wall_N_ghost_cell(
     $     i,j,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     T_guess,
     $     md_guess,
     $     grad_type_x)

          implicit none

          integer(ikind)               , intent(in)    :: i
          integer(ikind)               , intent(in)    :: j
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: T_guess
          real(rkind)                  , intent(in)    :: md_guess
          integer                      , intent(in)    :: grad_type_x


          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux


          md_x = [nodes(i  ,j-2,1), nodes(i  ,j-1,1)]

          select case(grad_type_x)

            case(sd_interior_type,sd_L1_type,sd_R1_type)
               md_y = [nodes(i-1,j-1,1), nodes(i+1,j-1,1)]

            case(sd_L0_type)
               md_y = [debug_real, nodes(i+1,j-1,1)]

            case(sd_R0_type)
               md_y = [nodes(i-1,j-1,1), debug_real]

            case default
               call error_grad_type(
     $              'wall_xy_equilibrium_module',
     $              'compute_wall_N_ghost_cell',
     $              grad_type_x)
               
          end select

          velocity_x     = nodes(i,j-1,3)/nodes(i,j-1,1)
          velocity_y     = nodes(i,j-1,2)/nodes(i,j-1,1)
          Ed             = nodes(i,j-1,4)
          
          delta_x        = y_map(2)-y_map(1)
          delta_y        = x_map(2)-x_map(1)
          micro_angle    = get_wall_micro_contact_angle(x_map(i),y_map(j))
          wall_heat_flux = get_wall_heat_flux(t,x_map(i),y_map(j))


          ! compute the mass density of the ghost cell (i,j)
          nodes(i,j,1) = get_wall_x_md_ghost_cell(
     $            T_guess,
     $            md_guess,
     $            md_x,
     $            md_y,
     $            velocity_x,
     $            velocity_y,
     $            Ed,
     $            delta_x,
     $            delta_y,
     $            micro_angle,
     $            wall_heat_flux,
     $            sd_interior_type,
     $            grad_type_x)

          ! compute the mass density of the ghost cell (i,j+1)
          nodes(i,j+1,1) = get_wall_md_extra_ghost_cell(
     $         [nodes(i,j-2,1),nodes(i,j-1,1),nodes(i,j,1)])

          ! reflect the velocity_x for the ghost cell (i,j)
          nodes(i,j,2) = sign_trans*nodes(i,j-1,2)/nodes(i,j-1,1)*nodes(i,j,1)

          ! reflect the velocity_y for the ghost cell (i,j)
          nodes(i,j,3) = sign_normal*nodes(i,j-1,3)/nodes(i,j-1,1)*nodes(i,j,1)
          
        end subroutine compute_wall_N_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine compute_wall_W_ghost_cell(
     $     i,j,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     T_guess,
     $     md_guess,
     $     grad_type_y)

          implicit none

          integer(ikind)               , intent(in)    :: i
          integer(ikind)               , intent(in)    :: j
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: T_guess
          real(rkind)                  , intent(in)    :: md_guess
          integer                      , intent(in)    :: grad_type_y


          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux


          md_x = [nodes(i+2,j,1), nodes(i+1,j,1)]

          select case(grad_type_y)

            case(sd_interior_type,sd_L1_type,sd_R1_type)
               md_y = [nodes(i+1,j-1,1), nodes(i+1,j+1,1)]

            case(sd_L0_type)
               md_y = [debug_real, nodes(i+1,j+1,1)]

            case(sd_R0_type)
               md_y = [nodes(i+1,j-1,1), debug_real]

            case default
               call error_grad_type(
     $              'wall_xy_equilibrium_module',
     $              'compute_wall_W_ghost_cell',
     $              grad_type_y)
               
          end select

          
          velocity_x     =-nodes(i+1,j,2)/nodes(i+1,j,1)
          velocity_y     = nodes(i+1,j,3)/nodes(i+1,j,1)
          Ed             = nodes(i+1,j,4)
          
          delta_x        = x_map(2)-x_map(1)
          delta_y        = y_map(2)-y_map(1)
          micro_angle    = get_wall_micro_contact_angle(x_map(i),y_map(j))
          wall_heat_flux = get_wall_heat_flux(t,x_map(i),y_map(j))


          ! compute the mass density of the ghost cell (i,j)
          nodes(i,j,1) = get_wall_x_md_ghost_cell(
     $            T_guess,
     $            md_guess,
     $            md_x,
     $            md_y,
     $            velocity_x,
     $            velocity_y,
     $            Ed,
     $            delta_x,
     $            delta_y,
     $            micro_angle,
     $            wall_heat_flux,
     $            sd_interior_type,
     $            grad_type_y)

          ! compute the mass density of the ghost cell (i-1,j)
          nodes(i-1,j,1) = get_wall_md_extra_ghost_cell(
     $         [nodes(i+2,j,1),nodes(i+1,j,1),nodes(i,j,1)])

          ! reflect the velocity_x for the ghost cell (i,j)
          nodes(i,j,2) = sign_normal*nodes(i+1,j,2)/nodes(i+1,j,1)*nodes(i,j,1)

          ! reflect the velocity_y for the ghost cell (i,j)
          nodes(i,j,3) = sign_trans*nodes(i+1,j,3)/nodes(i+1,j,1)*nodes(i,j,1)

          
        end subroutine compute_wall_W_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine compute_wall_E_ghost_cell(
     $     i,j,
     $     t,
     $     x_map,
     $     y_map,
     $     nodes,
     $     T_guess,
     $     md_guess,
     $     grad_type_y)

          implicit none

          integer(ikind)               , intent(in)    :: i
          integer(ikind)               , intent(in)    :: j
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: T_guess
          real(rkind)                  , intent(in)    :: md_guess
          integer                      , intent(in)    :: grad_type_y


          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux


          md_x = [nodes(i-2,j,1), nodes(i-1,j,1)]
          
          select case(grad_type_y)

            case(sd_interior_type,sd_L1_type,sd_R1_type)
               md_y = [nodes(i-1,j-1,1), nodes(i-1,j+1,1)]

            case(sd_L0_type)
               md_y = [debug_real, nodes(i-1,j+1,1)]

            case(sd_R0_type)
               md_y = [nodes(i-1,j-1,1), debug_real]

            case default
               call error_grad_type(
     $              'wall_xy_equilibrium_module',
     $              'compute_wall_E_ghost_cell',
     $              grad_type_y)
               
          end select

          velocity_x     = nodes(i-1,j,2)/nodes(i-1,j,1)
          velocity_y     = nodes(i-1,j,3)/nodes(i-1,j,1)
          Ed             = nodes(i-1,j,4)
          
          delta_x        = x_map(2)-x_map(1)
          delta_y        = y_map(2)-y_map(1)
          micro_angle    = get_wall_micro_contact_angle(x_map(i),y_map(j))
          wall_heat_flux = get_wall_heat_flux(t,x_map(i),y_map(j))


          ! compute the mass density of the ghost cell (i,j)
          nodes(i,j,1) = get_wall_x_md_ghost_cell(
     $            T_guess,
     $            md_guess,
     $            md_x,
     $            md_y,
     $            velocity_x,
     $            velocity_y,
     $            Ed,
     $            delta_x,
     $            delta_y,
     $            micro_angle,
     $            wall_heat_flux,
     $            sd_interior_type,
     $            grad_type_y)

          ! compute the mass density of the ghost cell (i+1,j)
          nodes(i+1,j,1) = get_wall_md_extra_ghost_cell(
     $         [nodes(i-2,j,1),nodes(i-1,j,1),nodes(i,j,1)])

          ! reflect the velocity_x for the ghost cell (i,j)
          nodes(i,j,2) = sign_normal*nodes(i-1,j,2)/nodes(i-1,j,1)*nodes(i,j,1)

          ! reflect the velocity_y for the ghost cell (i,j)
          nodes(i,j,3) = sign_trans*nodes(i-1,j,3)/nodes(i-1,j,1)*nodes(i,j,1)

          
        end subroutine compute_wall_E_ghost_cell


        function get_wall_md_extra_ghost_cell(wall_md)
     $     result(wall_md_ghost_cell)

          implicit none

          real(rkind), dimension(3), intent(in) :: wall_md
          real(rkind)                           :: wall_md_ghost_cell

          select case(wall_md_ghost_cell2_choice)

            case(wall_dmddx_high_order)
               wall_md_ghost_cell = wall_md(1) + 3.0d0*(-wall_md(2) + wall_md(3))

            case(wall_md_reflection)
               wall_md_ghost_cell = wall_md(1)

            case(wall_md_copy)
               wall_md_ghost_cell = wall_md(3)

            case default
               print '(''wall_xy_equilibrium_module'')'
               print '(''get_wall_md_extra_ghost_cell'')'
               print '(''choice for md_ghost_cell2 not recognized'')'
               stop ''
          end select

        end function get_wall_md_extra_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the micro contact angle at the wall
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param x
        !> abscissa in $x$-direction
        !
        !>@param y
        !> abscissa in $y$-direction
        !
        !>@param wall_micro_contact_angle
        !> static micro contact angle at the wall
        !--------------------------------------------------------------
        function get_wall_micro_contact_angle(x,y)
     $     result(contact_angle)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: contact_angle

          real(rkind) :: s

          contact_angle = wall_micro_contact_angle

          s = x+y

        end function get_wall_micro_contact_angle


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the heat flux at the wall
        !
        !> @date
        !> 03_06_2015 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x
        !> abscissa in $x$-direction
        !
        !>@param y
        !> abscissa in $y$-direction
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        function get_wall_heat_flux(t,x,y)
     $     result(wall_heat_flux)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: wall_heat_flux

          real(rkind) :: s
          real(rkind) :: x_center
          real(rkind) :: variance

          select case(wall_heat_source_choice)

            case(no_heat_source)
               wall_heat_flux = 0.0d0
               
            case(constant_heat_source)
               wall_heat_flux = wall_maximum_heat_flux

            case(gaussian_heat_source)
               x_center = 0.0d0
               variance = 0.2d0
               wall_heat_flux = wall_maximum_heat_flux*Exp(-0.5d0*((x-x_center)/variance)**2)

            case default
               print '(''wall_xy_equilibrium_module'')'
               print '(''wall_heat_source_choice not recognized'')'
               print '(''wall_heat_source_choice: '',I1)', wall_heat_source_choice
               stop ''

          end select               

          s = (t+x+y)

        end function get_wall_heat_flux


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the heat flux at the wall for
        !> the extra source term
        !
        !> @date
        !> 05_06_2016 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x
        !> abscissa in $x$-direction
        !
        !>@param y
        !> abscissa in $y$-direction
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        function get_wall_extra_heat_flux(t,x,y)
     $     result(wall_extra_heat_flux)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: wall_extra_heat_flux

          real(rkind) :: s
          real(rkind) :: x_center
          real(rkind) :: variance

          
          select case(wall_extra_heat_source_choice)

            case(no_heat_source)
               wall_extra_heat_flux = 0.0d0

            case(constant_heat_source)
               wall_extra_heat_flux = wall_maximum_extra_heat_flux

            case(gaussian_heat_source)
               x_center = 0.0d0
               variance = 0.02d0
               wall_extra_heat_flux = wall_maximum_extra_heat_flux*Exp(-0.5d0*((x-x_center)/variance)**2)

            case default
               print '(''wall_xy_equilibrium_module'')'
               print '(''wall_extra_heat_source_choice not recognized'')'
               print '(''wall_extra_heat_source_choice: '',I1)', wall_extra_heat_source_choice
               stop ''

          end select               

          s = (t+x+y)

        end function get_wall_extra_heat_flux

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        function get_wall_x_md_ghost_cell(
     $     T_guess,
     $     md_guess,
     $     md_x,
     $     md_y,
     $     velocity_x,
     $     velocity_y,
     $     Ed,
     $     delta_x,
     $     delta_y,
     $     micro_angle,
     $     wall_heat_flux,
     $     grad_type_x,
     $     grad_type_y)
     $     result(md)

          implicit none

          real(rkind)              , intent(in) :: T_guess
          real(rkind)              , intent(in) :: md_guess
          real(rkind), dimension(2), intent(in) :: md_x          
          real(rkind), dimension(2), intent(in) :: md_y          
          real(rkind)              , intent(in) :: velocity_x    
          real(rkind)              , intent(in) :: velocity_y    
          real(rkind)              , intent(in) :: Ed            
          real(rkind)              , intent(in) :: delta_x       
          real(rkind)              , intent(in) :: delta_y       
          real(rkind)              , intent(in) :: micro_angle   
          real(rkind)              , intent(in) :: wall_heat_flux
          integer                  , intent(in) :: grad_type_x
          integer                  , intent(in) :: grad_type_y
          real(rkind)                           :: md


          type(wall_x_root_fct)     :: wall_x_root_fct_used
          real(rkind), dimension(2) :: md_brackets


          ! initialize the parameters configurating 
          ! the wall equilibrium function whose root
          ! determines the mass density in the ghost
          ! cell
          call wall_x_root_fct_used%ini(
     $         md_x,
     $         md_y,
     $         velocity_x,
     $         velocity_y,
     $         Ed,
     $         delta_x,
     $         delta_y,
     $         micro_angle,
     $         wall_heat_flux,
     $         grad_type_x,
     $         grad_type_y)


          ! determine the brackets for the root
          ! solution of the wall equilibrium
          ! function
          md_brackets = get_wall_x_root_brackets(
     $         wall_x_root_fct_used,
     $         T_guess,
     $         md_guess)


          ! determine the root of the wall
          ! equilibrium function using
          ! Ridder's method
          md = get_root_ridder_method(
     $         wall_x_root_fct_used,
     $         md_brackets(1),
     $         md_brackets(2),
     $         1.0d0*1e-12)

        end function get_wall_x_md_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the brackets around the root
        !> of the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param wall_x_root_fct_used
        !> wall equilibirum function
        !
        !>@param T_guess
        !> guess for the temperature at the wall
        !
        !>@param md_guess
        !> guess for the mass density at the wall
        !
        !>@param md_brackets
        !> brackets for the root
        !--------------------------------------------------------------
        function get_wall_x_root_brackets(
     $     wall_x_root_fct_used,
     $     T_guess,
     $     md_guess)
     $     result(md_brackets)

          implicit none

          class(root_fct_abstract) , intent(in) :: wall_x_root_fct_used
          real(rkind)              , intent(in) :: T_guess
          real(rkind)              , intent(in) :: md_guess
          real(rkind), dimension(2)             :: md_brackets

          
          real(rkind) :: xl
          real(rkind) :: xh
          real(rkind) :: fl
          real(rkind) :: fh
          real(rkind) :: fmd

          real(rkind) :: dmd_l
          real(rkind) :: dmd_r
          integer     :: j

          integer, parameter :: MAXIT = 99


          ! use as first potential brackets the
          ! vapor and liquid mass densities at the
          ! temperature T_guess
          xl = get_mass_density_vapor(T_guess)
          xh = get_mass_density_liquid(T_guess)
          
          fl = wall_x_root_fct_used%f(xl)
          fh = wall_x_root_fct_used%f(xh)


          ! if the brackets proposed do not lead to
          ! opposite signs of the wall equilibrium 
          ! function, new brackets should be found
          if( have_same_sign(fl,fh) ) then

             dmd_l = (md_guess)/real(MAXIT+1)
             dmd_r = (3.0d0-md_guess)/real(MAXIT+1)
             fmd   = wall_x_root_fct_used%f(md_guess)

             do j=1, MAXIT

                xl = md_guess - j*dmd_l
                xh = md_guess + j*dmd_r

                fl = wall_x_root_fct_used%f(xl)
                fh = wall_x_root_fct_used%f(xh)

                if(.not.have_same_sign(fl,fmd)) then

                   md_brackets(1) = xl
                   md_brackets(2) = md_guess
                   return

                else if (.not.have_same_sign(fh,fmd)) then

                   md_brackets(1) = md_guess
                   md_brackets(2) = xh
                   return

                end if

             end do

             print '(''wall_xy_equilibrium_module'')'
             print '(''get_wall_x_root_brackets'')'
             print '(''exceeded maximum iterations'')'

             print '(''T_guess: '',F10.4)', T_guess
             print '(''md_guess: '',F10.4)', md_guess

             xl = get_mass_density_vapor(T_guess)
             xh = get_mass_density_liquid(T_guess)
          
             fl = wall_x_root_fct_used%f(xl)
             fh = wall_x_root_fct_used%f(xh)

             print *, '(xl,fl) f(md_vap(T_guess)):', xl,fl
             print *, '(xh,fh) f(md_liq(T_guess)):', xh,fh

             stop ''

          else

             md_brackets(1) = xl
             md_brackets(2) = xh

          end if

        end function get_wall_x_root_brackets


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether two reals have a same sign
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param x1
        !> first real tested
        !>
        !>@param x2
        !> second real tested
        !>
        !>@return same_sign
        !> check whether x1 and x2 have same sign
        !--------------------------------------------------------------
        function have_same_sign(x1,x2) result(same_sign)

          implicit none

          real(rkind), intent(in) :: x1
          real(rkind), intent(in) :: x2
          logical                 :: same_sign

          same_sign = (sign(1.0d0,x1)*sign(1.0d0,x2)).gt.0

        end function have_same_sign


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the additional parameters
        !> determining the wall equilibrium
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> additional parameters to determine the wall
        !> equilibrium
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine wall_x_root_ini(
     $     this,
     $     md_x,
     $     md_y,
     $     velocity_x,
     $     velocity_y,
     $     Ed,
     $     delta_x,
     $     delta_y,
     $     micro_angle,
     $     wall_heat_flux,
     $     grad_type_x,
     $     grad_type_y)

          implicit none

          class(wall_x_root_fct)   , intent(inout) :: this
          real(rkind), dimension(2), intent(in)    :: md_x          
          real(rkind), dimension(2), intent(in)    :: md_y          
          real(rkind)              , intent(in)    :: velocity_x    
          real(rkind)              , intent(in)    :: velocity_y    
          real(rkind)              , intent(in)    :: Ed            
          real(rkind)              , intent(in)    :: delta_x       
          real(rkind)              , intent(in)    :: delta_y       
          real(rkind)              , intent(in)    :: micro_angle   
          real(rkind)              , intent(in)    :: wall_heat_flux
          integer                  , intent(in)    :: grad_type_x
          integer                  , intent(in)    :: grad_type_y

          
          this%md_x           = md_x          
          this%md_y           = md_y          
          this%velocity_x     = velocity_x    
          this%velocity_y     = velocity_y    
          this%Ed             = Ed            
          this%delta_x        = delta_x       
          this%delta_y        = delta_y       
          this%micro_angle    = micro_angle   
          this%wall_heat_flux = wall_heat_flux
          this%grad_type_x    = grad_type_x
          this%grad_type_y    = grad_type_y

        end subroutine wall_x_root_ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function whose root determines the
        !> equilibrium at the wall
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> additional parameters needed to evaluate the
        !> wall equilibrium
        !
        !>@param x
        !> parameter for the non-linear equation
        !
        !>@param fx
        !> evaluation of the non-linear equation at x
        !--------------------------------------------------------------
        function wall_x_root_f(this,x) result(fx)

          implicit none

          class(wall_x_root_fct), intent(in) :: this
          real(rkind)           , intent(in) :: x
          real(rkind)                        :: fx

          fx = wall_x_equilibrium_root_fct(
     $         this%md_x(2),
     $         [this%md_x(1),x],
     $         this%md_y,
     $         this%velocity_x,
     $         this%velocity_y,
     $         this%Ed,
     $         this%delta_x,
     $         this%delta_y,
     $         this%micro_angle,
     $         this%wall_heat_flux,
     $         this%grad_type_x,
     $         this%grad_type_y)

        end function wall_x_root_f


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the wall equilibrium function controlling
        !> the value of the mass density for the ghost cell
        !> \f$ \frac{1}{We} \frac{\partial \rho}{\partial x} 
        !> + \frac{d E^S_{\textrm{wall}}(T,\rho)}{d \rho} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density at {i,j}
        !>
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !>
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !>
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !>
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !>
        !>@param Ed
        !> total energy density at {i,j}
        !>
        !>@param delta_x
        !> space step along the x-direction
        !>
        !>@param delta_y
        !> space step along the y-direction
        !>
        !>@param micro_angle
        !> micro contact angle
        !>
        !>@param wall_heat_flux
        !> heat flux at the wall
        !
        !>@param equilibrium_fct
        !> function controlling the equilibrium at the wall
        !--------------------------------------------------------------
        function wall_x_equilibrium_root_fct(
     $       md,
     $       md_x,
     $       md_y,
     $       velocity_x,
     $       velocity_y,
     $       Ed,
     $       delta_x,
     $       delta_y,
     $       micro_angle,
     $       wall_heat_flux,
     $       grad_type_x,
     $       grad_type_y)
     $       result(equilibrium_fct)

          implicit none

          real(rkind)              , intent(in) :: md
          real(rkind), dimension(2), intent(in) :: md_x
          real(rkind), dimension(2), intent(in) :: md_y
          real(rkind)              , intent(in) :: velocity_x
          real(rkind)              , intent(in) :: velocity_y
          real(rkind)              , intent(in) :: Ed
          real(rkind)              , intent(in) :: delta_x
          real(rkind)              , intent(in) :: delta_y
          real(rkind)              , intent(in) :: micro_angle
          real(rkind)              , intent(in) :: wall_heat_flux
          integer                  , intent(in) :: grad_type_x
          integer                  , intent(in) :: grad_type_y
          real(rkind)                           :: equilibrium_fct

          real(rkind) :: dmddx_half
          real(rkind) :: md_av
          real(rkind) :: md_grad_squared
          real(rkind) :: temperature1
          real(rkind) :: T_average

          dmddx_half = (-md+md_x(2))/(delta_x)

          md_av = md_average(md,md_x(2))

          md_grad_squared =
     $         dmddx(delta_x,[md_x(1),md,md_x(2)],grad_type_x)**2 +
     $         dmddx(delta_y,[md_y(1),md,md_y(2)],grad_type_y)**2

          temperature1 = temperature(
     $         md,
     $         md_grad_squared,
     $         velocity_x,
     $         velocity_y,
     $         Ed)

          T_average = temperature_average(
     $         temperature1,
     $         delta_x,
     $         wall_heat_flux)

          equilibrium_fct =
     $         1.0d0/We*dmddx_half +
     $         dwallInternalEnergy_dmd(md_av,T_average,micro_angle)

        end function wall_x_equilibrium_root_fct


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the derivative of the wall internal energy
        !> per unit surface
        !> \f$ \frac{d E_{\textrm{wall}}}{d \rho} =
        !> \frac{6 (\md - \md_{\textrm{l}})(\md - \md_{\textrm{v}})}
        !> { \left( \rho_{\textrm{l}} - \rho_{\textrm{v}} \right)^3}
        !> \sigma \Cos(\theta_m) \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !>
        !>@param temperature
        !> temperature
        !>
        !>@param micro_angle
        !> microscopic contact angle
        !--------------------------------------------------------------
        function dwallInternalEnergy_dmd(md,temperature,micro_angle)

          implicit none

          real(rkind), intent(in) :: md
          real(rkind), intent(in) :: temperature
          real(rkind), intent(in) :: micro_angle
          real(rkind)             :: dwallInternalEnergy_dmd

          real(rkind) :: md_vap
          real(rkind) :: md_liq
          real(rkind) :: surface_tension
          real(rkind) :: angle

          md_vap = get_mass_density_vapor(temperature)
          md_liq = get_mass_density_liquid(temperature)
          
          surface_tension = get_surface_tension(temperature)

          if(is_real_validated(angle,90.0d0,.false.)) then
             dwallInternalEnergy_dmd = 0.0d0
          else
             
             if(rkind.eq.8) then
                angle = micro_angle*ACOS(-1.0d0)/180.0d0
                dwallInternalEnergy_dmd = 6.0d0*(md-md_liq)*(md-md_vap)/((md_liq-md_vap)**3)*surface_tension*Cos(angle)
             else
                angle = micro_angle*ACOS(-1.0)/180.0
                dwallInternalEnergy_dmd =   6.0*(md-md_liq)*(md-md_vap)/((md_liq-md_vap)**3)*surface_tension*Cos(angle)
             end if
          end if

        end function dwallInternalEnergy_dmd


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass-density at midpoint $i+1/2$:
        !> \f$ \rho_{i+1/2} = \frac{\rho_i+\rho_{i+1}}{2} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param temperature0
        !> temperature at {i}
        !>
        !>@param delta_x
        !> space step
        !>
        !>@param wall heat flux
        !> heat flux at the wall
        !
        !>@param temperature_average
        !> temperature at the wall
        !--------------------------------------------------------------
        function temperature_average(temperature0,delta_x,wall_heat_flux)

          implicit none

          real(rkind), intent(in) :: temperature0
          real(rkind), intent(in) :: delta_x
          real(rkind), intent(in) :: wall_heat_flux
          real(rkind)             :: temperature_average

          
          if(rkind.eq.8) then
             temperature_average =
     $            temperature0 +
     $            0.5d0*delta_x*Pr*wall_heat_flux
          else
             temperature_average =
     $            temperature0 +
     $            0.5*delta_x*Pr*wall_heat_flux
          end if

        end function temperature_average


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass-density at midpoint $i+1/2$:
        !> \f$ \rho_{i+1/2} = \frac{\rho_i+\rho_{i+1}}{2} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md0
        !> mass density at {i}
        !>
        !>@param md1
        !> mass density at {i+1}
        !>
        !>@param md_average
        !> mass density at mid-point
        !--------------------------------------------------------------
        function md_average(md0,md1)

          implicit none

          real(rkind), intent(in) :: md0
          real(rkind), intent(in) :: md1
          real(rkind)             :: md_average

          if(rkind.eq.8) then
             md_average = 0.5d0*(md0+md1)
          else
             md_average = 0.5*(md0+md1)
          end if

        end function md_average


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the temperature at $i$:
        !> \f$ T_i = \frac{3}{8 cv}\left[
        !> \frac{\rho E_i}{\rho_i} - \frac{1}{2} \left(u_i^2 + v_i^2
        !> \right) - \frac{1}{2 We} \frac{(\nabla \md)^2}{\md_i}
        !> - 3 \rho_i \right] \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density at {i,j}
        !>
        !>@param md_grad_squared
        !> mass density gradient at {i,j} squared
        !>
        !>@param velocity_x
        !> velocity along x-direction at {i,j}
        !>
        !>@param velocity-y
        !> velocity along y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !--------------------------------------------------------------
        function temperature(md,md_grad_squared,velocity_x,velocity_y,Ed)
        
          implicit none

          real(rkind), intent(in) :: md
          real(rkind), intent(in) :: md_grad_squared
          real(rkind), intent(in) :: velocity_x
          real(rkind), intent(in) :: velocity_y
          real(rkind), intent(in) :: Ed
          real(rkind)             :: temperature
          
          if(rkind.eq.8) then
             
             temperature = 3.0d0/(8.0d0*cv_r)*(
     $            Ed/md - 0.5d0*(velocity_x**2+velocity_y**2)
     $            - 0.5d0/We*md_grad_squared/md + 3.0d0*md)

          else
             
             temperature = 3.0/(8.0*cv_r)*(
     $            Ed/md - 0.5*(velocity_x**2+velocity_y**2)
     $            - 0.5/We*md_grad_squared/md + 3.0*md)

          end if

        end function temperature
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the gradient of mass density at $i$:
        !> \f$ \frac{d \rho}{d x}\bigg|_i = \frac{-\rho_{i-1} + \rho_{i+1}}
        !> {2 \Delta x} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param delta_x
        !> space step
        !>
        !>@param md0
        !> mass density evaluated at i-1
        !>
        !>@param md2
        !> mass_density evaluated at i+1
        !>
        !>@param drhodx
        !> mass density gradient at $i$
        !--------------------------------------------------------------
        function dmddx(delta_x,md_tab,grad_type)

          implicit none

          real(rkind)              , intent(in) :: delta_x
          real(rkind), dimension(3), intent(in) :: md_tab
          integer                  , intent(in) :: grad_type
          real(rkind)                           :: dmddx

          select case(grad_type)

            case(sd_interior_type,sd_L1_type,sd_R1_type)
               dmddx = (md_tab(3)-md_tab(1))/(2.0d0*delta_x)

            case(sd_L0_type)
               dmddx = (md_tab(3)-md_tab(2))/delta_x

            case(sd_R0_type)
               dmddx = (md_tab(2)-md_tab(1))/delta_x

            case default
               call error_grad_type(
     $              'wall_xy_equilibrium_module',
     $              'dmdx',
     $              grad_type)

          end select

        end function dmddx

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> error printed if grad_type not recognized
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param file_name
        !> name of the file where the exception is caught
        !
        !>@param fct_name
        !> name of the subroutine where the exception is caught
        !
        !>@param grad_type
        !> value of the parameter trigerring the exception
        !--------------------------------------------------------------
        subroutine error_grad_type(
     $     file_name, fct_name, grad_type)

          implicit none

          character(*), intent(in) :: file_name
          character(*), intent(in) :: fct_name
          integer     , intent(in) :: grad_type

          write(ERROR_UNIT, '(A)') file_name
          write(ERROR_UNIT, '(A)') fct_name
          write(ERROR_UNIT, '(A)') 'grad_type not recognized'
          write(ERROR_UNIT, '(''grad_type: '',I2)') grad_type

          stop 'error_grad_type'

        end subroutine error_grad_type

      end module wall_xy_equilibrium_module
