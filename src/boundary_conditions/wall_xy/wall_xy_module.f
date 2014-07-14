      !> @file
      !> module encapsulating subroutines for the computation
      !> of the prefactors in the wall boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation
      !> of the prefactors in the wall boundary conditions
      !
      !> @date
      ! 24_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_xy_module
      
        use sd_operators_class , only : sd_operators
        use pmodel_eq_class    , only : pmodel_eq
        use dim2d_parameters   , only : rho_c, u_c, length_c, time_c,
     $                                  re, we, viscous_r, cv_r
        use dim2d_prim_module  , only : mass_density,
     $                                  velocity_x, velocity_y,
     $                                  capillarity_pressure
        use parameters_constant, only : vector_x, vector_y
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind,rkind
        use wall_prim_module   , only : wall_pressure

        implicit none


        private
        public :: wall_prefactor,
     $            wall_fx_momentum_x,
     $            wall_fx_momentum_y,
     $            wall_fy_momentum_x,
     $            wall_fy_momentum_y,
     $            wall_heat_flux,
     $            compute_wall_flux_x,
     $            compute_wall_flux_y


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the prefactor in determining the
        !> ghost cells : the prefactor is equal to -1 or +1
        !> depending if the variable type is a vector_x,
        !> vector_y or not
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param prefactor
        !> table containing the prefactor for the reflection
        !--------------------------------------------------------------
        function wall_prefactor(p_model) result(prefactor)

          implicit none

          type(pmodel_eq)        , intent(in) :: p_model
          integer, dimension(ne)             :: prefactor
        
          integer, dimension(ne) :: var_type
          integer :: k

          var_type = p_model%get_var_type()

          do k=1,ne
             if((var_type(k).eq.vector_x).or.(
     $           var_type(k).eq.vector_y)) then
                prefactor(k)=-1
             else
                prefactor(k)= 1
             end if
          end do

        end function wall_prefactor


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the flux for the momentum-x
        !> along the x-direction at a wall
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ {\left(f_{q_x}\right)}_x \f$ evaluated
        !> at $(i+\frac{1}{2},j)$
        !--------------------------------------------------------------
        function wall_fx_momentum_x(nodes,dx,dy,s,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          type(sd_operators)           , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%f(nodes,i,j,wall_pressure)
     $            -1.0d0/re*(2.0d0+viscous_r)*
     $                s%dfdx(nodes,i,j,velocity_x,dx)
     $            -1.0d0/we*(s%dfdy(nodes,i,j,mass_density,dy))**2*(
     $                0.5d0+1.5d0/cv_r*
     $                s%f(nodes,i,j,capillarity_pressure))
     $            -1.0d0/we*s%f(nodes,i,j,mass_density)*(
     $                s%d2fdx2(nodes,i,j,mass_density,dx)+
     $                s%d2fdy2(nodes,i,j,mass_density,dy))

c$$$             print *, 'wall_pressure', s%f(f_used,i,j,wall_pressure)
c$$$             print *, 'viscosity', -1.0d0/re*(2.0d0+viscous_r)*
c$$$     $                s%dfdx(f_used,i,j,velocity_x)
c$$$             print *, 'capillarity_pressure', -1.0d0/we*
c$$$     $                (s%dfdy(f_used,i,j,mass_density))**2*(
c$$$     $                0.5d0+1.5d0/cv_r*
c$$$     $                s%f(f_used,i,j,capillarity_pressure))
c$$$             print *, 'capillarity', -1.0d0/we*
c$$$     $                s%f(f_used,i,j,mass_density)*(
c$$$     $                s%d2fdx2(f_used,i,j,mass_density)+
c$$$     $                s%d2fdy2(f_used,i,j,mass_density))

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%f(nodes,i,j,wall_pressure)
     $            -1.0/re*(2.0+viscous_r)*
     $                s%dfdx(nodes,i,j,velocity_x,dx)
     $            +1.0/we*(s%dfdy(nodes,i,j,mass_density,dy))**2*(
     $                0.5-1.5/cv_r*
     $                s%f(nodes,i,j,capillarity_pressure))
     $            -1.0/we*s%f(nodes,i,j,mass_density)*(
     $                s%d2fdx2(nodes,i,j,mass_density,dx)+
     $                s%d2fdy2(nodes,i,j,mass_density,dy))

          end if

        end function wall_fx_momentum_x
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the flux for the momentum-y
        !> along the x-direction at a wall
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ {\left(f_{q_y}\right)}_x \f$ evaluated
        !> at $(i+\frac{1}{2},j)$
        !--------------------------------------------------------------
        function wall_fx_momentum_y(nodes,dx,s,i,j) result(var)

          implicit none

          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          type(sd_operators)           , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0d0/re*s%dfdx(nodes,i,j,velocity_y,dx)

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0/re*s%dfdx(nodes,i,j,velocity_y,dx)

          end if

        end function wall_fx_momentum_y        
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the flux for the momentum-x
        !> along the y-direction at a wall
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ {\left(f_{q_x}\right)}_y \f$ evaluated
        !> at $(i,j+\frac{1}{2})$
        !--------------------------------------------------------------
        function wall_fy_momentum_x(nodes,dy,s,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dy
          type(sd_operators)           , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0d0/re*s%dgdy(nodes,i,j,velocity_x,dy)

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0/re*s%dgdy(nodes,i,j,velocity_x,dy)

          end if

        end function wall_fy_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the flux for the momentum-x
        !> along the y-direction at a wall
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ {\left(f_{q_y}\right)}_y \f$ evaluated
        !> at $(i,j+\frac{1}{2})$
        !--------------------------------------------------------------
        function wall_fy_momentum_y(nodes,dx,dy,s,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          type(sd_operators)           , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%g(nodes,i,j,wall_pressure)
     $            -1.0d0/re*(2.0d0+viscous_r)*
     $                s%dgdy(nodes,i,j,velocity_y,dy)
     $            -1.0d0/we*(s%dgdx(nodes,i,j,mass_density,dx))**2*(
     $                0.5d0+1.5d0/cv_r*
     $                s%g(nodes,i,j,capillarity_pressure))
     $            -1.0d0/we*s%g(nodes,i,j,mass_density)*(
     $                s%d2gdx2(nodes,i,j,mass_density,dx)+
     $                s%d2gdy2(nodes,i,j,mass_density,dy))

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%g(nodes,i,j,wall_pressure)
     $            -1.0/re*(2.0+viscous_r)*
     $                s%dgdy(nodes,i,j,velocity_y,dy)
     $            -1.0/we*(s%dgdx(nodes,i,j,mass_density,dx))**2*(
     $                0.5+1.5/cv_r*
     $                s%g(nodes,i,j,capillarity_pressure))
     $            -1.0/we*s%g(nodes,i,j,mass_density)*(
     $                s%d2gdx2(nodes,i,j,mass_density,dx)+
     $                s%d2gdy2(nodes,i,j,mass_density,dy))

          end if

        end function wall_fy_momentum_y

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the heat flux at the wall
        !
        !> @date
        !> 24_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> heat flux evaluated at (x,y)
        !--------------------------------------------------------------
        function wall_heat_flux(nodes,i,j) result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var


          real(rkind) :: node_s
          integer(ikind) :: i_s, j_s          
          real(rkind) :: heater_x
          real(rkind) :: heater_y
          real(rkind) :: heater_sigma
          real(rkind) :: heater_power

          node_s = nodes(1,1,1)
          i_s = i
          j_s = j


          !< the heater is localized at the bottom of the
          !> system. It is a square of 1mm by 1mm. In 2D, it
          !> is modelled as a segment of length 1mm.
          heater_x=0.0
          heater_y=10.0
          heater_sigma=(1.0e-3)/length_c
          heater_power=1.0*time_c/(rho_c*u_c**2)

          var = 0.0d0

          !if(rkind.eq.8) then
          !
          !   if(f_used%y_map(j).eq.(heater_y+f_used%dy/2.0d0)) then
          !      var = heater_power/
     $    !           (heater_sigma*sqrt(2.0d0*acos(-1.0d0)))*
     $    !           exp(-(f_used%x_map(i)-heater_x)**2/
     $    !           (2.0d0*heater_sigma**2))
          !   else
          !      var=0.0d0
          !   end if
          !
          !else
          !
          !   if(f_used%y_map(j).eq.(heater_y+f_used%dy/2.0d0)) then
          !      var = heater_power/
     $    !           (heater_sigma*sqrt(2.0*acos(-1.0)))*
     $    !           exp(-(f_used%x_map(i)-heater_x)**2/
     $    !           (2.0*heater_sigma**2))
          !   else
          !      var=0.0
          !   end if
          !   
          !end if

        end function wall_heat_flux

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the heat flux at the wall
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param s
        !> spatial discretisation operators
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param flux_x
        !> modified table for the fluxes along the x-direction
        !--------------------------------------------------------------
        subroutine compute_wall_flux_x(nodes,dx,dy,s,i,flux_x)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s
          integer(ikind)                    , intent(in)    :: i
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          integer(ikind) :: j


          do j=bc_size+1, ny-bc_size
             

             !< no mass entering the system
             flux_x(i,j,1) = 0.0d0
                
             !< b.c. for the momentum along the x-direction
             flux_x(i,j,2) = wall_fx_momentum_x(nodes,dx,dy,s,i-1,j)

             !< b.c. for the momentum along the y-direction
             flux_x(i,j,3) = wall_fx_momentum_y(nodes,dx,s,i-1,j)

             !< constant heat flux entering the system
             flux_x(i,j,4) = wall_heat_flux(nodes,i-1,j)

          end do        

        end subroutine compute_wall_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the heat flux at the wall
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f_used
        !> object encapsulating the conservative variables
        !> and the coordinates
        !>
        !>@param s
        !> spatial discretisation operators
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param flux_y
        !> modified table for the fluxes along the y-direction
        !--------------------------------------------------------------
        subroutine compute_wall_flux_y(nodes,dx,dy,s,j,flux_y)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s
          integer(ikind)                    , intent(in)    :: j
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer(ikind) :: i


          do i=bc_size+1, nx-bc_size

            !< no mass entering the system
            flux_y(i,j,1)= 0.0d0
            
            !< b.c. for the momentum along the x-direction
            flux_y(i,j,2)= wall_fy_momentum_x(nodes,dy,s,i,j-1)

            !< b.c. for the momentum along the y-direction
            flux_y(i,j,3)= wall_fy_momentum_y(nodes,dx,dy,s,i,j-1)

            !< constant heat flux entering the system
            flux_y(i,j,4)= wall_heat_flux(nodes,i,j-1)

          end do

        end subroutine compute_wall_flux_y

      end module wall_xy_module
