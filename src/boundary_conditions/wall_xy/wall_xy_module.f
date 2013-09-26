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
      
        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use dim2d_eq_class     , only : dim2d_eq
        use dim2d_parameters   , only : rho_c, u_c, length_c, time_c,
     $                                  re, we, viscous_r, cv_r
        use dim2d_prim_module  , only : mass_density,
     $                                  velocity_x, velocity_y,
     $                                  capillarity_pressure
        use parameters_constant, only : vector_x, vector_y
        use parameters_input   , only : nx,ny,ne
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

          type(dim2d_eq)        , intent(in) :: p_model
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
        function wall_fx_momentum_x(f_used,s,i,j) result(var)

          implicit none

          class(field)      , intent(in) :: f_used
          type(cg_operators), intent(in) :: s
          integer(ikind)    , intent(in) :: i
          integer(ikind)    , intent(in) :: j
          real(rkind)                    :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%f(f_used,i,j,wall_pressure)
     $            -1.0d0/re*(2.0d0+viscous_r)*
     $                s%dfdx(f_used,i,j,velocity_x)
     $            -1.0d0/we*(s%dfdy(f_used,i,j,mass_density))**2*(
     $                0.5d0+1.5d0/cv_r*
     $                s%f(f_used,i,j,capillarity_pressure))
     $            -1.0d0/we*s%f(f_used,i,j,mass_density)*(
     $                s%d2fdx2(f_used,i,j,mass_density)+
     $                s%d2fdy2(f_used,i,j,mass_density))

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
             var=s%f(f_used,i,j,wall_pressure)
     $            -1.0/re*(2.0+viscous_r)*
     $                s%dfdx(f_used,i,j,velocity_x)
     $            +1.0/we*(s%dfdy(f_used,i,j,mass_density))**2*(
     $                0.5-1.5/cv_r*
     $                s%f(f_used,i,j,capillarity_pressure))
     $            -1.0/we*s%f(f_used,i,j,mass_density)*(
     $                s%d2fdx2(f_used,i,j,mass_density)+
     $                s%d2fdy2(f_used,i,j,mass_density))

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
        function wall_fx_momentum_y(f_used,s,i,j) result(var)

          implicit none

          class(field)      , intent(in) :: f_used
          type(cg_operators), intent(in) :: s
          integer(ikind)    , intent(in) :: i
          integer(ikind)    , intent(in) :: j
          real(rkind)                    :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0d0/re*s%dfdx(f_used,i,j,velocity_y)

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0/re*s%dfdx(f_used,i,j,velocity_y)

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
        function wall_fy_momentum_x(f_used,s,i,j) result(var)

          implicit none

          class(field)      , intent(in) :: f_used
          type(cg_operators), intent(in) :: s
          integer(ikind)    , intent(in) :: i
          integer(ikind)    , intent(in) :: j
          real(rkind)                    :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0d0/re*s%dgdy(f_used,i,j,velocity_x)

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var= -1.0/re*s%dgdy(f_used,i,j,velocity_x)

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
        function wall_fy_momentum_y(f_used,s,i,j) result(var)

          implicit none

          class(field)      , intent(in) :: f_used
          type(cg_operators), intent(in) :: s
          integer(ikind)    , intent(in) :: i
          integer(ikind)    , intent(in) :: j
          real(rkind)                    :: var


          if(rkind.eq.8) then
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%g(f_used,i,j,wall_pressure)
     $            -1.0d0/re*(2.0d0+viscous_r)*
     $                s%dgdy(f_used,i,j,velocity_y)
     $            -1.0d0/we*(s%dgdx(f_used,i,j,mass_density))**2*(
     $                0.5d0+1.5d0/cv_r*
     $                s%g(f_used,i,j,capillarity_pressure))
     $            -1.0d0/we*s%g(f_used,i,j,mass_density)*(
     $                s%d2gdx2(f_used,i,j,mass_density)+
     $                s%d2gdy2(f_used,i,j,mass_density))

          else
             
             !DEC$ FORCEINLINE RECURSIVE
             var=s%g(f_used,i,j,wall_pressure)
     $            -1.0/re*(2.0+viscous_r)*
     $                s%dgdy(f_used,i,j,velocity_y)
     $            -1.0/we*(s%dgdx(f_used,i,j,mass_density))**2*(
     $                0.5+1.5/cv_r*
     $                s%g(f_used,i,j,capillarity_pressure))
     $            -1.0/we*s%g(f_used,i,j,mass_density)*(
     $                s%d2gdx2(f_used,i,j,mass_density)+
     $                s%d2gdy2(f_used,i,j,mass_density))

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
        function wall_heat_flux(f_used,i,j) result(var)

          implicit none
          
          class(field)      , intent(in) :: f_used
          integer(ikind)    , intent(in) :: i
          integer(ikind)    , intent(in) :: j
          real(rkind)                    :: var


          real(rkind) :: heater_x
          real(rkind) :: heater_y
          real(rkind) :: heater_sigma
          real(rkind) :: heater_power


          !< the heater is localized at the bottom of the
          !> system. It is a square of 1mm by 1mm. In 2D, it
          !> is modelled as a segment of length 1mm.
          !> 
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
        subroutine compute_wall_flux_x(f_used,s,i,flux_x)

          implicit none

          class(field)                      , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s
          integer(ikind)                    , intent(in)    :: i
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x

          integer        :: bc_size
          integer(ikind) :: j


          !< get the size of the boundary layer
          bc_size = s%get_bc_size()


          do j=bc_size+1, ny-bc_size
             

             !< no mass entering the system
             flux_x(i,j,1) = 0.0d0
                
             !< b.c. for the momentum along the x-direction
             flux_x(i,j,2) = wall_fx_momentum_x(f_used,s,i-1,j)

             !< b.c. for the momentum along the y-direction
             flux_x(i,j,3) = wall_fx_momentum_y(f_used,s,i-1,j)

             !< constant heat flux entering the system
             flux_x(i,j,4) = wall_heat_flux(f_used,i-1,j)

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
        subroutine compute_wall_flux_y(f_used,s,j,flux_y)

          implicit none

          class(field)                      , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s
          integer(ikind)                    , intent(in)    :: j
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer        :: bc_size
          integer(ikind) :: i


          !< get the size of the boundary layer
          bc_size = s%get_bc_size()


          do i=bc_size+1, nx-bc_size

            !< no mass entering the system
            flux_y(i,j,1)= 0.0d0
            
            !< b.c. for the momentum along the x-direction
            flux_y(i,j,2)= wall_fy_momentum_x(f_used,s,i,j-1)

            !< b.c. for the momentum along the y-direction
            flux_y(i,j,3)= wall_fy_momentum_y(f_used,s,i,j-1)

            !< constant heat flux entering the system
            flux_y(i,j,4)= wall_heat_flux(f_used,i,j-1)

          end do

        end subroutine compute_wall_flux_y

      end module wall_xy_module
