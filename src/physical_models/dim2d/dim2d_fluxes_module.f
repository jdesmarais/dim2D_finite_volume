      !> @file
      !> module encapsulating the computation of the fluxes
      !> for the Diffuse Interface Model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the computation of the fluxes
      !> for the Diffuse Interface Model
      !
      !> @date
      !> - 09_08_2013 - initial version               - J.L. Desmarais
      !> - 11_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_fluxes_module

        use dim2d_parameters, only :
     $       viscous_r,
     $       re,pr,we,
     $       cv_r

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use sd_operators_class, only :
     $       sd_operators

        use dim2d_prim_module , only :
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       velocity_x,
     $       velocity_y,
     $       classical_pressure,
     $       classical_temperature_eff,
     $       capillarity_temperature_eff,
     $       temperature_eff,
     $       classical_pressure_xwork,
     $       classical_pressure_ywork,
     $       qx_transport_x,
     $       qy_transport_x,
     $       qx_transport_y,
     $       qy_transport_y,
     $       energy_transport_x,
     $       energy_transport_y,
     $       capillarity_pressure,
     $       capillarity_pressure_xwork,
     $       capillarity_pressure_ywork

        implicit none

        private
        public ::
     $       flux_x_mass_density, flux_y_mass_density,
     $       flux_x_momentum_x,   flux_y_momentum_x,
     $       flux_x_momentum_y,   flux_y_momentum_y,
     $       flux_x_total_energy, flux_y_total_energy,
     $       
     $       flux_x_inviscid_momentum_x,    flux_y_inviscid_momentum_x,    
     $       flux_x_viscid_momentum_x,      flux_y_viscid_momentum_x,      
     $       flux_x_capillarity_momentum_x, flux_y_capillarity_momentum_x,
     $       
     $       flux_x_inviscid_momentum_y,    flux_y_inviscid_momentum_y,   
     $       flux_x_viscid_momentum_y,      flux_y_viscid_momentum_y,     
     $       flux_x_capillarity_momentum_y, flux_y_capillarity_momentum_y,
     $       
     $       flux_x_inviscid_total_energy,    flux_y_inviscid_total_energy,
     $       flux_x_viscid_total_energy,      flux_y_viscid_total_energy,
     $       flux_x_capillarity_total_energy, flux_y_capillarity_total_energy


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the mass density \f$ \rho \f$
        !> along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ F^x_{\rho} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_mass_density(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,momentum_x)

        end function flux_x_mass_density


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the mass density \f$ \rho \f$
        !> along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@return
        !> \f$ F^y_{\rho} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_mass_density(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,momentum_y)

        end function flux_y_mass_density


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the x-axis
        !> \f$ \rho u\f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho u} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var


          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var =  flux_x_inviscid_momentum_x(nodes,s,i,j)
     $            -1.0d0/re*flux_x_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $            -1.0d0/we*flux_x_capillarity_momentum_x(nodes,s,i,j,dx,dy)
          else

             !DEC$ FORCEINLINE RECURSIVE
             var =  flux_x_inviscid_momentum_x(nodes,s,i,j)
     $            -1.0/re*flux_x_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $            -1/we*flux_x_capillarity_momentum_x(nodes,s,i,j,dx,dy)

          end if

        end function flux_x_momentum_x
        

        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the x-axis
        !> \f$ \rho u\f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho u} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
          
          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var =  flux_y_inviscid_momentum_x(nodes,s,i,j)
     $            -1.0d0/re*flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $            -1.0d0/we*flux_y_capillarity_momentum_x(nodes,s,i,j,dx,dy)

          else

            !DEC$ FORCEINLINE RECURSIVE
             var =  flux_y_inviscid_momentum_x(nodes,s,i,j)
     $            -1.0/re*flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $            -1.0/we*flux_y_capillarity_momentum_x(nodes,s,i,j,dx,dy)

          end if

        end function flux_y_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the y-axis
        !> \f$ \rho v\f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho v} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_x_inviscid_momentum_y(nodes,s,i,j)
     $         -1.0d0/re*flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $         -1.0d0/we*flux_x_capillarity_momentum_y(nodes,s,i,j,dx,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_x_inviscid_momentum_y(nodes,s,i,j)
     $            -1.0/re*flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $            -1.0/we*flux_x_capillarity_momentum_y(nodes,s,i,j,dx,dy)

          end if

        end function flux_x_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the y-axis
        !> \f$ \rho v\f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho v} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_y_inviscid_momentum_y(nodes,s,i,j)
     $            -1.0d0/re*flux_y_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $            -1.0d0/we*flux_y_capillarity_momentum_y(nodes,s,i,j,dx,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_y_inviscid_momentum_y(nodes,s,i,j)
     $            -1.0/re*flux_y_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $            -1.0/we*flux_y_capillarity_momentum_y(nodes,s,i,j,dx,dy)

          end if

        end function flux_y_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the total energy density
        !> \f$ \rho E \f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho E} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var


          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var= 
     $            flux_x_inviscid_total_energy(nodes,s,i,j)
     $            -1.0d0/re*flux_x_viscid_total_energy(nodes,s,i,j,dx,dy)
     $            -1.0d0/we*flux_x_capillarity_total_energy(nodes,s,i,j,dx,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var= 
     $            flux_x_inviscid_total_energy(nodes,s,i,j)
     $            -1.0/re*flux_x_viscid_total_energy(nodes,s,i,j,dx,dy)
     $            -1.0/we*flux_x_capillarity_total_energy(nodes,s,i,j,dx,dy)

          end if

        end function flux_x_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the total energy density
        !> \f$ \rho E \f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 09_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho E} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var
          
          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_y_inviscid_total_energy(nodes,s,i,j)
     $            -1.0d0/Re*flux_y_viscid_total_energy(nodes,s,i,j,dx,dy)
     $            -1.0d0/We*flux_y_capillarity_total_energy(nodes,s,i,j,dx,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = flux_y_inviscid_total_energy(nodes,s,i,j)
     $            -1.0/Re*flux_y_viscid_total_energy(nodes,s,i,j,dx,dy)
     $            -1.0/We*flux_y_capillarity_total_energy(nodes,s,i,j,dx,dy)

          end if

        end function flux_y_total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^x_{\rho u} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_momentum_x(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var =
     $         s%f(nodes,i,j,qx_transport_x) +
     $         s%f(nodes,i,j,classical_pressure)

        end function flux_x_inviscid_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho u} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var =
     $         (2.0d0+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $         viscous_r*s%dfdy(nodes,i,j,velocity_y,dy)

        end function flux_x_viscid_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho u} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_capillarity_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var =
     $            1.5d0/cv_r*
     $               s%f(nodes,i,j,capillarity_pressure)*(
     $                  (s%dfdx(nodes,i,j,mass_density,dx))**2 +
     $                  (s%dfdy(nodes,i,j,mass_density,dy))**2
     $               )
     $            +
     $            s%f(nodes,i,j,mass_density)*(
     $                   s%d2fdx2(nodes,i,j,mass_density,dx) +
     $                   s%d2fdy2(nodes,i,j,mass_density,dy))
     $            +
     $            0.5d0*(
     $                 -(s%dfdx(nodes,i,j,mass_density,dx))**2
     $                 +(s%dfdy(nodes,i,j,mass_density,dy))**2)
     $            

          else

             !DEC$ FORCEINLINE RECURSIVE
             var =
     $            1.5/cv_r*
     $               s%f(nodes,i,j,capillarity_pressure)*(
     $                  (s%dfdx(nodes,i,j,mass_density,dx))**2 +
     $                  (s%dfdy(nodes,i,j,mass_density,dy))**2
     $               )
     $            +(
     $               s%f(nodes,i,j,mass_density)*(
     $                   s%d2fdx2(nodes,i,j,mass_density,dx) +
     $                   s%d2fdy2(nodes,i,j,mass_density,dy)) +
     $               0.5*(
     $                 -(s%dfdx(nodes,i,j,mass_density,dx))**2
     $                 +(s%dfdy(nodes,i,j,mass_density,dy))**2)
     $            )

          end if             

        end function flux_x_capillarity_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^y_{\rho u} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_inviscid_momentum_x(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,qx_transport_y)

        end function flux_y_inviscid_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho u} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%dgdy(nodes,i,j,velocity_x,dy) +
     $          s%dgdx(nodes,i,j,velocity_y,dx)

        end function flux_y_viscid_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the momentum
        !> density along the x-axis \f$ \rho u\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho u} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_capillarity_momentum_x(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = - s%dgdx(nodes,i,j,mass_density,dx)*
     $            s%dgdy(nodes,i,j,mass_density,dy)

        end function flux_y_capillarity_momentum_x
        

        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^x_{\rho v} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_momentum_y(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,qy_transport_x)

        end function flux_x_inviscid_momentum_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho v} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%dfdy(nodes,i,j,velocity_x,dy) +
     $          s%dfdx(nodes,i,j,velocity_y,dx)

        end function flux_x_viscid_momentum_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the x-axis
        !> at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho v} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_capillarity_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = - s%dfdx(nodes,i,j,mass_density,dx)*
     $            s%dfdy(nodes,i,j,mass_density,dy)

        end function flux_x_capillarity_momentum_y
        
      
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^y_{\rho v} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_inviscid_momentum_y(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,qy_transport_y) +
     $          s%g(nodes,i,j,classical_pressure)

        end function flux_y_inviscid_momentum_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho v} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = viscous_r*s%dgdx(nodes,i,j,velocity_x,dx) +
     $          (2+viscous_r)*s%dgdy(nodes,i,j,velocity_y,dy)

        end function flux_y_viscid_momentum_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the momentum
        !> density along the y-axis \f$ \rho v\f$ along the y-axis
        !> at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho v} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_capillarity_momentum_y(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = 
     $         1.5d0/cv_r*s%g(nodes,i,j,capillarity_pressure)*(
     $            (s%dgdx(nodes,i,j,mass_density,dx))**2 +
     $            (s%dgdy(nodes,i,j,mass_density,dy))**2
     $         ) +
     $         s%g(nodes,i,j,mass_density)*(
     $             s%d2gdx2(nodes,i,j,mass_density,dx) +
     $             s%d2gdy2(nodes,i,j,mass_density,dy)
     $         ) +
     $         0.5d0*(
     $             (s%dgdx(nodes,i,j,mass_density,dx))**2 -
     $             (s%dgdy(nodes,i,j,mass_density,dy))**2
     $         )

        end function flux_y_capillarity_momentum_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the total energy
        !> density \f$ \rho E\f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^x_{\rho E} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_total_energy(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,energy_transport_x) +
     $          s%f(nodes,i,j,classical_pressure_xwork)

        end function flux_x_inviscid_total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the total energy
        !> density \f$ \rho E\f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho E} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_viscid_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          real(rkind) :: ux,uy

          !DEC$ FORCEINLINE RECURSIVE
          ux     = s%f(nodes,i,j,velocity_x)
          !DEC$ FORCEINLINE RECURSIVE
          uy     = s%f(nodes,i,j,velocity_y)
          
          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = ux*( (2+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $                      viscous_r*s%dfdy(nodes,i,j,velocity_y,dy))
     $             +
     $             uy*( s%dfdy(nodes,i,j,velocity_x,dy) + 
     $                  s%dfdx(nodes,i,j,velocity_y,dx))
     $             +
     $             1.0d0/Pr*s%dfdx(nodes,i,j,classical_temperature_eff,dx)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = ux*( (2+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $                      viscous_r*s%dfdy(nodes,i,j,velocity_y,dy))
     $             +
     $             uy*( s%dfdy(nodes,i,j,velocity_x,dy) + 
     $                  s%dfdx(nodes,i,j,velocity_y,dx))
     $             +
     $             1.0/Pr*s%dfdx(nodes,i,j,classical_temperature_eff,dx)

          end if
             
        end function flux_x_viscid_total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the total
        !> energy density \f$ \rho E\f$ along the x-axis at [i-1/2,j]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^x_{\rho E} \f$ evaluated at [i-1/2,j]
        !---------------------------------------------------------------
        function flux_x_capillarity_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          real(rkind) :: ux,uy,duxdx,duydy,drhodx,drhody

          !DEC$ FORCEINLINE RECURSIVE
          ux     = s%f(nodes,i,j,velocity_x)
          !DEC$ FORCEINLINE RECURSIVE
          uy     = s%f(nodes,i,j,velocity_y)
          !DEC$ FORCEINLINE RECURSIVE
          duxdx  = s%dfdx(nodes,i,j,velocity_x,dx)
          !DEC$ FORCEINLINE RECURSIVE
          duydy  = s%dfdy(nodes,i,j,velocity_y,dy)
          !DEC$ FORCEINLINE RECURSIVE
          drhodx = s%dfdx(nodes,i,j,mass_density,dx)
          !DEC$ FORCEINLINE RECURSIVE
          drhody = s%dfdy(nodes,i,j,mass_density,dy)


          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = 
     $           1.5d0/cv_r*
     $           s%f(nodes,i,j,capillarity_pressure_xwork)*(
     $              (drhodx)**2
     $              + (drhody)**2)
     $              
     $         + s%f(nodes,i,j,momentum_x)*(
     $               s%d2fdx2(nodes,i,j,mass_density,dx)+
     $               s%d2fdy2(nodes,i,j,mass_density,dy))
     $           
     $         + 0.5d0*ux*(
     $             - (drhodx)**2
     $             + (drhody)**2)
     $           
     $         - drhodx*drhody*uy
     $           
     $         - s%f(nodes,i,j,mass_density)*drhodx*(
     $               duxdx + duydy)
     $            
     $         - 1.0d0/(Re*Pr)*s%dfdx_nl(
     $               nodes,i,j,capillarity_temperature_eff,dx,dy)
               
          else

             !DEC$ FORCEINLINE RECURSIVE
             var = 
     $           1.5/cv_r*
     $           s%f(nodes,i,j,capillarity_pressure_xwork)*(
     $              (drhodx)**2
     $              + (drhody)**2)
     $              
     $         + s%f(nodes,i,j,momentum_x)*(
     $               s%d2fdx2(nodes,i,j,mass_density,dx)+
     $               s%d2fdy2(nodes,i,j,mass_density,dy))
     $           
     $         + 0.5*ux*(
     $             - (drhodx)**2
     $             + (drhody)**2)
     $           
     $         - drhodx*drhody*uy
     $           
     $         - s%f(nodes,i,j,mass_density)*drhodx*(
     $               duxdx + duydy)
     $            
     $         - 1.0/(Re*Pr)*s%dfdx_nl(
     $               nodes,i,j,capillarity_temperature_eff,dx,dy)

          end if

        end function flux_x_capillarity_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid part of the flux for the total energy
        !> density \f$ \rho E\f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@return
        !> \f$ F^y_{\rho E} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_inviscid_total_energy(nodes,s,i,j)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,energy_transport_y)
     $         +s%g(nodes,i,j,classical_pressure_ywork)

        end function flux_y_inviscid_total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid part of the flux for the total energy
        !> density \f$ \rho E\f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho E} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_viscid_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          real(rkind) :: ux,uy,duxdx,duydy

          !DEC$ FORCEINLINE RECURSIVE
          ux     = s%g(nodes,i,j,velocity_x)
          !DEC$ FORCEINLINE RECURSIVE
          uy     = s%g(nodes,i,j,velocity_y)
          !DEC$ FORCEINLINE RECURSIVE
          duxdx  = s%dgdx(nodes,i,j,velocity_x,dx)
          !DEC$ FORCEINLINE RECURSIVE
          duydy  = s%dgdy(nodes,i,j,velocity_y,dy)


          if(rkind.eq.8) then

             !DEC$ FORCEINLINE RECURSIVE
             var = ux*( s%dgdy(nodes,i,j,velocity_x,dy)
     $                + s%dgdx(nodes,i,j,velocity_y,dx)) +
     $             uy*(         viscous_r*duxdx
     $                + (2.0d0+viscous_r)*duydy ) +
     $             1.0d0/Pr*s%dgdy(nodes,i,j,classical_temperature_eff,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = ux*( s%dgdy(nodes,i,j,velocity_x,dy)
     $                + s%dgdx(nodes,i,j,velocity_y,dx)) +
     $             uy*(       viscous_r*duxdx
     $                + (2.0+viscous_r)*duydy ) +
     $             1.0/Pr*s%dgdy(nodes,i,j,classical_temperature_eff,dy)

          end if

        end function flux_y_viscid_total_energy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the capillarity part of the flux for the total energy
        !> density \f$ \rho E\f$ along the y-axis at [i,j-1/2]
        !
        !> @date
        !> 11_12_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param dx
        !> grid spacing for the x-axis
        !
        !>@param dy
        !> grid spacing for the y-axis
        !
        !>@return
        !> \f$ F^y_{\rho E} \f$ evaluated at [i,j-1/2]
        !---------------------------------------------------------------
        function flux_y_capillarity_total_energy(nodes,s,i,j,dx,dy)
     $       result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          real(rkind) :: ux,uy,duxdx,duydy,drhodx,drhody

          !DEC$ FORCEINLINE RECURSIVE
          ux     = s%g(nodes,i,j,velocity_x)
          !DEC$ FORCEINLINE RECURSIVE
          uy     = s%g(nodes,i,j,velocity_y)
          !DEC$ FORCEINLINE RECURSIVE
          duxdx  = s%dgdx(nodes,i,j,velocity_x,dx)
          !DEC$ FORCEINLINE RECURSIVE
          duydy  = s%dgdy(nodes,i,j,velocity_y,dy)
          !DEC$ FORCEINLINE RECURSIVE
          drhodx = s%dgdx(nodes,i,j,mass_density,dx)
          !DEC$ FORCEINLINE RECURSIVE
          drhody = s%dgdy(nodes,i,j,mass_density,dy)

          if(rkind.eq.8) then

            !DEC$ FORCEINLINE RECURSIVE
            var = 
     $         1.5d0/cv_r*s%g(nodes,i,j,capillarity_pressure_ywork)*(
     $               (s%dgdx(nodes,i,j,mass_density,dx))**2 +
     $               (s%dgdy(nodes,i,j,mass_density,dy))**2)
     $         
     $         + s%g(nodes,i,j,momentum_y)*(
     $               s%d2gdx2(nodes,i,j,mass_density,dx) +
     $               s%d2gdy2(nodes,i,j,mass_density,dy))
     $         
     $         + 0.5d0*uy*(
     $               (drhodx)**2
     $             - (drhody)**2)
     $         
     $         - drhodx*drhody*ux
     $         
     $         - s%g(nodes,i,j,mass_density)*drhody*(
     $               duxdx + duydy)
     $         
     $         -1.0d0/(Re*Pr)*s%dgdy_nl(
     $               nodes,i,j,capillarity_temperature_eff,dx,dy)

         else

            !DEC$ FORCEINLINE RECURSIVE
            var = 
     $         1.5d0/cv_r*s%g(nodes,i,j,capillarity_pressure_ywork)*(
     $               (s%dgdx(nodes,i,j,mass_density,dx))**2 +
     $               (s%dgdy(nodes,i,j,mass_density,dy))**2)
     $         
     $         + s%g(nodes,i,j,momentum_y)*(
     $               s%d2gdx2(nodes,i,j,mass_density,dx) +
     $               s%d2gdy2(nodes,i,j,mass_density,dy))
     $         
     $         + 0.5d0*uy*(
     $               (drhodx)**2
     $             - (drhody)**2)
     $         
     $         - drhodx*drhody*ux
     $         
     $         - s%g(nodes,i,j,mass_density)*drhody*(
     $               duxdx + duydy)
     $         
     $         -1.0d0/(Re*Pr)*s%dgdy_nl(
     $               nodes,i,j,capillarity_temperature_eff,dx,dy)

         end if

        end function flux_y_capillarity_total_energy


      end module dim2d_fluxes_module
