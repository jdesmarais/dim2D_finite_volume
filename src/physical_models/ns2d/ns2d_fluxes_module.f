      !> @file
      !> module encapsulating the computation of the fluxes
      !> for the Navier-Stokes equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the computation of the fluxes
      !> for the Navier-Stokes equations
      !
      !> @date
      ! 08_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ns2d_fluxes_module

        use ns2d_parameters   , only : gamma,
     $                                 viscous_r,
     $                                 mach_infty,
     $                                 Pr, epsilon

        use ns2d_prim_module  , only : momentum_x,
     $                                 momentum_y,
     $                                 velocity_x,
     $                                 velocity_y,
     $                                 temperature,
     $                                 qx_inviscid_x_flux,
     $                                 qy_inviscid_y_flux,
     $                                 qxy_transport,
     $                                 energy_inviscid_x_flux,
     $                                 energy_inviscid_y_flux

        use parameters_kind   , only : ikind, rkind

        use sd_operators_class, only : sd_operators
        

        implicit none

        private
        public ::
     $       flux_x_mass_density,
     $       flux_x_momentum_x  ,
     $       flux_x_momentum_y  ,
     $       flux_x_total_energy,
     $       flux_y_mass_density,
     $       flux_y_momentum_x,  
     $       flux_y_momentum_y,  
     $       flux_y_total_energy,
     $       
     $       flux_x_inviscid_momentum_x,
     $       flux_x_inviscid_momentum_y,
     $       flux_x_inviscid_total_energy,
     $       flux_y_inviscid_momentum_x,
     $       flux_y_inviscid_momentum_y,
     $       flux_y_inviscid_total_energy,
     $       
     $       flux_x_viscid_momentum_x,
     $       flux_x_viscid_momentum_y,
     $       flux_x_viscid_total_energy,
     $       flux_y_viscid_momentum_x,
     $       flux_y_viscid_momentum_y,
     $       flux_y_viscid_total_energy

        contains

        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the mass density \f$ \rho \f$
        !> along the x-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !>@param var
        !> \f$ fx_{\rho} \f$ evaluated at [i,j]
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
        !> along the y-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !>@param var
        !> \f$ fy_{\rho} \f$ evaluated at [i,j]
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
        !> \f$ \rho u_x\f$ along the x-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !
        !>@param var
        !> \f$ fx_{\rho u_x} \f$ evaluated at [i,j]
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


          !DEC$ FORCEINLINE RECURSIVE
          var =  flux_x_inviscid_momentum_x(nodes,s,i,j)
     $         - epsilon*flux_x_viscid_momentum_x(nodes,s,i,j,dx,dy)

        end function flux_x_momentum_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the x-axis
        !> \f$ \rho u_x\f$ along the y-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !
        !>@param var
        !> \f$ fy_{\rho u_x} \f$ evaluated at [i,j]
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


          !DEC$ FORCEINLINE RECURSIVE
          var =  flux_y_inviscid_momentum_x(nodes,s,i,j)
     $         - epsilon*flux_y_viscid_momentum_x(nodes,s,i,j,dx,dy)

        end function flux_y_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the y-axis
        !> \f$ \rho u_y\f$ along the x-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{\rho u_y} \f$ evaluated at [i,j]
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


          !DEC$ FORCEINLINE RECURSIVE
          var = flux_x_inviscid_momentum_y(nodes,s,i,j)
     $         -epsilon*flux_x_viscid_momentum_y(nodes,s,i,j,dx,dy)

        end function flux_x_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the momentum density along the y-axis
        !> \f$ \rho u_y\f$ along the y-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{\rho u_y} \f$ evaluated at [i,j]
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

          
          !DEC$ FORCEINLINE RECURSIVE
          var = flux_y_inviscid_momentum_y(nodes,s,i,j)
     $         -epsilon*flux_y_viscid_momentum_y(nodes,s,i,j,dx,dy)

        end function flux_y_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the total energy density
        !> \f$ \rho E \f$ along the x-axis
        !
        !> @date
        !> 08_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{\rho E} \f$ evaluated at [i,j]
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
      

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,energy_inviscid_x_flux)
     $         -epsilon*flux_x_viscid_total_energy(nodes,s,i,j,dx,dy)

        end function flux_x_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the flux for the total energy density
        !> \f$ \rho E \f$ along the x-axis
        !
        !> @date
        !> 09_08_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{\rho E} \f$ evaluated at [i,j]
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


          var = s%g(nodes,i,j,energy_inviscid_y_flux)
     $         -epsilon*flux_y_viscid_total_energy(nodes,s,i,j,dx,dy)

        end function flux_y_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{qx}^I \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_momentum_x(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,qx_inviscid_x_flux)

        end function flux_x_inviscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{qx}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_x_viscid_momentum_x(nodes,s,i,j,dx,dy)
     $     result(var)

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
             var =(2.0d0+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $            viscous_r*s%dfdy(nodes,i,j,velocity_y,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var =(2.0+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $            viscous_r*s%dfdy(nodes,i,j,velocity_y,dy)
             
          end if

        end function flux_x_viscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{qx}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_y_inviscid_momentum_x(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,qxy_transport)

        end function flux_y_inviscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{qx}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
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


          !DEC$ FORCEINLINE RECURSIVE
          var =s%dgdy(nodes,i,j,velocity_x,dy) +
     $         s%dgdx(nodes,i,j,velocity_y,dx)

        end function flux_y_viscid_momentum_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{qy}^I \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_momentum_y(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,qxy_transport)

        end function flux_x_inviscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{qy}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
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


          !DEC$ FORCEINLINE RECURSIVE
          var =s%dfdy(nodes,i,j,velocity_x,dy)+
     $         s%dfdx(nodes,i,j,velocity_y,dx)

        end function flux_x_viscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{qy}^I \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_y_inviscid_momentum_y(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,qy_inviscid_y_flux)

        end function flux_y_inviscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{qy}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_y_viscid_momentum_y(nodes,s,i,j,dx,dy)
     $     result(var)

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
             var =viscous_r*s%dgdx(nodes,i,j,velocity_x,dx)+
     $            (2.0d0+viscous_r)*s%dgdy(nodes,i,j,velocity_y,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var = viscous_r*s%dgdx(nodes,i,j,velocity_x,dx)+
     $            (2.0+viscous_r)*s%dgdy(nodes,i,j,velocity_y,dy)

          end if
             
        end function flux_y_viscid_momentum_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{qy}^I \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_x_inviscid_total_energy(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%f(nodes,i,j,energy_inviscid_x_flux)

        end function flux_x_inviscid_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fx_{\rho E}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_x_viscid_total_energy(nodes,s,i,j,dx,dy)
     $     result(var)

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
             var =((2.0d0+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $            viscous_r*s%dfdy(nodes,i,j,velocity_y,dy))*
     $            s%f(nodes,i,j,velocity_x)
     $            +
     $            (s%dfdy(nodes,i,j,velocity_x,dy)+s%dfdx(nodes,i,j,velocity_y,dx))*
     $            s%f(nodes,i,j,velocity_y)
     $            +
     $            1.0d0/((gamma-1.0d0)*mach_infty**2*Pr)*
     $            s%dfdx(nodes,i,j,temperature,dx)

          else

             var =((2.0+viscous_r)*s%dfdx(nodes,i,j,velocity_x,dx) +
     $            viscous_r*s%dfdy(nodes,i,j,velocity_y,dy))*
     $            s%f(nodes,i,j,velocity_x)
     $            +
     $            (s%dfdy(nodes,i,j,velocity_x,dy)+s%dfdx(nodes,i,j,velocity_y,dx))*
     $            s%f(nodes,i,j,velocity_y)
     $            +
     $            1.0/((gamma-1.0)*mach_infty**2*Pr)*
     $            s%dfdx(nodes,i,j,temperature,dx)

          end if             

        end function flux_x_viscid_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{\rho E}^I \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_y_inviscid_total_energy(nodes,s,i,j)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          class(sd_operators)          , intent(in) :: s
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          !DEC$ FORCEINLINE RECURSIVE
          var = s%g(nodes,i,j,energy_inviscid_y_flux)

        end function flux_y_inviscid_total_energy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the viscid flux for the momentum_x along the x-axis
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
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
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> \f$ fy_{\rho E}^V \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function flux_y_viscid_total_energy(nodes,s,i,j,dx,dy)
     $     result(var)

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
             var =(s%dgdy(nodes,i,j,velocity_x,dy)+
     $             s%dgdx(nodes,i,j,velocity_y,dx)
     $            )*s%g(nodes,i,j,velocity_x)
     $            +
     $            (viscous_r*s%dgdx(nodes,i,j,velocity_x,dx)+
     $            (2.0d0+viscous_r)*s%dgdy(nodes,i,j,velocity_y,dy)
     $            )*s%g(nodes,i,j,velocity_y)
     $            +
     $            1.0d0/((gamma-1.0d0)*mach_infty**2*Pr)*
     $            s%dgdy(nodes,i,j,temperature,dy)

          else

             !DEC$ FORCEINLINE RECURSIVE
             var =(s%dgdy(nodes,i,j,velocity_x,dy)+
     $             s%dgdx(nodes,i,j,velocity_y,dx)
     $            )*s%g(nodes,i,j,velocity_x)
     $            +
     $            (viscous_r*s%dgdx(nodes,i,j,velocity_x,dx)+
     $            (2.0d0+viscous_r)*s%dgdy(nodes,i,j,velocity_y,dy)
     $            )*s%g(nodes,i,j,velocity_y)
     $            +
     $            1.0d0/((gamma-1.0d0)*mach_infty**2*Pr)*
     $            s%dgdy(nodes,i,j,temperature,dy)

          end if
             
        end function flux_y_viscid_total_energy

      end module ns2d_fluxes_module
