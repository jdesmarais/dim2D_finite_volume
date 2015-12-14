      !> @file
      !> module encapsulating subroutines for the computation of
      !> multiphase quantities: saturated water liquid and vapor mass
      !> densities, interface length...
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of
      !> multiphase quantities: saturated water liquid and vapor mass
      !> densities, interface length...
      !
      !> @date
      !> 14_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_state_eq_module

        use dim2d_parameters, only :
     $       we

        use parameters_input, only :
     $       dim2d_lowTemperature

        use parameters_kind , only :
     $       rkind

        implicit none

        private
        public :: 
     $       get_surface_tension,
     $       get_interface_length,
     $       get_mass_density_liquid,
     $       get_mass_density_vapor


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the surface tension at a given temperature
        !> \f[ \sigma = \displaystyle{\frac{2 \sqrt{2}}{3 We} \left(9.69 {\left(1-T\right)}^{\frac{3}{2}} \right)} \f]
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature, \f$ T \f$
        !
        !>@return
        !> surface tension, \f$ \sigma \f$
        !---------------------------------------------------------------
        function get_surface_tension(temperature)
     $       result(surface_tension)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: surface_tension


          if(rkind.eq.8) then
             surface_tension = 2.0d0*sqrt(2.0d0)/(3.0d0*sqrt(We))*(9.691686257306483d0*(sqrt(1.0d0-temperature))**3)
          else
             surface_tension = 2.0*sqrt(2.0)/(3.0*sqrt(We))*(9.691686257306483*(sqrt(1-temperature))**3)
          end if          

        end function get_surface_tension


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the width of the interface at a given temperature
        !> \f[ L_i = \frac{2}{We} (-0.19 + 1.65 \sqrt{1-T}) \f]
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature, \f$T \f$
        !
        !>@return
        !> width of the interface, \f$ L_i \f$
        !---------------------------------------------------------------
        function get_interface_length(temperature)
     $       result(interface_lgh)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: interface_lgh


c$$$          if(dim2d_lowTemperature) then
c$$$           
c$$$             if(rkind.eq.8) then
c$$$
c$$$                interface_lgh=2.0d0/SQRT(we)*(-0.18904846589641416d0 + 1.6488520842641357d0/(SQRT(1.0d0-temperature)))
c$$$
c$$$             else
c$$$
c$$$                interface_lgh=2.0/SQRT(we)*(-0.18904846589641416d0 + 1.6488520842641357d0/(SQRT(1.0-temperature)))
c$$$
c$$$             end if
c$$$
c$$$          else
           
             if(rkind.eq.8) then

                interface_lgh=1.0d0/SQRT(2.0*we)*(1.6309400153322933d0/(SQRT(1.0d0-temperature)))
                 
             else

                interface_lgh=1.0d0/SQRT(2.0*We)*(1.6309400153322933/(SQRT(1.-temperature)))

             end if

c$$$          end if
              
        end function get_interface_length


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get mass density of saturated liquid water
        !> for the Van der Waals equation
        !> \f[ \rho_{\textrm{liq}} = 1.0 + 2.08 \sqrt{1-T} \f]
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature, \f$ T\f$
        !
        !>@return
        !> mass density of saturated liquid
        !---------------------------------------------------------------
        function get_mass_density_liquid(temperature) result(md_liquid)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: md_liquid
          
          if(dim2d_lowTemperature) then

c$$$             if(rkind.eq.8) then
c$$$                md_liquid = 0.9938571064249365d0 + 2.090919049992133d0*SQRT(1.0d0-temperature)
c$$$             else
c$$$                md_liquid = 0.9938571064249365   + 2.090919049992133*SQRT(1.0-temperature)
c$$$             end if

             if(rkind.eq.8) then
                md_liquid = 1.0d0+2.08d0*SQRT(1.0d0-temperature)
             else
                md_liquid = 1.+2.08*SQRT(1.-temperature)
             end if

          else             

             if(rkind.eq.8) then
                md_liquid = 1.0d0+2.06d0*SQRT(1.0d0-temperature)
             else
                md_liquid = 1.+2.06*SQRT(1.-temperature)
             end if

          end if

        end function get_mass_density_liquid


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get mass density of saturated vapor water
        !> for the Van der Waals equation
        !> \f[ \rho_{\textrm{vap}} = 0.98 - 1.82 \sqrt{1-T} \f]
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature, \f$ T\f$
        !
        !>@return
        !> mass density of saturated vapor
        !---------------------------------------------------------------
        function get_mass_density_vapor(temperature) result(md_vapor)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: md_vapor
          
          if(dim2d_lowTemperature) then

c$$$             if(rkind.eq.8) then
c$$$                md_vapor = 0.9846620916424216d0 - 1.8210340140036445d0*SQRT(1.0d0-temperature)
c$$$             else
c$$$                md_vapor = 0.9846620916424216 - 1.8210340140036445*SQRT(1.0-temperature)
c$$$             end if

             if(rkind.eq.8) then
                md_vapor = 1.0d0-1.86d0*SQRT(1.0d0-temperature)
             else
                md_vapor = 1.-1.86*SQRT(1.-temperature)
             end if

          else

             if(rkind.eq.8) then
                md_vapor = 1.0d0-1.91d0*SQRT(1.0d0-temperature)
             else
                md_vapor = 1.-1.91*SQRT(1.-temperature)
             end if

          end if

        end function get_mass_density_vapor

      end module dim2d_state_eq_module
