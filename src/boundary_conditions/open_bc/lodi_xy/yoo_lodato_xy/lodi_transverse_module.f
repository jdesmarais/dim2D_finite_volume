      !module for the computation of the transverse and viscous terms
      !needed when deriving the LODI vector
      module lodi_transverse_module

        use sd_operators_class, only :
     $     sd_operators

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: compute_edge_fluxes

        
        abstract interface
          function inviscid_flux(nodes,s,i,j) result(var)

            import sd_operators
            import ikind
            import rkind

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            class(sd_operators)          , intent(in) :: s
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            real(rkind)                               :: var

          end function inviscid_flux

          function viscid_flux(nodes,s,i,j,dx,dy) result(var)

            import sd_operators
            import ikind
            import rkind

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            class(sd_operators)          , intent(in) :: s
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            real(rkind)                  , intent(in) :: dx
            real(rkind)                  , intent(in) :: dy
            real(rkind)                               :: var
          end function viscid_flux

        end interface


        contains

        
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the inviscid and viscid fluxes at the edge of the computational domain
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid points data
        !
        !>@param s
        !> object encapsulating the spatial discretization operators
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i_min
        !> index along the x-direction identifying the first index filled
        !> in the edge_inviscid_flux_x table
        !
        !>@param i_max
        !> index along the x-direction identifying the last index filled
        !> in the edge_inviscid_flux_x table
        !
        !>@param i_offset
        !> index along the x-direction identifying the corresponding index
        !> computed in the nodes table: i -> i_offset+i
        !
        !>@param j_min
        !> index along the y-direction identifying the first index filled
        !> in the edge_inviscid_flux_x table
        !
        !>@param j_max
        !> index along the y-direction identifying the last index filled
        !> in the edge_inviscid_flux_x table
        !
        !>@param j_offset
        !> index along the y-direction identifying the corresponding index
        !> computed in the nodes table: j -> j_offset+j
        !
        !>@param epsilon
        !> dissipation constant
        !
        !>@param flux_inviscid_mass_density
        !> procedure computing the inviscid flux of the mass density
        !
        !>@param flux_inviscid_momentum_x
        !> procedure computing the inviscid flux of the momentum_x
        !
        !>@param flux_inviscid_momentum_y
        !> procedure computing the inviscid flux of the momentum_y
        !
        !>@param flux_inviscid_total_energy
        !> procedure computing the inviscid flux of the total energy
        !
        !>@param flux_viscid_mass_density
        !> procedure computing the viscid flux of the mass density
        !
        !>@param flux_viscid_momentum_x
        !> procedure computing the viscid flux of the momentum_x
        !
        !>@param flux_viscid_momentum_y
        !> procedure computing the viscid flux of the momentum_y
        !
        !>@param flux_viscid_total_energy
        !> procedure computing the viscid flux of the total energy
        !
        !>@param edge_inviscid_flux
        !> inviscid fluxes at the edge of the computational
        !> domain
        !
        !>@param edge_viscid_flux
        !> viscid fluxes at the edge of the computational
        !> domain
        !
        !>@param flux
        !> fluxes at the edge of the computational domain
        !---------------------------------------------------------------
        subroutine compute_edge_fluxes(
     $       nodes,
     $       s,
     $       dx, dy,
     $       i_min, i_max, i_offset,
     $       j_min, j_max, j_offset,
     $       epsilon,
     $       flux_mass_density,
     $       flux_inviscid_momentum_x,
     $       flux_inviscid_momentum_y,
     $       flux_inviscid_total_energy,
     $       flux_viscid_momentum_x,
     $       flux_viscid_momentum_y,
     $       flux_viscid_total_energy,
     $       edge_inviscid_flux,
     $       edge_viscid_flux,
     $       flux)

          implicit none

          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          class(sd_operators)          , intent(in)    :: s
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          integer(ikind)               , intent(in)    :: i_min
          integer(ikind)               , intent(in)    :: i_max
          integer(ikind)               , intent(in)    :: i_offset
          integer(ikind)               , intent(in)    :: j_min
          integer(ikind)               , intent(in)    :: j_max
          integer(ikind)               , intent(in)    :: j_offset
          real(rkind)                  , intent(in)    :: epsilon
                                                       
          procedure(inviscid_flux)                     :: flux_mass_density
                                                    
          procedure(inviscid_flux)                     :: flux_inviscid_momentum_x
          procedure(inviscid_flux)                     :: flux_inviscid_momentum_y
          procedure(inviscid_flux)                     :: flux_inviscid_total_energy
                                                       
          procedure(viscid_flux)                       :: flux_viscid_momentum_x
          procedure(viscid_flux)                       :: flux_viscid_momentum_y
          procedure(viscid_flux)                       :: flux_viscid_total_energy

          real(rkind), dimension(:,:,:), intent(out)   :: edge_inviscid_flux
          real(rkind), dimension(:,:,:), intent(out)   :: edge_viscid_flux
          real(rkind), dimension(:,:,:), intent(inout) :: flux


          integer(ikind) :: i,j
          integer(ikind) :: i_nodes,j_nodes


          do j=j_min, j_max
             j_nodes = j_offset+j

             do i=i_min, i_max
                i_nodes = i_offset + i

                edge_inviscid_flux(i,j,1) = flux_mass_density(nodes,s,i_nodes,j_nodes)
                edge_inviscid_flux(i,j,2) = flux_inviscid_momentum_x(nodes,s,i_nodes,j_nodes)
                edge_inviscid_flux(i,j,3) = flux_inviscid_momentum_y(nodes,s,i_nodes,j_nodes)
                edge_inviscid_flux(i,j,4) = flux_inviscid_total_energy(nodes,s,i_nodes,j_nodes)

                if(rkind.eq.8) then
                   edge_viscid_flux(i,j,1)   = 0.0d0
                else
                   edge_viscid_flux(i,j,1)   = 0.0
                end if
                edge_viscid_flux(i,j,2)   = flux_viscid_momentum_x(nodes,s,i_nodes,j_nodes,dx,dy)
                edge_viscid_flux(i,j,3)   = flux_viscid_momentum_y(nodes,s,i_nodes,j_nodes,dx,dy)
                edge_viscid_flux(i,j,4)   = flux_viscid_total_energy(nodes,s,i_nodes,j_nodes,dx,dy)

                flux(i_nodes,j_nodes,1)   = edge_inviscid_flux(i,j,1)-epsilon*edge_viscid_flux(i,j,1)
                flux(i_nodes,j_nodes,2)   = edge_inviscid_flux(i,j,2)-epsilon*edge_viscid_flux(i,j,2)
                flux(i_nodes,j_nodes,3)   = edge_inviscid_flux(i,j,3)-epsilon*edge_viscid_flux(i,j,3)
                flux(i_nodes,j_nodes,4)   = edge_inviscid_flux(i,j,4)-epsilon*edge_viscid_flux(i,j,4)

             end do
          end do

        end subroutine compute_edge_fluxes

      end module lodi_transverse_module
