      !module for the computation of the transverse and viscous terms
      !needed when deriving the LODI vector
      module lodi_transverse_module

        use sd_operators_class, only :
     $       sd_operators

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public ::
     $       get_enhanced_lodi
c$$$     $       compute_edge_fluxes,
c$$$     $       compute_lodi_terms,
c$$$     $       compute_dev_from_flux_x,
c$$$     $       compute_dev_from_flux_y,


        
c$$$        abstract interface
c$$$          function inviscid_flux(nodes,s,i,j) result(var)
c$$$
c$$$            import sd_operators
c$$$            import ikind
c$$$            import rkind
c$$$
c$$$            real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$            class(sd_operators)          , intent(in) :: s
c$$$            integer(ikind)               , intent(in) :: i
c$$$            integer(ikind)               , intent(in) :: j
c$$$            real(rkind)                               :: var
c$$$
c$$$          end function inviscid_flux
c$$$
c$$$          function viscid_flux(nodes,s,i,j,dx,dy) result(var)
c$$$
c$$$            import sd_operators
c$$$            import ikind
c$$$            import rkind
c$$$
c$$$            real(rkind), dimension(:,:,:), intent(in) :: nodes
c$$$            class(sd_operators)          , intent(in) :: s
c$$$            integer(ikind)               , intent(in) :: i
c$$$            integer(ikind)               , intent(in) :: j
c$$$            real(rkind)                  , intent(in) :: dx
c$$$            real(rkind)                  , intent(in) :: dy
c$$$            real(rkind)                               :: var
c$$$          end function viscid_flux
c$$$
c$$$
c$$$          function lodi_matrix(nodes) result(var)
c$$$
c$$$            import ne
c$$$            import ikind
c$$$            import rkind
c$$$
c$$$            real(rkind), dimension(ne)   , intent(in) :: nodes
c$$$            real(rkind), dimension(ne,ne)             :: var
c$$$
c$$$          end function lodi_matrix
c$$$
c$$$      
c$$$          function dev_from_flux(flux,i,j,ds) result(dev)
c$$$
c$$$            import ne
c$$$            import ikind
c$$$            import rkind
c$$$
c$$$            real(rkind), dimension(:,:,:), intent(in) :: flux
c$$$            integer(ikind)               , intent(in) :: i
c$$$            integer(ikind)               , intent(in) :: j
c$$$            real(rkind)                  , intent(in) :: ds
c$$$            real(rkind), dimension(ne)                :: dev
c$$$
c$$$          end function dev_from_flux
c$$$
c$$$        end interface


        contains

        
c$$$        !> @author 
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the inviscid and viscid fluxes at the edge of the computational domain
c$$$        !
c$$$        !> @date
c$$$        !> 02_09_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param nodes
c$$$        !> array with the grid points data
c$$$        !
c$$$        !>@param s
c$$$        !> object encapsulating the spatial discretization operators
c$$$        !
c$$$        !>@param dx
c$$$        !> space step along the x-direction
c$$$        !
c$$$        !>@param dy
c$$$        !> space step along the y-direction
c$$$        !
c$$$        !>@param i_min
c$$$        !> index along the x-direction identifying the first index filled
c$$$        !> in the edge_inviscid_flux_x table
c$$$        !
c$$$        !>@param i_max
c$$$        !> index along the x-direction identifying the last index filled
c$$$        !> in the edge_inviscid_flux_x table
c$$$        !
c$$$        !>@param i_offset
c$$$        !> index along the x-direction identifying the corresponding index
c$$$        !> computed in the nodes table: i -> i_offset+i
c$$$        !
c$$$        !>@param j_min
c$$$        !> index along the y-direction identifying the first index filled
c$$$        !> in the edge_inviscid_flux_x table
c$$$        !
c$$$        !>@param j_max
c$$$        !> index along the y-direction identifying the last index filled
c$$$        !> in the edge_inviscid_flux_x table
c$$$        !
c$$$        !>@param j_offset
c$$$        !> index along the y-direction identifying the corresponding index
c$$$        !> computed in the nodes table: j -> j_offset+j
c$$$        !
c$$$        !>@param epsilon
c$$$        !> dissipation constant
c$$$        !
c$$$        !>@param flux_inviscid_mass_density
c$$$        !> procedure computing the inviscid flux of the mass density
c$$$        !
c$$$        !>@param flux_inviscid_momentum_x
c$$$        !> procedure computing the inviscid flux of the momentum_x
c$$$        !
c$$$        !>@param flux_inviscid_momentum_y
c$$$        !> procedure computing the inviscid flux of the momentum_y
c$$$        !
c$$$        !>@param flux_inviscid_total_energy
c$$$        !> procedure computing the inviscid flux of the total energy
c$$$        !
c$$$        !>@param flux_viscid_mass_density
c$$$        !> procedure computing the viscid flux of the mass density
c$$$        !
c$$$        !>@param flux_viscid_momentum_x
c$$$        !> procedure computing the viscid flux of the momentum_x
c$$$        !
c$$$        !>@param flux_viscid_momentum_y
c$$$        !> procedure computing the viscid flux of the momentum_y
c$$$        !
c$$$        !>@param flux_viscid_total_energy
c$$$        !> procedure computing the viscid flux of the total energy
c$$$        !
c$$$        !>@param edge_inviscid_flux
c$$$        !> inviscid fluxes at the edge of the computational
c$$$        !> domain
c$$$        !
c$$$        !>@param edge_viscid_flux
c$$$        !> viscid fluxes at the edge of the computational
c$$$        !> domain
c$$$        !
c$$$        !>@param flux
c$$$        !> fluxes at the edge of the computational domain
c$$$        !---------------------------------------------------------------
c$$$        subroutine compute_edge_fluxes(
c$$$     $       nodes,
c$$$     $       s,
c$$$     $       dx, dy,
c$$$     $       i_min, i_max, i_offset,
c$$$     $       j_min, j_max, j_offset,
c$$$     $       epsilon,
c$$$     $       flux_mass_density,
c$$$     $       flux_inviscid_momentum_x,
c$$$     $       flux_inviscid_momentum_y,
c$$$     $       flux_inviscid_total_energy,
c$$$     $       flux_viscid_momentum_x,
c$$$     $       flux_viscid_momentum_y,
c$$$     $       flux_viscid_total_energy,
c$$$     $       edge_inviscid_flux,
c$$$     $       edge_viscid_flux,
c$$$     $       flux)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(in)    :: nodes
c$$$          class(sd_operators)          , intent(in)    :: s
c$$$          real(rkind)                  , intent(in)    :: dx
c$$$          real(rkind)                  , intent(in)    :: dy
c$$$          integer(ikind)               , intent(in)    :: i_min
c$$$          integer(ikind)               , intent(in)    :: i_max
c$$$          integer(ikind)               , intent(in)    :: i_offset
c$$$          integer(ikind)               , intent(in)    :: j_min
c$$$          integer(ikind)               , intent(in)    :: j_max
c$$$          integer(ikind)               , intent(in)    :: j_offset
c$$$          real(rkind)                  , intent(in)    :: epsilon
c$$$                                                       
c$$$          procedure(inviscid_flux)                     :: flux_mass_density
c$$$                                                    
c$$$          procedure(inviscid_flux)                     :: flux_inviscid_momentum_x
c$$$          procedure(inviscid_flux)                     :: flux_inviscid_momentum_y
c$$$          procedure(inviscid_flux)                     :: flux_inviscid_total_energy
c$$$                                                       
c$$$          procedure(viscid_flux)                       :: flux_viscid_momentum_x
c$$$          procedure(viscid_flux)                       :: flux_viscid_momentum_y
c$$$          procedure(viscid_flux)                       :: flux_viscid_total_energy
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(out)   :: edge_inviscid_flux
c$$$          real(rkind), dimension(:,:,:), intent(out)   :: edge_viscid_flux
c$$$          real(rkind), dimension(:,:,:), intent(inout) :: flux
c$$$
c$$$
c$$$          integer(ikind) :: i,j
c$$$          integer(ikind) :: i_nodes,j_nodes
c$$$
c$$$
c$$$          do j=j_min, j_max
c$$$             j_nodes = j_offset+j
c$$$
c$$$             do i=i_min, i_max
c$$$                i_nodes = i_offset + i
c$$$
c$$$                edge_inviscid_flux(i,j,1) = flux_mass_density(nodes,s,i_nodes,j_nodes)
c$$$                edge_inviscid_flux(i,j,2) = flux_inviscid_momentum_x(nodes,s,i_nodes,j_nodes)
c$$$                edge_inviscid_flux(i,j,3) = flux_inviscid_momentum_y(nodes,s,i_nodes,j_nodes)
c$$$                edge_inviscid_flux(i,j,4) = flux_inviscid_total_energy(nodes,s,i_nodes,j_nodes)
c$$$
c$$$                if(rkind.eq.8) then
c$$$                   edge_viscid_flux(i,j,1)   = 0.0d0
c$$$                else
c$$$                   edge_viscid_flux(i,j,1)   = 0.0
c$$$                end if
c$$$                edge_viscid_flux(i,j,2)   = flux_viscid_momentum_x(nodes,s,i_nodes,j_nodes,dx,dy)
c$$$                edge_viscid_flux(i,j,3)   = flux_viscid_momentum_y(nodes,s,i_nodes,j_nodes,dx,dy)
c$$$                edge_viscid_flux(i,j,4)   = flux_viscid_total_energy(nodes,s,i_nodes,j_nodes,dx,dy)
c$$$
c$$$                flux(i_nodes,j_nodes,1)   = edge_inviscid_flux(i,j,1)-epsilon*edge_viscid_flux(i,j,1)
c$$$                flux(i_nodes,j_nodes,2)   = edge_inviscid_flux(i,j,2)-epsilon*edge_viscid_flux(i,j,2)
c$$$                flux(i_nodes,j_nodes,3)   = edge_inviscid_flux(i,j,3)-epsilon*edge_viscid_flux(i,j,3)
c$$$                flux(i_nodes,j_nodes,4)   = edge_inviscid_flux(i,j,4)-epsilon*edge_viscid_flux(i,j,4)
c$$$
c$$$             end do
c$$$          end do
c$$$
c$$$        end subroutine compute_edge_fluxes
c$$$
c$$$
c$$$        !> @author 
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the transverse and viscous LODI vectors
c$$$        !
c$$$        !> @date
c$$$        !> 02_09_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param nodes
c$$$        !> array with the grid points data
c$$$        !
c$$$        !>@param ds
c$$$        !> space step
c$$$        !
c$$$        !>@param i_offset
c$$$        !> index along the x-direction identifying the corresponding index
c$$$        !> computed in the nodes table: i -> i_offset+i
c$$$        !
c$$$        !>@param j_offset
c$$$        !> index along the y-direction identifying the corresponding index
c$$$        !> computed in the nodes table: j -> j_offset+j
c$$$        !
c$$$        !>@param epsilon
c$$$        !> dissipation constant
c$$$        !
c$$$        !>@param edge_inviscid_flux
c$$$        !> inviscid fluxes at the edge of the computational
c$$$        !> domain
c$$$        !
c$$$        !>@param edge_viscid_flux
c$$$        !> viscid fluxes at the edge of the computational
c$$$        !> domain
c$$$        !
c$$$        !>@param compute_conservative_lodi_matrix
c$$$        !> compute the conservative LODI matrix in the j-direction
c$$$        !
c$$$        !>@param transverse_lodi
c$$$        !> transverse LODI vector at the edge of the computational domain
c$$$        !> in the j-direction
c$$$        !
c$$$        !>@param viscous_lodi
c$$$        !> viscous LODI vector at the edge of the computational domain
c$$$        !> in the j-direction
c$$$        !---------------------------------------------------------------
c$$$        subroutine compute_lodi_terms(
c$$$     $     nodes,
c$$$     $     ds,
c$$$     $     i_offset, j_offset,
c$$$     $     epsilon,
c$$$     $     edge_inviscid_flux,
c$$$     $     edge_viscid_flux,
c$$$     $     compute_conservative_lodi_matrix,
c$$$     $     compute_dev_from_flux,
c$$$     $     transverse_lodi,
c$$$     $     viscous_lodi)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(in)  :: nodes
c$$$          real(rkind)                  , intent(in)  :: ds
c$$$          integer(ikind)               , intent(in)  :: i_offset
c$$$          integer(ikind)               , intent(in)  :: j_offset
c$$$          real(rkind)                  , intent(in)  :: epsilon
c$$$          real(rkind), dimension(:,:,:), intent(in)  :: edge_inviscid_flux
c$$$          real(rkind), dimension(:,:,:), intent(in)  :: edge_viscid_flux
c$$$          procedure(lodi_matrix)                     :: compute_conservative_lodi_matrix
c$$$          procedure(dev_from_flux)                   :: compute_dev_from_flux
c$$$          real(rkind), dimension(:,:,:), intent(out) :: transverse_lodi
c$$$          real(rkind), dimension(:,:,:), intent(out) :: viscous_lodi
c$$$
c$$$          integer(ikind)                :: i,j
c$$$          integer(ikind)                :: i_nodes,j_nodes
c$$$          real(rkind), dimension(ne,ne) :: cons_lodi_matrix
c$$$          real(rkind), dimension(ne)    :: dev
c$$$          
c$$$
c$$$          do j=1, size(transverse_lodi,2)
c$$$             j_nodes = j_offset + j
c$$$
c$$$             do i=1, size(transverse_lodi,1)
c$$$                i_nodes = i_offset + i
c$$$
c$$$                cons_lodi_matrix = compute_conservative_lodi_matrix(nodes(i_nodes,j_nodes,:))
c$$$
c$$$                dev = compute_dev_from_flux(edge_inviscid_flux,i,j,ds)
c$$$                transverse_lodi(i,j,:) = - MATMUL(dev,cons_lodi_matrix)
c$$$
c$$$                dev = compute_dev_from_flux(edge_viscid_flux,i,j,ds)
c$$$                viscous_lodi(i,j,:) = epsilon*MATMUL(dev,cons_lodi_matrix)
c$$$
c$$$             end do
c$$$
c$$$          end do
c$$$
c$$$        end subroutine compute_lodi_terms
c$$$
c$$$
c$$$        !> @author 
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the time derivative from the fluxes along
c$$$        !> the x-direction
c$$$        !
c$$$        !> @date
c$$$        !> 03_09_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param flux
c$$$        !> fluxes along the x-direction
c$$$        !
c$$$        !>@param i
c$$$        !> index along the x-direction
c$$$        !
c$$$        !>@param j
c$$$        !> index along the y-direction
c$$$        !
c$$$        !>@param ds
c$$$        !> space step
c$$$        !
c$$$        !>@return dev
c$$$        !> derivative
c$$$        !---------------------------------------------------------------
c$$$        function compute_dev_from_flux_x(flux,i,j,ds) result(dev)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(in) :: flux
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          real(rkind)                  , intent(in) :: ds
c$$$          real(rkind), dimension(ne)                :: dev
c$$$
c$$$          integer :: k
c$$$
c$$$          do k=1,ne
c$$$             dev(k) = (flux(i+1,j,k)-flux(i,j,k))/ds
c$$$          end do
c$$$
c$$$        end function compute_dev_from_flux_x
c$$$
c$$$
c$$$        !> @author 
c$$$        !> Julien L. Desmarais
c$$$        !
c$$$        !> @brief
c$$$        !> compute the time derivative from the fluxes along
c$$$        !> the x-direction
c$$$        !
c$$$        !> @date
c$$$        !> 03_09_2014 - initial version - J.L. Desmarais
c$$$        !
c$$$        !>@param flux
c$$$        !> fluxes along the x-direction
c$$$        !
c$$$        !>@param i
c$$$        !> index along the x-direction
c$$$        !
c$$$        !>@param j
c$$$        !> index along the y-direction
c$$$        !
c$$$        !>@param ds
c$$$        !> space step
c$$$        !
c$$$        !>@return dev
c$$$        !> derivative
c$$$        !---------------------------------------------------------------
c$$$        function compute_dev_from_flux_y(flux,i,j,ds) result(dev)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(in) :: flux
c$$$          integer(ikind)               , intent(in) :: i
c$$$          integer(ikind)               , intent(in) :: j
c$$$          real(rkind)                  , intent(in) :: ds
c$$$          real(rkind), dimension(ne)                :: dev
c$$$
c$$$          integer :: k
c$$$
c$$$          do k=1,ne
c$$$             dev(k) = (flux(i,j+1,k)-flux(i,j,k))/ds
c$$$          end do
c$$$
c$$$        end function compute_dev_from_flux_y


         !> Julien L. Desmarais
         !
         !> @brief
         !> compute the component for the enhancement of the LODI
         !> relations
         !
         !> @date
         !> 04_09_2014 - initial version - J.L. Desmarais
         !
         !>@param relaxation_lodiT
         !> relaxation coefficient for the LODI transverse terms
         !
         !>@param transverse_lodi
         !> transverse lodi component
         !
         !>@param viscous_lodi
         !> viscous lodi component
         !
         !>@return enhanced_lodi
         !> enhanced lodi component
         !---------------------------------------------------------------
         function get_enhanced_lodi(
     $     relaxation_lodiT, transverse_lodi,viscous_lodi)
     $     result(enhanced_lodi)
 
           implicit none
 
           real(rkind), intent(in) :: relaxation_lodiT
           real(rkind), intent(in) :: transverse_lodi
           real(rkind), intent(in) :: viscous_lodi
           real(rkind)             :: enhanced_lodi
 
           if(rkind.eq.8) then
              enhanced_lodi =
     $             (1.0d0-relaxation_lodiT)*transverse_lodi +
     $             viscous_lodi
           else
              enhanced_lodi =
     $             (1.0-relaxation_lodiT)*transverse_lodi +
     $             viscous_lodi
           end if
 
         end function get_enhanced_lodi

      end module lodi_transverse_module
