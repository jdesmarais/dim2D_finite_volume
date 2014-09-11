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
     $       compute_edge_fluxes,
     $       compute_lodi_terms,
     $       compute_dev_from_flux_x,
     $       compute_dev_from_flux_y,
     $       get_enhanced_lodi

        
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


          function lodi_matrix(nodes) result(var)

            import ne
            import ikind
            import rkind

            real(rkind), dimension(ne)   , intent(in) :: nodes
            real(rkind), dimension(ne,ne)             :: var

          end function lodi_matrix

      
          function dev_from_flux(flux,i,j,ds) result(dev)

            import ne
            import ikind
            import rkind

            real(rkind), dimension(:,:,:), intent(in) :: flux
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            real(rkind)                  , intent(in) :: ds
            real(rkind), dimension(ne)                :: dev

          end function dev_from_flux

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


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the transverse and viscous LODI vectors
        !
        !> @date
        !> 02_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid points data
        !
        !>@param ds
        !> space step
        !
        !>@param i_offset
        !> index along the x-direction identifying the corresponding index
        !> computed in the nodes table: i -> i_offset+i
        !
        !>@param j_offset
        !> index along the y-direction identifying the corresponding index
        !> computed in the nodes table: j -> j_offset+j
        !
        !>@param epsilon
        !> dissipation constant
        !
        !>@param edge_inviscid_flux
        !> inviscid fluxes at the edge of the computational
        !> domain
        !
        !>@param edge_viscid_flux
        !> viscid fluxes at the edge of the computational
        !> domain
        !
        !>@param compute_conservative_lodi_matrix
        !> compute the conservative LODI matrix in the j-direction
        !
        !>@param transverse_lodi
        !> transverse LODI vector at the edge of the computational domain
        !> in the j-direction
        !
        !>@param viscous_lodi
        !> viscous LODI vector at the edge of the computational domain
        !> in the j-direction
        !---------------------------------------------------------------
        subroutine compute_lodi_terms(
     $     nodes,
     $     ds,
     $     i_offset, j_offset,
     $     epsilon,
     $     edge_inviscid_flux,
     $     edge_viscid_flux,
     $     compute_conservative_lodi_matrix,
     $     compute_dev_from_flux,
     $     transverse_lodi,
     $     viscous_lodi)

          implicit none

          real(rkind), dimension(:,:,:), intent(in)  :: nodes
          real(rkind)                  , intent(in)  :: ds
          integer(ikind)               , intent(in)  :: i_offset
          integer(ikind)               , intent(in)  :: j_offset
          real(rkind)                  , intent(in)  :: epsilon
          real(rkind), dimension(:,:,:), intent(in)  :: edge_inviscid_flux
          real(rkind), dimension(:,:,:), intent(in)  :: edge_viscid_flux
          procedure(lodi_matrix)                     :: compute_conservative_lodi_matrix
          procedure(dev_from_flux)                   :: compute_dev_from_flux
          real(rkind), dimension(:,:,:), intent(out) :: transverse_lodi
          real(rkind), dimension(:,:,:), intent(out) :: viscous_lodi

          integer(ikind)                :: i,j
          integer(ikind)                :: i_nodes,j_nodes
          real(rkind), dimension(ne,ne) :: cons_lodi_matrix
          real(rkind), dimension(ne)    :: dev
          

          do j=1, size(transverse_lodi,2)
             j_nodes = j_offset + j

             do i=1, size(transverse_lodi,1)
                i_nodes = i_offset + i

                cons_lodi_matrix = compute_conservative_lodi_matrix(nodes(i_nodes,j_nodes,:))

                dev = compute_dev_from_flux(edge_inviscid_flux,i,j,ds)
                transverse_lodi(i,j,:) = - MATMUL(dev,cons_lodi_matrix)

                dev = compute_dev_from_flux(edge_viscid_flux,i,j,ds)
                viscous_lodi(i,j,:) = epsilon*MATMUL(dev,cons_lodi_matrix)

             end do

          end do

        end subroutine compute_lodi_terms


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivative from the fluxes along
        !> the x-direction
        !
        !> @date
        !> 03_09_2014 - initial version - J.L. Desmarais
        !
        !>@param flux
        !> fluxes along the x-direction
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param ds
        !> space step
        !
        !>@return dev
        !> derivative
        !---------------------------------------------------------------
        function compute_dev_from_flux_x(flux,i,j,ds) result(dev)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: flux
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: ds
          real(rkind), dimension(ne)                :: dev

          integer :: k

          do k=1,ne
             dev(k) = (flux(i+1,j,k)-flux(i,j,k))/ds
          end do

        end function compute_dev_from_flux_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivative from the fluxes along
        !> the x-direction
        !
        !> @date
        !> 03_09_2014 - initial version - J.L. Desmarais
        !
        !>@param flux
        !> fluxes along the x-direction
        !
        !>@param i
        !> index along the x-direction
        !
        !>@param j
        !> index along the y-direction
        !
        !>@param ds
        !> space step
        !
        !>@return dev
        !> derivative
        !---------------------------------------------------------------
        function compute_dev_from_flux_y(flux,i,j,ds) result(dev)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: flux
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: ds
          real(rkind), dimension(ne)                :: dev

          integer :: k

          do k=1,ne
             dev(k) = (flux(i,j+1,k)-flux(i,j,k))/ds
          end do

        end function compute_dev_from_flux_y


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
