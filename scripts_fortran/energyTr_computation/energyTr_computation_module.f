      !---------------------------------------------------------------------------
      !
      ! MODULE: energyTr_computation_module
      !
      !> @author 
      !> Julien L. Desmarais
      !
      ! DESCRIPTION:
      !  useful subroutines to extract the energy flowing outside
      !  the computational domain
      !> @brief
      !  useful subroutines to extract the energy flowing outside
      !  the computational domain
      !
      ! REVISION HISTORY:
      !>@date
      !> 11_08_2015 - initial version - J.L. Desmarais
      !---------------------------------------------------------------------------
      module energyTr_computation_module
      
        use parameters_cst, only :
     $       N,S,E,W

        use parameters_kind, only :
     $     rkind


        implicit none

        private
        public ::
     $       compute_energyTr,
     $       compute_energyTr_across_x_edge,
     $       compute_energyTr_across_y_edge


        contains


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> computation of the energy flowing outside the computational
        !> domain
        !> @brief
        !> computation of the energy flowing outside the computational
        !> domain
        !
        ! REVISION HISTORY:
        ! 11_08_2015 - initial version - J.L. Desmarais
        !
        !>@param[in]     : x_map    : coordinates along the x-direction
        !>@param[in]     : y_map    : coordinates along the y-direction
        !>@param[in]     : nodes    : grid-points [mass,momentum_x,momentum_y,energy]
        !>@param[in,opt] : borders  : [N,S,E,W] to say check to compute the energy for this border
        !>@param[out]    : energyTr : energy transported flowing through the borders
        !---------------------------------------------------------------------------
        function compute_EnergyTr(
     $       x_map,
     $       y_map,
     $       nodes,
     $       borders)
     $       result(energyTr)

          implicit none

          real(rkind), dimension(:)              , intent(in) :: x_map
          real(rkind), dimension(:)              , intent(in) :: y_map
          real(rkind), dimension(:,:,:)          , intent(in) :: nodes
          logical    , dimension(4)    , optional, intent(in) :: borders
          real(rkind)                                         :: energyTr

          logical, dimension(4) :: borders_opt


          if(present(borders)) then
             borders_opt = borders
          else
             borders_opt = [.true.,.true.,.true.,.true.]
          end if

          energyTr = 0.0d0

        end function compute_EnergyTr        



        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> computation of the energy flowing outside an x-edge
        !> @brief
        !> computation of the energy flowing outside an x-edge
        !
        ! REVISION HISTORY:
        ! 11_08_2015 - initial version - J.L. Desmarais
        !
        !>@param[in]     : x_map    : coordinates along the x-direction
        !>@param[in]     : nodes    : grid-points [mass,momentum_x,momentum_y,energy]
        !>@param[in]     : i_min    : min x-index identifying the x-edge location
        !>@param[in]     : i_max    : max x-index identifying the x-edge location
        !>@param[in]     : j        : y-index identifying the x-edge location at [j+1/2]
        !>@param[out]    : energyTr : energy transported flowing through the borders
        !---------------------------------------------------------------------------
        function compute_EnergyTr_across_x_edge(
     $     x_map,
     $     nodes,
     $     i_min,
     $     i_max,
     $     j)
     $     result(energyTr)

          implicit none

          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer                      , intent(in) :: i_min
          integer                      , intent(in) :: i_max
          integer                      , intent(in) :: j
          real(rkind)                               :: energyTr

          integer     :: i
          real(rkind) :: energy_av
          real(rkind) :: velocity_av

          !---------------------------
          ! nodes(i,j,1) : mass
          ! nodes(i,j,2) : momentum_x
          ! nodes(i,j,3) : momentum_y
          ! nodes(i,j,4) : energy
          ! 
          !     +    +    +     j
          ! ------------------- j+1/2
          !     +    +    +     j-1
          !          i
          !---------------------------
          energyTr = 0.0d0

          do i=i_min, i_max

             energy_av   = 0.5d0*(nodes(i,j,4) + nodes (i,j+1,4))
             velocity_av = 0.5d0*(nodes(i,j,3)/nodes(i,j,1) + nodes(i,j+1,3)/nodes(i,j+1,1))

             energyTr = energyTr + 0.5d0*(energy_av*velocity_av)*(x_map(i+1)-x_map(i-1))

          end do


        end function compute_EnergyTr_across_x_edge


        !---------------------------------------------------------------------------  
        !> @author 
        !> Julien L. Desmarais
        !
        ! DESCRIPTION: 
        !> computation of the energy flowing across an y-edge
        !> @brief
        !> computation of the energy flowing across an y-edge
        !
        ! REVISION HISTORY:
        ! 11_08_2015 - initial version - J.L. Desmarais
        !
        !>@param[in]     : y_map    : coordinates along the x-direction
        !>@param[in]     : nodes    : grid-points [mass,momentum_x,momentum_y,energy]
        !>@param[in]     : i        : x-index identifying the y-edge location at [i+1/2]
        !>@param[in]     : j_min    : min y-index identifying the y-edge location
        !>@param[in]     : j_max    : max y-index identifying the y-edge location
        !>@param[out]    : energyTr : energy transported flowing through the borders
        !---------------------------------------------------------------------------
        function compute_EnergyTr_across_y_edge(
     $     y_map,
     $     nodes,
     $     i,
     $     j_min,
     $     j_max)
     $     result(energyTr)

          implicit none

          real(rkind), dimension(:)    , intent(in) :: y_map
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer                      , intent(in) :: i
          integer                      , intent(in) :: j_min
          integer                      , intent(in) :: j_max
          real(rkind)                               :: energyTr

          integer     :: j
          real(rkind) :: energy_av
          real(rkind) :: velocity_av

          !---------------------------
          ! nodes(i,j,1) : mass
          ! nodes(i,j,2) : momentum_x
          ! nodes(i,j,3) : momentum_y
          ! nodes(i,j,4) : energy
          ! 
          !   i+1/2
          !   |  
          !   | +
          !   |
          ! j | +
          !   |
          !   | +
          !   |
          !---------------------------
          energyTr = 0.0d0

          do j=j_min, j_max

             energy_av   = 0.5d0*(nodes(i,j,4) + nodes (i+1,j,4))
             velocity_av = 0.5d0*(nodes(i,j,2)/nodes(i,j,1) + nodes(i+1,j+1,2)/nodes(i+1,j,1))

             energyTr = energyTr + 0.5d0*(energy_av*velocity_av)*(y_map(j+1)-y_map(j-1))

          end do


        end function compute_EnergyTr_across_y_edge


      end module energyTr_computation_module
