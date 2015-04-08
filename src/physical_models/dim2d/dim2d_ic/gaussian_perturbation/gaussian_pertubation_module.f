      module gaussian_perturbation_module

        use parameters_kind, only :
     $     ikind,
     $     rkind

        implicit none


        private
        public :: add_gaussian_perturbation


        contains

        subroutine add_gaussian_perturbation(
     $       x_map,
     $       y_map,
     $       nodes)

          implicit none

          real(rkind), dimension(nx)   , intent(in)    :: x_map
          real(rkind), dimension(ny)   , intent(in)    :: y_map
          real(rkind), dimension(nx,ny), intent(inout) :: nodes


          integer                                :: Nx_s
          integer                                :: Ny_s
          real(rkind), dimension(:), allocatable :: Spx
          real(rkind), dimension(:), allocatable :: Spy
          real(rkind), dimension(nx)             :: P_x
          real(rkind), dimension(ny)             :: P_y
          
          
          ! compute the number of sinusoidial perturbations
          ! added to the profile alogn the x-direction
          Nx_s = nint((x_map(nx-bc_size)-x_map(1+bc_size))/()
          

          ! compute the spectrum intensity for the perturbation
          ! along the x-direction
          allocate(SPx(Nx_s))
          

          ! compute the perturbation along the x-direction
          do i=1, size(x_map,1)
             
             P_x(i) = 

          end do

          
        end subroutine add_gaussian_perturbation


        function 

      end module gaussian_perturbation_module
