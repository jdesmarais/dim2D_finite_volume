      !> @file
      !> random noise modelling
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> application of a random noise modelled by the approximation
      !> of a stationary Gaussian signal, i.e. a large number of
      !> sinusoidal functions of random phase angles
      !
      !> @date
      !> 09_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module gaussian_perturbation_module

        use parameters_input, only :
     $     nx,ny,bc_size

        use parameters_kind, only :
     $     ikind,
     $     rkind

        implicit none


        private
        public ::
     $       add_gaussian_perturbation,
     $       compute_gaussian_perturbation_smooth,
     $       compute_wave_numbers,
     $       compute_amplitudes,
     $       power_spectrum,
     $       compute_random_phases,
     $       compute_gaussian_intensity,
     $       smooth


        integer    , parameter :: p=8
        real(rkind), parameter :: pi = ACOS(-1.0d0)


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the random noise to the field
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param y_map
        !> map of y-coordinates
        !
        !>@param nodes
        !> array with the grid point data
        !---------------------------------------------------------------
        subroutine add_gaussian_perturbation(
     $       x_map,
     $       y_map,
     $       nodes,
     $       amplitude)

          implicit none

          real(rkind), dimension(nx)   , intent(in)    :: x_map
          real(rkind), dimension(ny)   , intent(in)    :: y_map
          real(rkind), dimension(nx,ny), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: amplitude


          integer, parameter :: seed_x = 1384973980
          integer, parameter :: seed_y = 2179427892


          real(rkind) :: x_min
          real(rkind) :: x_max
          real(rkind) :: y_min
          real(rkind) :: y_max          

          real(rkind), dimension(nx) :: P_x
          real(rkind), dimension(ny) :: P_y

          real(rkind) :: P_x_max
          real(rkind) :: P_y_max
          
          integer     :: i,j
          real(rkind) :: x,y
          real(rkind) :: scaling
          real(rkind) :: perturbation


          ! computation of the border coordinates
          x_min = x_map(bc_size+1)
          x_max = x_map(size(x_map,1)-bc_size)
          y_min = y_map(bc_size+1)
          y_max = y_map(size(y_map,1)-bc_size)


          ! computation of the gaussian perturbation
          ! along the x-direction
          call compute_gaussian_perturbation_smooth(x_min,x_max,x_map, seed_x, P_x, P_x_max)


          ! computation of the gaussian perturbation
          ! along the y-direction
          call compute_gaussian_perturbation_smooth(y_min,y_max,y_map, seed_y, P_y, P_y_max)


          ! add the perturbation to the field
          scaling = amplitude/(P_x_max*P_y_max)

          do j=1+bc_size, size(y_map,1)-bc_size

             y = y_map(j)

             do i=1+bc_size, size(x_map,1)-bc_size

                x = x_map(i)

                perturbation = scaling*P_x(i)*P_y(j)

                nodes(i,j) = nodes(i,j) + perturbation

             end do

          end do

        end subroutine add_gaussian_perturbation


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the amplitude of the gaussian perturbation
        !> over the coordinate map
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param seed_x
        !> random seed used to compute the random phase angles
        !
        !>@param P_x
        !> gaussian perturbation over the coordinates
        !
        !>@param P_x_max
        !> maximum of the gaussian perturbation over the coordinates
        !---------------------------------------------------------------
        subroutine compute_gaussian_perturbation_smooth(x_min,x_max,x_map, seed_x, P_x, P_x_max)

          implicit none

          real(rkind)              , intent(in)  :: x_min
          real(rkind)              , intent(in)  :: x_max
          real(rkind), dimension(:), intent(in)  :: x_map
          integer                  , intent(in)  :: seed_x
          real(rkind), dimension(:), intent(out) :: P_x
          real(rkind)              , intent(out) :: P_x_max


          real(rkind) :: dx
          integer     :: Nx_s

          real(rkind), dimension(:), allocatable :: kx
          real(rkind), dimension(:), allocatable :: Ax
          real(rkind), dimension(:), allocatable :: Phix

          integer :: i


          ! compute the number of sinusoidial perturbations
          ! added to the profile along the x-direction
          dx   = x_map(2) - x_map(1)
          Nx_s = nint((x_max-x_min)/(p*dx))
          

          ! compute the wave numbers for the perturbation
          ! along the x-direction
          allocate(kx(Nx_s))
          call compute_wave_numbers(x_min,x_max,kx)


          ! compute the amplitude for the perturbation
          ! along the x-direction
          allocate(Ax(Nx_s))
          call compute_amplitudes(kx,Ax)


          ! compute the random phases for the pertubation
          ! along the x-direction
          allocate(Phix(Nx_s))
          call compute_random_phases(seed_x,Phix)
          

          ! compute the perturbation along the x-direction
          P_x_max = 0.0d0

          do i=1, size(x_map,1)
             
             P_x(i) = compute_gaussian_intensity(Ax,kx,x_map(i),Phix)

             P_x(i) = P_x(i)*smooth(x_min,x_max,x_map(i))

             if(abs(P_x(i)).gt.P_x_max) then
                P_x_max = abs(P_x(i))
             end if

          end do


          ! deallocate unneeded arrays
          deallocate(kx)
          deallocate(Ax)
          deallocate(Phix)
          
        end subroutine compute_gaussian_perturbation_smooth


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the wave number for the sinusoidal functions
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param kx
        !> wave numbers for the different phase angles
        !---------------------------------------------------------------
        subroutine compute_wave_numbers(x_min,x_max,kx)

          implicit none

          real(rkind)              , intent(in)  :: x_min
          real(rkind)              , intent(in)  :: x_max
          real(rkind), dimension(:), intent(out) :: kx

          real(rkind) :: dk
          integer     :: i

          
          ! compute the wave number step
          dk = 2.0d0*pi/(x_max-x_min)


          ! compute the wave numbers
          do i=1, size(kx,1)

             kx(i) = i*dk

          end do

        end subroutine compute_wave_numbers


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the amplitudes of the sinusoidal fucntions
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param kx
        !> wave numbers for the different phase angles
        !
        !>@param Ax
        !> amplitudes of the sinusoidal functions
        !---------------------------------------------------------------
        subroutine compute_amplitudes(kx,Ax)

          implicit none

          real(rkind), dimension(:), intent(in)  :: kx
          real(rkind), dimension(:), intent(out) :: Ax

          
          integer     :: i
          real(rkind) :: dk
          real(rkind) :: kx_max_spectrum


          dk              = kx(2)-kx(1)
          kx_max_spectrum = (size(Ax,1)*dk)/4


          do i=1, size(Ax,1)

             Ax(i) = SQRT(2.0d0*power_spectrum(kx(i),kx_max_spectrum)*dk)

          end do

        end subroutine compute_amplitudes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the value of the power spectrum at a wave number
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param kx
        !> wave numbers for the different phase angles
        !
        !>@param kx_max_spectrum
        !> wave numbers where the power spectrum reaches its maximum
        !
        !>@return power_spectrum
        !> amplitude of the power spectrum
        !---------------------------------------------------------------
        function power_spectrum(kx,kx_max_spectrum)

          implicit none

          real(rkind), intent(in) :: kx
          real(rkind), intent(in) :: kx_max_spectrum
          real(rkind)             :: power_spectrum


          power_spectrum = (2.0d0*kx_max_spectrum)/(kx**2)*
     $                     EXP(-(2.0d0*kx_max_spectrum)/kx)

        end function power_spectrum


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the random phase angles
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param seed_x
        !> random seed used to initialize the random generator
        !
        !>@param Phix
        !> random phase angles
        !---------------------------------------------------------------
        subroutine compute_random_phases(seed_x,Phix)

          implicit none

          integer                  , intent(in)  :: seed_x
          real(rkind), dimension(:), intent(out) :: Phix


          integer, dimension(:), allocatable :: seed

          integer     :: seed_size
          integer     :: i
          real(rkind) :: x


          ! generate the correct seed to have random numbers
          ! over the size of Phix
          call random_seed(size=seed_size)
          allocate(seed(seed_size))
          call seedgen(0,seed_x,seed)


          !initialize the random generator with the seed
          call random_seed(put=seed)


          ! create the random phases
          do i=1, size(Phix,1)

             call RANDOM_NUMBER(x)
             Phix(i) = 2.0d0*pi*x
             
          end do

          deallocate(seed)

        end subroutine compute_random_phases


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> generate a random seed
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param pid
        !> pid of the processor
        !
        !>@param s
        !> seed used to initialize the random generator
        !
        !>@param seed
        !> random seed created for the random generator
        !---------------------------------------------------------------
        subroutine seedgen(pid,s,seed)
          
          implicit none

          integer              , intent(in)    :: pid
          integer              , intent(in)    :: s
          integer, dimension(:), intent(inout) :: seed

          integer :: i
          !call system_clock(s)
          !print *, s
          !s= 1384973980 !seconds since Epoch

          do i=1, size(seed)
             seed(i) = abs( mod((s*181)*((pid-83)*359), 104729) )
             seed(i) = seed(i) + 37*(i-1)
          end do

        end subroutine seedgen


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> sum the sinusoidal contributions to create the gaussian
        !> contribution at a coordinate
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param Ax
        !> amplitudes of the sinusoidal functions
        !
        !>@param kx
        !> wave numbers of the sinusoidal functions
        !
        !>@param x
        !> x-coordinate
        !
        !>@param Phix
        !> phase angles of the sinusoidal functions
        !---------------------------------------------------------------
        function compute_gaussian_intensity(Ax,kx,x,Phix)
     $     result(Px)

          implicit none
          
          real(rkind), dimension(:), intent(in) :: Ax
          real(rkind), dimension(:), intent(in) :: kx
          real(rkind)              , intent(in) :: x
          real(rkind), dimension(:), intent(in) :: Phix
          real(rkind)                           :: Px

          
          integer :: i

          Px = 0

          do i=1, size(Ax,1)

             Px = Px + Ax(i)*COS(kx(i)*x+Phix(i))

          end do

        end function compute_gaussian_intensity


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> smoothing factor at a coordinate
        !
        !> @date
        !> 09_04_2015 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param x_map
        !> x-coordinate map
        !
        !>@return smooth
        !> smoothing factor
        !---------------------------------------------------------------
        function smooth(x_min,x_max,x)

          implicit none

          real(rkind) , intent(in) :: x_min
          real(rkind) , intent(in) :: x_max
          real(rkind) , intent(in) :: x
          real(rkind)              :: smooth

          smooth = 1.0d0 - (2.0d0*(x-x_min)/(x_max-x_min)-1.0d0)**2

        end function smooth

      end module gaussian_perturbation_module
