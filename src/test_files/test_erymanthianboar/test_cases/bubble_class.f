      module bubble_class
      
        use parameters_kind, only : ikind, rkind

        implicit none

        private
        public :: bubble


        type, abstract :: bubble

          real(rkind)               :: l_interface
          real(rkind)               :: radius
          real(rkind)               :: d_vap
          real(rkind)               :: d_liq
          real(rkind), dimension(2) :: center

          contains

          procedure               , pass           :: ini
          procedure               , pass           :: set_center
          procedure               , pass           :: set_radius
          procedure(update_proc)  , pass, deferred :: update
          procedure               , pass           :: get_mass_profile
          procedure               , pass           :: get_mass

        end type bubble


        abstract interface
          subroutine update_proc(this, velocity, dx, dy)

            import bubble
            import rkind

            class(bubble)            , intent(inout) :: this
            real(rkind), dimension(2), intent(in)    :: velocity
            real(rkind)              , intent(in)    :: dx
            real(rkind)              , intent(in)    :: dy

          end subroutine update_proc
        end interface        


        contains


        subroutine ini(this, center)
        
          implicit none

          class(bubble)                        , intent(inout) :: this
          real(rkind)  , dimension(2), optional, intent(in)    :: center


          this%l_interface = 0.1
          this%radius      = 0.3
          this%d_liq       = 1.1
          this%d_vap       = 0.1
          
          if(present(center)) then
             this%center = center
          else
             this%center = [0.0d0, 0.0d0]
          end if
          
        end subroutine ini


        subroutine set_center(this, center)

          implicit none

          class(bubble)            , intent(inout) :: this
          real(rkind), dimension(2), intent(in)    :: center

          this%center = center

        end subroutine set_center


        subroutine set_radius(this, radius)

          implicit none

          class(bubble), intent(inout) :: this
          real(rkind)  , intent(in)    :: radius

          this%radius = radius

        end subroutine set_radius


        function get_mass_profile(this, coords) result(mass)

          implicit none

          class(bubble)            , intent(in) :: this
          real(rkind), dimension(2), intent(in) :: coords
          real(rkind)                           :: mass

          real(rkind) :: r

          r = SQRT((coords(1)-this%center(1))**2+(coords(2)-this%center(2))**2)

          mass = (this%d_liq+this%d_vap)/2.0d0 +
     $           (this%d_liq-this%d_vap)/2.0d0*
     $               Tanh(-2*(r-this%radius)/this%l_interface)

        end function get_mass_profile


        function get_mass(this, coords) result(mass)
          
          implicit none

          class(bubble)            , intent(in) :: this
          real(rkind), dimension(2), intent(in) :: coords
          real(rkind)                           :: mass

          mass = this%get_mass_profile(coords)

        end function get_mass

      end module bubble_class
