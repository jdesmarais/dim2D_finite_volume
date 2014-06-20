      module rising_bubble_class

        use bubble_class   , only : bubble
        use parameters_kind, only : ikind, rkind

        implicit none

        private
        public :: rising_bubble


        type, extends(bubble) :: rising_bubble

          integer :: nb_bubbles

          contains

          procedure, pass :: set_nb_bubbles
          procedure, pass :: get_nb_bubbles
          procedure, pass :: update
          procedure, pass :: get_mass

        end type rising_bubble


        contains


        subroutine set_nb_bubbles(this, nb_bubbles)

          implicit none
          
          class(rising_bubble), intent(inout) :: this
          integer             , intent(in)    :: nb_bubbles

          this%nb_bubbles = nb_bubbles

        end subroutine set_nb_bubbles


        function get_nb_bubbles(this) result(nb_bubbles)

          implicit none

          class(rising_bubble), intent(inout) :: this
          real(rkind)                         :: nb_bubbles

          nb_bubbles = this%nb_bubbles

        end function get_nb_bubbles


        subroutine update(this, velocity, dx, dy)

          implicit none

          class(rising_bubble)     , intent(inout) :: this
          real(rkind), dimension(2), intent(in)    :: velocity
          real(rkind)              , intent(in)    :: dx
          real(rkind)              , intent(in)    :: dy
          
          this%center(1) = this%center(1) + velocity(1)*dx
          this%center(2) = this%center(2) + velocity(2)*dy
          
        end subroutine update


        function get_mass(this, coords) result(mass)
          
          implicit none

          class(rising_bubble)     , intent(in) :: this
          real(rkind), dimension(2), intent(in) :: coords
          real(rkind)                           :: mass

          real(rkind), dimension(2) :: coords_temp

          select case(this%nb_bubbles)
            case(1)
               mass = this%get_mass_profile(coords)
            case(2)
               coords_temp(1) =  abs(coords(1))
               coords_temp(2) =  coords(2)
               mass = this%get_mass_profile(coords_temp)
            case default
               stop 'rising_bubble: get_mass: nb_bubble < 3'
          end select
          
        end function get_mass

      end module rising_bubble_class
