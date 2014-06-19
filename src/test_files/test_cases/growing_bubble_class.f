      module growing_bubble_class

        use bubble_class   , only : bubble
        use parameters_kind, only : rkind

        implicit none

        private
        public :: growing_bubble


        type, extends(bubble) :: growing_bubble

          contains

          procedure, pass :: update

        end type growing_bubble


        contains


        subroutine update(this, velocity, dx, dy)

          implicit none

          class(growing_bubble)  , intent(inout) :: this
          real(rkind), dimension(2), intent(in)    :: velocity
          real(rkind)              , intent(in)    :: dx
          real(rkind)              , intent(in)    :: dy
          
          this%radius = this%radius + 
     $         SQRT(velocity(1)**2+velocity(2)**2)*
     $         SQRT(0.5d0*(dx**2+dy**2))
          
        end subroutine update        

      end module growing_bubble_class
