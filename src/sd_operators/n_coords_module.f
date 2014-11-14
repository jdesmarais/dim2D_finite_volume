      module n_coords_module

        use parameters_kind, only :
     $     rkind

        implicit none

        private
        public :: 
     $       get_x_coord,
     $       get_y_coord,
     $       get_n1_coord,
     $       get_n2_coord
        
        contains


        function get_x_coord(n1,n2) result(x)

          implicit none
          
          real(rkind), intent(in) :: n1
          real(rkind), intent(in) :: n2
          real(rkind)             :: x

          x = SQRT(2.0d0)*0.25d0*(n1+n2)

        end function get_x_coord

        
        function get_y_coord(n1,n2) result(y)

          implicit none

          real(rkind), intent(in) :: n1
          real(rkind), intent(in) :: n2
          real(rkind)             :: y

          y = SQRT(2.0d0)*0.25d0*(n2-n1)

        end function get_y_coord


        function get_n1_coord(x,y) result(n1)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: n1

          n1 = SQRT(2.0d0)*0.5d0*(x-y)

        end function get_n1_coord


        function get_n2_coord(x,y) result(n2)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: n2

          n2 = SQRT(2.0d0)*0.5d0*(x+y)

        end function get_n2_coord

      end module n_coords_module
