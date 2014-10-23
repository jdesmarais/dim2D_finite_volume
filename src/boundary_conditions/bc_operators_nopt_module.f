      module bc_operators_nopt_module

        use parameters_input, only :
     $       x_min, x_max,
     $       y_min, y_max,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public ::
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W

        contains

        !check whether the N edge should be computed depending
        !on the type of boundary conditions used on the north
        !layer
        function compute_edge_N(j,y_map,bc_type_choice) result(compute_edge)

          implicit none

          integer(ikind)           , intent(in) :: j
          real(rkind), dimension(:), intent(in) :: y_map
          integer                  , intent(in) :: bc_type_choice
          logical                               :: compute_edge
          
          compute_edge = .not.(
     $         (y_map(j).ge.y_max).and.
     $         (bc_N_type_choice.ne.bc_type_choice))

        end function compute_edge_N


        !check whether the S edge should be computed depending
        !on the type of boundary conditions used on the south
        !layer
        function compute_edge_S(j,y_map,bc_type_choice) result(compute_edge)

          implicit none

          integer(ikind)           , intent(in) :: j
          real(rkind), dimension(:), intent(in) :: y_map
          integer                  , intent(in) :: bc_type_choice
          logical                               :: compute_edge
          
          compute_edge = .not.(
     $         (y_map(j+1).le.y_min).and.
     $         (bc_S_type_choice.ne.bc_type_choice))

        end function compute_edge_S


        !check whether the E edge should be computed depending
        !on the type of boundary conditions used on the east
        !layer
        function compute_edge_E(i,x_map,bc_type_choice) result(compute_edge)

          implicit none

          integer(ikind)           , intent(in) :: i
          real(rkind), dimension(:), intent(in) :: x_map
          integer                  , intent(in) :: bc_type_choice
          logical                               :: compute_edge
          
          compute_edge = .not.(
     $         (x_map(i).ge.y_max).and.
     $         (bc_E_type_choice.ne.bc_type_choice))

        end function compute_edge_E


        !check whether the W edge should be computed depending
        !on the type of boundary conditions used on the west
        !layer
        function compute_edge_W(i,x_map,bc_type_choice) result(compute_edge)

          implicit none

          integer(ikind)           , intent(in) :: i
          real(rkind), dimension(:), intent(in) :: x_map
          integer                  , intent(in) :: bc_type_choice
          logical                               :: compute_edge
          
          compute_edge = .not.(
     $         (x_map(i+1).le.x_min).and.
     $         (bc_W_type_choice.ne.bc_type_choice))

        end function compute_edge_W

      end module bc_operators_nopt_module
