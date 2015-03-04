      !> @file
      !> annex subroutines to decide whether boundary sections should
      !> be computed in the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to decide whether boundary
      !> layers should be computed
      !
      !> @date
      ! 28_01_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
       module bf_layer_bc_checks_module

        use parameters_input, only :
     $       nx,ny,ne,
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

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the N edge should be computed depending
        !> on the type of boundary conditions used on the north
        !> layer
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param j
        !> y-index identifying the grid point position in the y_map
        !
        !> @param y_map
        !> map of the y-coordinates along the y-axis
        !
        !> @param bc_type_choice
        !> integer identifying the type of boundary condition
        !
        !> @return compute_edge
        !> logical indicating whether the boundary layer should be
        !> computed
        !--------------------------------------------------------------
        function compute_edge_N(y,bc_type_choice) result(compute_edge)

          implicit none

          real(rkind), intent(in) :: y
          integer    , intent(in) :: bc_type_choice
          logical                 :: compute_edge
          
          compute_edge = .not.(
     $         (y.ge.y_max).and.
     $         (bc_N_type_choice.ne.bc_type_choice))

        end function compute_edge_N


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the S edge should be computed depending
        !> on the type of boundary conditions used on the south
        !> layer
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param j
        !> y-index identifying the grid point position in the y_map
        !
        !> @param y_map
        !> map of the y-coordinates along the y-axis
        !
        !> @param bc_type_choice
        !> integer identifying the type of boundary condition
        !
        !> @return compute_edge
        !> logical indicating whether the boundary layer should be
        !> computed
        !--------------------------------------------------------------
        function compute_edge_S(y,bc_type_choice) result(compute_edge)

          implicit none

          real(rkind), intent(in) :: y
          integer    , intent(in) :: bc_type_choice
          logical                 :: compute_edge
          
          compute_edge = .not.(
     $         (y.le.y_min).and.
     $         (bc_S_type_choice.ne.bc_type_choice))

        end function compute_edge_S


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the E edge should be computed depending
        !> on the type of boundary conditions used on the east
        !> layer
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param i
        !> x-index identifying the grid point position in the x_map
        !
        !> @param x_map
        !> map of the x-coordinates along the x-axis
        !
        !> @param bc_type_choice
        !> integer identifying the type of boundary condition
        !
        !> @return compute_edge
        !> logical indicating whether the boundary layer should be
        !> computed
        !--------------------------------------------------------------
        function compute_edge_E(x,bc_type_choice) result(compute_edge)

          implicit none

          real(rkind), intent(in) :: x
          integer    , intent(in) :: bc_type_choice
          logical                 :: compute_edge
          
          compute_edge = .not.(
     $         (x.ge.x_max).and.
     $         (bc_E_type_choice.ne.bc_type_choice))

        end function compute_edge_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the W edge should be computed depending
        !> on the type of boundary conditions used on the west
        !> layer
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param i
        !> x-index identifying the grid point position in the x_map
        !
        !> @param x_map
        !> map of the x-coordinates along the x-axis
        !
        !> @param bc_type_choice
        !> integer identifying the type of boundary condition
        !
        !> @return compute_edge
        !> logical indicating whether the boundary layer should be
        !> computed
        !--------------------------------------------------------------
        function compute_edge_W(x,bc_type_choice) result(compute_edge)

          implicit none

          real(rkind), intent(in) :: x
          integer    , intent(in) :: bc_type_choice
          logical                 :: compute_edge
          
          compute_edge = .not.(
     $         (x.le.x_min).and.
     $         (bc_W_type_choice.ne.bc_type_choice))

        end function compute_edge_W

      end module bf_layer_bc_checks_module
