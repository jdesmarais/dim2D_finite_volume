      !> @file
      !> class implementing the computation of the contribution
      !> of the hyperbolic terms to the time derivative of the
      !> conservative variables for the corner boundaries
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class implementing the computation of the contribution
      !> of the hyperbolic terms to the time derivative of the
      !> conservative variables for the corner boundaries
      !
      !> @date
      ! 08_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_corner_ns2d_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use lodi_corner_abstract_class, only :
     $       lodi_corner_abstract

        use ns2d_prim_module, only :
     $       compute_timedev_from_LODI_vectors

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_corner_ns2d


        !>@class lodi_corner_ns2d
        !> class encapsulating the interfaces for computing the
        !> contributions of the hyperbolic terms in the x- and
        !> y- directions at the corners of the computational domain
        !
        !>@param ini
        !> initialization of the functions describing the inlet
        !> or outlet flow (ex: u_in, v_in, P_out...)
        !
        !>@param compute_x_and_y_timedev
        !> compute the contributions to the time derivative of
        !> the LODI amplitudes in the x- and y- directions
        !---------------------------------------------------------------
        type, abstract, extends(lodi_corner_abstract) :: lodi_corner_ns2d

          contains

          procedure, pass :: compute_x_and_y_timedev

        end type lodi_corner_ns2d


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> conmpute the contribution of the hyperbolic terms in
        !> the x- and y- directions to the time derivatives
        !
        !> @date
        !> 08_09_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object with the edge conditions
        !
        !>@param p_model
        !> physical model needed to compute the eigenvalues
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param x_map
        !> map for the coordinates along the x-direction
        !
        !>@param y_map
        !> map for the coordinates along the y-direction
        !
        !>@param i
        !> index identifying the grid point along the x-direction
        !
        !>@param j
        !> index identifying the grid point along the y-direction
        !
        !>@param side_x
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain along the
        !> x-direction
        !
        !>@param side_y
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain along the
        !> y-direction
        !
        !>@param gradient_x
        !> procedure for the gradient computation along the x-direction
        !
        !>@param gradient_y
        !> procedure for the gradient computation along the y-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic tems along the x- and y-
        !> directions to the time derivatives at the grid point (i,j)
        !---------------------------------------------------------------
        function compute_x_and_y_timedev(
     $    this, p_model,
     $    t, nodes, x_map, y_map, i,j,
     $    side_x, side_y,
     $    gradient_x, gradient_y)
     $    result(timedev)

          implicit none
          
          class(lodi_corner_ns2d)      , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side_x
          logical                      , intent(in) :: side_y
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev


          real(rkind), dimension(ne) :: lodi_x
          real(rkind), dimension(ne) :: lodi_y


          !compute the lodi vectors in the x- and y-
          !directions
          call this%compute_x_and_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side_x, side_y,
     $         gradient_x, gradient_y,
     $         lodi_x, lodi_y)


          !compute the contributions of the hyperbolic terms in
          !the x- and y- directions to the time derivatives
          timedev = compute_timedev_from_LODI_vectors(
     $         nodes(i,j,:), lodi_x, lodi_y)


        end function compute_x_and_y_timedev

      end module lodi_corner_ns2d_class
