      !> @file
      !> class encapsulating the interfaces for the computation
      !> of the time derivatives from the LODI amplitudes in the
      !> x and y directions for the NS equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the interfaces for the computation
      !> of the time derivatives from the LODI amplitudes in the
      !> x and y directions for the NS equations
      !
      !> @date
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_ns2d_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use lodi_abstract_class, only :
     $       lodi_abstract

        use ns2d_prim_module, only :
     $       compute_x_timedev_from_LODI_vector,
     $       compute_y_timedev_from_LODI_vector

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_ns2d


        !>@class lodi_ns2d
        !> class encapsulating the interfaces for the computation
        !> of the time derivatives from the LODI amplitudes in the
        !> x and y directions for the NS equations
        !---------------------------------------------------------------
        type, abstract, extends(lodi_abstract) :: lodi_ns2d

          contains

          procedure, pass :: compute_x_timedev
          procedure, pass :: compute_y_timedev

        end type lodi_ns2d


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms of the
        !> governing equations in the x-direction to the time
        !> derivatives of the governing variables
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model governing equations
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param gradient
        !> procedure computing the gradient along the x-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic terms in the x-direction to
        !> the time derivatives
        !-------------------------------------------------------------
        function compute_x_timedev(
     $    this, p_model,
     $    t, nodes, x_map, y_map, i,j,
     $    side,
     $    gradient)
     $    result(timedev)

          implicit none
          
          class(lodi_ns2d)             , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_x_proc)                :: gradient
          real(rkind), dimension(ne)                :: timedev

          
          real(rkind), dimension(ne) :: lodi          


          !compute the lodi vector
          lodi = this%compute_x_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side,
     $         gradient)

          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the conservative variables
          timedev = compute_x_timedev_from_LODI_vector(nodes(i,j,:),lodi)

        end function compute_x_timedev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the contribution of the hyperbolic terms of the
        !> governing equations in the y-direction to the time
        !> derivatives of the governing variables
        !
        !> @date
        !> 13_08_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model governing equations
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid points
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param i
        !> index identifying the grid point position along the x-axis
        !
        !>@param j
        !> index identifying the grid point position along the y-axis
        !
        !>@param gradient
        !> procedure computing the gradient along the x-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic terms in the x-direction to
        !> the time derivatives
        !-------------------------------------------------------------
        function compute_y_timedev(
     $     this, p_model,
     $     t, nodes, x_map, y_map, i,j,
     $     side,
     $     gradient)
     $     result(timedev)

          implicit none
          
          class(lodi_ns2d)             , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side
          procedure(gradient_y_proc)                :: gradient
          real(rkind), dimension(ne)                :: timedev

          real(rkind), dimension(ne)    :: lodi

          !compute the lodi vector
          lodi = this%compute_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         side,
     $         gradient)
          
          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the conservative variables
          timedev = compute_y_timedev_from_LODI_vector(nodes(i,j,:),lodi)

        end function compute_y_timedev

      end module lodi_ns2d_class
