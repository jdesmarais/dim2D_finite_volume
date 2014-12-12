      !> @file
      !> class implementing the computation of the contribution
      !> of the hyperbolic terms to the time derivative of the
      !> conservative variables
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class implementing the computation of the contribution
      !> of the hyperbolic terms to the time derivative of the
      !> conservative variables
      !
      !> @date
      ! 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_edge_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use lodi_edge_abstract_class, only :
     $       lodi_edge_abstract

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_edge


        !>@class lodi_edge_abstract
        !> class encapsulating the interfaces for computing the
        !> LODI amplitudes in the x and y directions at the edge
        !> of the computational domain (boundary layers except the
        !> corners)
        !
        !>@param ini
        !> initialization of the functions describing the inlet
        !> or outlet flow (ex: u_in, v_in, P_out...)
        !
        !>@param compute_x_lodi
        !> compute the LODI amplitudes in the x-direction
        !
        !>@param compute_y_lodi
        !> compute the LODI amplitudes in the y-direction
        !
        !>@param compute_x_timedev
        !> compute the contribution to the time derivative of
        !> the LODI amplitudes in the x-direction
        !
        !>@param compute_y_timedev
        !> compute the contribution to the time derivative of
        !> the LODI amplitudes in the y-direction
        !---------------------------------------------------------------
        type, abstract, extends(lodi_edge_abstract) :: lodi_edge

          contains

          procedure, pass :: compute_x_timedev
          procedure, pass :: compute_y_timedev

        end type lodi_edge


        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> conmpute the contribution of the hyperbolic terms in
        !> the x-direction to the time derivatives
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
        !>@param transverse_lodi
        !> transverse LODI vector at (i,j)
        !
        !>@param viscous_lodi
        !> viscous LODI vector at (i,j)
        !
        !>@param side
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic tems along the x-direction to
        !> the time derivatives at the grid point (i,j)
        !---------------------------------------------------------------
        function compute_x_timedev(
     $       this, p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(timedev)

          implicit none
          
          class(lodi_edge)             , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(ne)   , intent(in) :: transverse_lodi
          real(rkind), dimension(ne)   , intent(in) :: viscous_lodi
          logical                      , intent(in) :: side
          procedure(gradient_x_proc)                :: gradient
          real(rkind), dimension(ne)                :: timedev


          real(rkind), dimension(ne) :: lodi          


          !compute the lodi vector
          lodi = this%compute_x_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         transverse_lodi, viscous_lodi,
     $         side,
     $         gradient)


          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the conservative
          !variables
          timedev = p_model%compute_x_timedev_from_LODI_vector(
     $         nodes(i,j,:),
     $         lodi)

        end function compute_x_timedev


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> conmpute the contribution of the hyperbolic terms in
        !> the y-direction to the time derivatives
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
        !>@param transverse_lodi
        !> transverse LODI vector at (i,j)
        !
        !>@param viscous_lodi
        !> viscous LODI vector at (i,j)
        !
        !>@param side
        !> boolean designating whether the procedure is applied at a
        !> left or right boundary of the computational domain
        !
        !>@param gradient
        !> procedure for the gradient computation along the x-direction
        !
        !>@return timedev
        !> contribution of the hyperbolic tems along the y-direction to
        !> the time derivatives at the grid point (i,j)
        !---------------------------------------------------------------
        function compute_y_timedev(
     $       this, p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(timedev)

          implicit none
          
          class(lodi_edge)             , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(ne)   , intent(in) :: transverse_lodi
          real(rkind), dimension(ne)   , intent(in) :: viscous_lodi
          logical                      , intent(in) :: side
          procedure(gradient_y_proc)                :: gradient
          real(rkind), dimension(ne)                :: timedev


          real(rkind), dimension(ne) :: lodi          


          !compute the lodi vector
          lodi = this%compute_y_lodi(
     $         p_model,
     $         t, nodes, x_map, y_map, i,j,
     $         transverse_lodi, viscous_lodi,
     $         side,
     $         gradient)


          !compute the contribution of the hyperbolic terms along the
          !x-direction to the time derivatives of the conservative
          !variables
          timedev = p_model%compute_y_timedev_from_LODI_vector(
     $         nodes(i,j,:),
     $         lodi)

        end function compute_y_timedev

      end module lodi_edge_class
