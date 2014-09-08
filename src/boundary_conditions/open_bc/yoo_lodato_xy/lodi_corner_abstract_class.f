      !> @file
      !> class encapsulating the interfaces for the computation
      !> of the LODI amplitudes in the x and y directions according
      !> to the procedure designed by Yoo et al. in "Characteristic
      !> boundary conditions for direct simulations of turbulent
      !> counterflow flames", Combustion Theory and Modelling,
      !> Vol 9, No. 4, pp 617 - 646, 2012
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the interfaces for the computation
      !> of the LODI amplitudes in the x and y directions
      !
      !> @date
      ! 08_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_corner_abstract_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        private
        public :: lodi_corner_abstract


        !>@class lodi_corner_abstract
        !> class encapsulating the interfaces for computing the
        !> LODI amplitudes in the x and y directions at the corner
        !> of the computational domain
        !
        !>@param ini
        !> initialization of the functions describing the inlet
        !> or outlet flow (ex: u_in, v_in, P_out...)
        !
        !>@param compute_x_and_y_lodi
        !> compute the LODI amplitudes in the x-and y-directions
        !
        !>@param compute_x_and_y_timedev
        !> compute the contributions to the time derivative of
        !> the LODI amplitudes in the x- and y- directions
        !---------------------------------------------------------------
        type, abstract :: lodi_corner_abstract

          contains

          procedure(ini_proc)      ,   pass, deferred :: ini
          procedure(xylodi_proc)   , nopass, deferred :: compute_x_and_y_lodi
          procedure(xytimedev_proc),   pass, deferred :: compute_x_and_y_timedev

        end type lodi_corner_abstract


        abstract interface

           !initialization of the lodi b.c.
           !(ex: inlet functions u_in,v_in,T_in)
           subroutine ini_proc(this)

             import lodi_corner_abstract

             class(lodi_corner_abstract), intent(inout) :: this

           end subroutine ini_proc


           !computation of the LODI amplitudes in the x-direction
           subroutine xylodi_proc(
     $       p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       side_x, side_y,
     $       gradient_x, gradient_y,
     $       lodi_x, lodi_y)

             import gradient_x_proc
             import gradient_y_proc
             import lodi_corner_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
             type(pmodel_eq)              , intent(in)  :: p_model
             real(rkind)                  , intent(in)  :: t
             real(rkind), dimension(:,:,:), intent(in)  :: nodes
             real(rkind), dimension(:)    , intent(in)  :: x_map
             real(rkind), dimension(:)    , intent(in)  :: y_map
             integer(ikind)               , intent(in)  :: i
             integer(ikind)               , intent(in)  :: j
             logical                      , intent(in)  :: side_x
             logical                      , intent(in)  :: side_y
             procedure(gradient_x_proc)                 :: gradient_x
             procedure(gradient_y_proc)                 :: gradient_y
             real(rkind), dimension(ne)   , intent(out) :: lodi_x
             real(rkind), dimension(ne)   , intent(out) :: lodi_y

           end subroutine xylodi_proc


           !computation of the time derivatives for a boundary normal
           !to the x-direction
           function xytimedev_proc(
     $       this, p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       side_x, side_y,
     $       gradient_x, gradient_y)
     $       result(timedev)

             import gradient_x_proc
             import gradient_y_proc
             import lodi_corner_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
             class(lodi_corner_abstract)  , intent(in) :: this
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

           end function xytimedev_proc

        end interface

      end module lodi_corner_abstract_class
