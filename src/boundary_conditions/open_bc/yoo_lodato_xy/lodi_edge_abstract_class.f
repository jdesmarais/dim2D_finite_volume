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
      ! 04_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module lodi_edge_abstract_class

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
        public :: lodi_edge_abstract


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
        type, abstract :: lodi_edge_abstract

          contains

          procedure(ini_proc)     ,   pass, deferred :: ini
          procedure(xlodi_proc)   , nopass, deferred :: compute_x_lodi
          procedure(ylodi_proc)   , nopass, deferred :: compute_y_lodi

          procedure(xtimedev_proc),   pass, deferred :: compute_x_timedev
          procedure(ytimedev_proc),   pass, deferred :: compute_y_timedev

        end type lodi_edge_abstract


        abstract interface

           !initialization of the lodi b.c.
           !(ex: inlet functions u_in,v_in,T_in)
           subroutine ini_proc(this)

             import lodi_edge_abstract

             class(lodi_edge_abstract), intent(inout) :: this

           end subroutine ini_proc


           !computation of the LODI amplitudes in the x-direction
           function xlodi_proc(
     $       p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(lodi)

             import gradient_x_proc
             import lodi_edge_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
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
             real(rkind), dimension(ne)                :: lodi

           end function xlodi_proc


           !computation of the LODI amplitudes in the y-direction
           function ylodi_proc(
     $       p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(lodi)

             import gradient_y_proc
             import lodi_edge_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
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
             real(rkind), dimension(ne)                :: lodi

           end function ylodi_proc


           !computation of the time derivatives for a boundary normal
           !to the x-direction
           function xtimedev_proc(
     $       this, p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(timedev)

             import gradient_x_proc
             import lodi_edge_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
             class(lodi_edge_abstract)    , intent(in) :: this
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

           end function xtimedev_proc


           !computation of the time derivatives for a boundary normal
           !to the y-direction
           function ytimedev_proc(
     $       this, p_model,
     $       t, nodes, x_map, y_map, i,j,
     $       transverse_lodi, viscous_lodi,
     $       side,
     $       gradient)
     $       result(timedev)

             import gradient_y_proc
             import lodi_edge_abstract
             import ikind
             import ne
             import pmodel_eq
             import rkind
             
             class(lodi_edge_abstract)    , intent(in) :: this
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

           end function ytimedev_proc

        end interface

      end module lodi_edge_abstract_class
