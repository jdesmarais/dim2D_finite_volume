      !> @file
      !> module encapsulating the buffer layer object where only
      !> its main attributes are implemented
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the buffer layer object where only
      !> its main attributes are implemented
      !
      !> @date
      ! 07_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_abstract_class

        use parameters_kind    , only : ikind, rkind

        implicit none


        private
        public :: bf_layer_abstract


        logical, parameter :: debug = .true.


        !> @class bf_layer_abstract
        !> class encapsulating the buffer layer object but
        !> only its main attributes are implemented
        !>
        !> @param localization
        !> localization of the buffer layer
        !> the first element gives the main layer (N,S,E,W,
        !> N_E,N_W,S_E,S_W), the second element gives its
        !> position inside the main layer
        !>
        !> @param alignment
        !> gives the buffer layer position compared to the
        !> interior grid points:
        !>     - alignment(1,1) : i_min
        !>     - alignment(1,2) : i_max
        !>     - alignment(2,1) : j_min
        !>     - alignment(2,2) : j_max
        !>
        !> @param nodes
        !> table where the data needed in the governing
        !> equations are saved
        !>
        !> @param grdpts_id
        !> table where the grid points ID are saved
        !> (the ID tells whether the grid point is not
        !> computed (no_pt) or it is an interior grid point
        !> computed using the normal stencil (interior_pt),
        !> or it is a point at the interface between the interior
        !> and the boundary points (bc_interior_pt) or it is a
        !> boundary point (bc_pt) or a grid point exchanged with
        !> the neighboring buffer layers (exchange_pt)
        !> 
        !> @param print_sizes
        !> procedure used to print the sizes of the main tables of the
        !> buffer layer
        !---------------------------------------------------------------
        type :: bf_layer_abstract

          integer                       , private :: localization
          integer(ikind), dimension(2,2), private :: alignment
          
          real(rkind)   , dimension(:,:,:), allocatable, private :: nodes
          integer       , dimension(:,:)  , allocatable, private :: grdpts_id

          contains

          procedure, pass :: ini

        end type bf_layer_abstract


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine print the nodes table in a binary
        !> file
        !
        !> @date
        !> 07_04_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_layer_abstract object encapsulating the main
        !> tables and the integer identifying the
        !> correspondance between the buffer layer and the
        !> interior grid points
        !
        !>@param localization
        !> localization of the buffer layer
        !>     - localization(1) : main layer (N,S,E,W,N_E,N_W,S_E,S_W)
        !>     - localization(2) : sub-layer
        !--------------------------------------------------------------        
        subroutine ini(this,localization)

          implicit none

          class(bf_layer_abstract)    , intent(inout) :: this
          integer(ikind)              , intent(in)    :: localization
          
          this%localization = localization

        end subroutine ini


        !!> @author
        !!> Julien L. Desmarais
        !!
        !!> @brief
        !!> subroutine print the nodes table in a binary
        !!> file
        !!
        !!> @date
        !!> 07_04_2013 - initial version - J.L. Desmarais
        !!
        !!>@param this
        !!> bf_layer_abstract object encapsulating the main
        !!> tables and the integer identifying the
        !!> correspondance between the buffer layer and the
        !!> interior grid points
        !!
        !!>@param localization
        !!> return the localization of the buffer layer
        !!--------------------------------------------------------------        
        !function get_localization(this)
        !
        !  implicit none
        !
        !  class(bf_layer_abstract)    , intent(in) :: this
        !  integer(ikind)                           :: localization
        !  
        !  localization = this%localization
        !
        !end function get_localization

      end module bf_layer_abstract_class
