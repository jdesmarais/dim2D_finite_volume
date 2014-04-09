      !> @file
      !> module implementing the object encapsulating data
      !> needed for the exchange of information between the
      !> interior domain and the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the abstract version of the
      !> interface between the interior domain and the buffer
      !> layers: only the attributes are implemented
      !
      !> @date
      ! 09_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module interface_abstract_class

        use parameters_kind, only : ikind, rkind
        use bf_layer_class , only : bf_layer    

        implicit none


        private
        public :: interface_abstract


        logical, parameter :: debug = .true.


        !> @class interface_abstract
        !> class encapsulating the bf_layer/interior interface
        !> object but only its main attributes are implemented
        !>
        !> @param N_layers
        !> allocatable table which contains the different
        !> sublayers for the north buffer layer
        !>
        !> @param S_layers
        !> allocatable table which contains the different
        !> sublayers for the south buffer layer
        !>
        !> @param E_layers
        !> allocatable table which contains the different
        !> sublayers for the east buffer layer
        !>
        !> @param W_layers
        !> allocatable table which contains the different
        !> sublayers for the west buffer layer
        !>
        !> @param NW_layer
        !> buffer layer representing the NW corner
        !>
        !> @param NE_layer
        !> buffer layer representing the NE corner
        !>
        !> @param SW_layer
        !> buffer layer representing the SW corner
        !>
        !> @param SE_layer
        !> buffer layer representing the SE corner
        !---------------------------------------------------------------
        type :: interface_abstract

          type(bf_layer), dimension(:), allocatable :: N_bf_layers
          type(bf_layer), dimension(:), allocatable :: S_bf_layers
          type(bf_layer), dimension(:), allocatable :: E_bf_layers
          type(bf_layer), dimension(:), allocatable :: W_bf_layers
          
          type(bf_layer) :: NW_bf_layer
          type(bf_layer) :: NE_bf_layer
          type(bf_layer) :: SW_bf_layer
          type(bf_layer) :: SE_bf_layer

        end type interface_abstract
      
