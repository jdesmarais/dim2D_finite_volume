      !> @file
      !> module containing the constants used in the
      !> program: they will be propagated at compilation
      !> time
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> definition of the constant used in the whole program
      !
      !> @date
      !> 08_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_constant

        use parameters_kind, only : rkind, ikind

        !<program version
        character*(*) , parameter :: prog_version = 'lerneanhydra V1.0'

        !<main variable types
        integer(ikind), parameter :: scalar=0
        integer(ikind), parameter :: vector_x=1
        integer(ikind), parameter :: vector_y=2

      end module parameters_constant
