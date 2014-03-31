      !> @file
      !> definition of the type of the integers and reals used
      !> in the program
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> definition of 'ikind' and 'rkind' deciding the type of the
      !> integers and reals at compilation time
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_kind

        implicit none

        private
        public :: rkind, ikind

        integer, parameter :: digits=8
        integer, parameter :: decades=9

        integer, parameter :: rkind = selected_real_kind(digits)
        integer, parameter :: ikind = selected_int_kind(decades)

      end module parameters_kind
