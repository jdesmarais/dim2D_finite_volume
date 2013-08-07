      module parameters_kind

        implicit none

        private
        public :: rkind, ikind

        integer, parameter :: digits=8
        integer, parameter :: decades=9

        integer, parameter :: rkind = selected_real_kind(digits)
        integer, parameter :: ikind = selected_int_kind(decades)

      end module parameters_kind
