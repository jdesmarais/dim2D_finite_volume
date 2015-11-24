      !> @file
      !> definition of 'ikind' and 'rkind' determining
      !> the type of the integers and reals at
      !> compilation time
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> definition of 'ikind' and 'rkind' determining
      !> the type of the integers and reals at
      !> compilation time
      !
      !> @date
      !> 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_kind

        implicit none

        private
        public :: rkind, ikind

        integer, parameter :: digits=8   !<@brief number of digits over which a real is saved
        integer, parameter :: decades=8  !<@brief number of digits over which an integer is saved

        integer, parameter :: rkind = selected_real_kind(digits) !<@brief Fortran real type for the program
        integer, parameter :: ikind = selected_int_kind(decades) !<@brief Fortran integer type for the program

      end module parameters_kind
