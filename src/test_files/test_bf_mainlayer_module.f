      module test_bf_mainlayer_module

        use bf_mainlayer_class, only : bf_mainlayer
        use bf_sublayer_class , only : bf_sublayer

        implicit none

        private
        public :: print_mainlayer

        contains

        subroutine print_mainlayer(bf_mainlayer_tested)

          implicit none

          class(bf_mainlayer), intent(inout) :: bf_mainlayer_tested

          type(bf_sublayer), pointer :: current_sublayer
          integer :: i
          
          print '(''mainlayer ID: '',I2)', bf_mainlayer_tested%get_mainlayer_id()
          print '(''nb_sublayers: '',I2)', bf_mainlayer_tested%get_nb_sublayers()
          print '('''')'

          i=1
          current_sublayer => bf_mainlayer_tested%get_head_sublayer()
          print '(''sublayer: '', I2, '', location: '', I3,
     $         '', alignment: '', 4I3)',
     $         i,
     $         current_sublayer%get_localization(),
     $         current_sublayer%get_alignment(1,1), current_sublayer%get_alignment(1,2),
     $         current_sublayer%get_alignment(2,1), current_sublayer%get_alignment(2,2)

          do while(associated(current_sublayer%get_next()))

             i=i+1
             current_sublayer => current_sublayer%get_next()
             print '(''sublayer: '', I2, '', location: '', I3,
     $         '', alignment: '', 4I3)',
     $            i,
     $            current_sublayer%get_localization(),
     $            current_sublayer%get_alignment(1,1), current_sublayer%get_alignment(1,2),
     $            current_sublayer%get_alignment(2,1), current_sublayer%get_alignment(2,2)

          end do
          print '('''')'

        end subroutine print_mainlayer

      end module test_bf_mainlayer_module
