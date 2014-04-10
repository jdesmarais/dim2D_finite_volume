      module bf_layer_pointer_class

        use bf_layer_class, only : bf_layer

        implicit none

        private
        public :: bf_layer_pointer


        type :: bf_layer_pointer

          type(bf_layer), pointer :: pointer_to_bf_layer

        end type bf_layer_pointer



      end module bf_layer_pointer_class
