      module bf_mainlayer_class

        use bf_layer_class, only : bf_layer

        implicit none

        private
        public :: bf_mainlayer
        
        
        type :: bf_mainlayer

          type(bf_layer), dimension(:), allocatable :: mainlayer

        end type :: bf_mainlayer


      end module bf_mainlayer_class
