      !> @file
      !> module implementing the abstract version of
      !> the buffer mainlayer object to encapsulate
      !> a double chained list of buffer sublayer
      !> elements
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the abstract version of
      !> the buffer mainlayer object to encapsulate
      !> a double chained list of buffer sublayer
      !> elements
      !
      !> @date
      ! 09_05_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_abstract_class

        use bf_sublayer_class, only : bf_sublayer

        implicit none

        private
        public  :: bf_mainlayer_abstract
        

        !> @class bf_mainlayer_abstract
        !> class encapsulating the buffer sublayers corresponding
        !> to the same cardinal point (N,S,E,W,NE,NW,SE,SW)
        !>
        !> @param mainlayer_id
        !> cardinal coordinate identifying the position of
        !> the mainlayer compared to the interior domain
        !>
        !> @param nb_sublayers
        !> number of sublayers saved in the main layer
        !>
        !> @param head_sublayer
        !> pointer of the head sublayer of the main layer
        !>
        !> @param tail_sublayer
        !> pointer of the tail sublayer of the main layer
        !---------------------------------------------------------------
        type :: bf_mainlayer_abstract

          integer :: mainlayer_id
          integer :: nb_sublayers

          type(bf_sublayer), pointer :: head_sublayer
          type(bf_sublayer), pointer :: tail_sublayer

        end type bf_mainlayer_abstract

      end module bf_mainlayer_abstract_class
