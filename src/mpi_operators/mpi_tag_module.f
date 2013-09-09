      !> @file
      !> module encapsulating subroutines to compute the tag
      !> need to recognize the messages sent between the tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute the tag
      !> need to recognize the messages sent between the tiles
      !
      !> @date
      ! 23_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_tag_module

        implicit none

        private
        public :: compute_mpi_tag

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> create a unique tag for an mpi request to be able to track
        !> it later
        !
        !> @date
        !> 03_07_2012 - initial version - J.L. Desmarais
        !> 23_08_2013 - adaptation      - J.L. Desmarais
        !
        !>@param rank_send
        !> rank of the processor sending the information
        !
        !>@param rank recv
        !> rank of the processor receiving the information
        !
        !>@param nb_procs
        !> total number of processors
        !
        !>@param tag
        !> tag for the MPI request
        !--------------------------------------------------------------
        function compute_mpi_tag(
     $     rank_send,
     $     rank_recv,
     $     nb_procs)
     $     result(tag)


          implicit none

          integer, intent(in) :: rank_send
          integer, intent(in) :: rank_recv
          integer, intent(in) :: nb_procs
          integer             :: tag


          !procedure
          !the tag has to be unique
          !as rank_send \in [0, nb_procs-1]
          !as rank_recv \in [0, nb_procs-1]
          !rank*nb_procs + rank_recv \in [0, nb_procs^2-1]
          !then the euclidian division of tag by nb_procs^2 defines
          !uniquely comm_id and (rank_send*nb_procs + rank_recv)
          !then the euclidian division of (rank_send*nb_procs + rank_recv)
          !by nb_procs gives uniquely rank_send and rank_recv

          !tag = comm_id*nb_procs**2 + (rank_send*nb_procs + rank_recv)
          tag = rank_send*nb_procs + rank_recv

        end function compute_mpi_tag

      end module mpi_tag_module
