      !> @file
      !> class encapsulating attribute to identify the procedures
      !> for the computation of the boundary layers
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating attributes to identify the procedures
      !> for the computation of the boundary layers : either, the
      !> tile will both compute some boundary layers and exchange
      !> the other ones or the tile will obtain the gridpoints of
      !> the boundary layers by exchanging the data with the
      !> neighbouring tiles
      !
      !> @date
      ! 22_08_2013  - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mpi_messenger_bc_ext_class

        use mpi_messenger_bc_class  , only : mpi_messenger_bc
        use mpi_mg_bc_ini_procedures, only : ini_bc_procedures
        use parameters_input        , only : nx,ny,ne

        implicit none

        private
        public :: mpi_messenger_bc_ext


        !> @class mpi_messenger_bc_ext
        !> class encapsulating attributes to identify the procedures
        !> for the computation of the boundary layers
        !>
        !> @param proc_x_choice
        !> integer identifying the procedure for the computation
        !> of boundary layers in the x direction: it will either
        !> compute the E and W boundary layers using the interior
        !> points, or compute one of them and exchange the other,
        !> or exchange both of them
        !>
        !> @param proc_y_choice
        !> integer identifying the procedure for the computation
        !> of boundary layers in the y direction: it will either
        !> compute the N and S boundary layers using the interior
        !> points, or compute one of them and exchange the other,
        !> or exchange both of them
        !>
        !> @param exchange_id
        !> if the procedure chosen for the boundary layers in the
        !> x-direction compute and exchange data, we need to know
        !> in which cardinal direction (N,S,E,W) the data are
        !> exchanged and in which direction the data are computed
        !>
        !> @param initialize
        !> subroutine initializing the attributes of the
        !> 'mpi_messenger_bc_ext' object
        !---------------------------------------------------------------
        type, extends(mpi_messenger_bc) :: mpi_messenger_bc_ext

          integer               :: proc_x_choice
          integer               :: proc_y_choice
          integer, dimension(2) :: exchange_id

          contains
          
          procedure, pass :: initialize

        end type mpi_messenger_bc_ext


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the attributes of
        !> the 'mpi_messenger_bc_ext' object
        !
        !> @date
        !> 21_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> 'mpi_messenger_bc' object initialized
        !
        !>@param f_used
        !> object containing the MPI cartesian communicator
        !
        !>@param s
        !> space discretisation method to know the boundary
        !> layer size
        !--------------------------------------------------------------
        subroutine initialize(this,f_used,s)

          implicit none

          class(mpi_messenger_bc_ext), intent(inout) :: this
          class(field_par)           , intent(in)    :: f_used
          type(cg_operators)         , intent(in)    :: s


          !< initialize the attributes of 'mpi_messenger_bc'
          call this%mpi_messenger_bc%initialize()

          
          !< initialize the attributes of 'mpi_messenger_bc_ext'
          call 



        end subroutine initialize


      end module mpi_messenger_bc_ext_class
