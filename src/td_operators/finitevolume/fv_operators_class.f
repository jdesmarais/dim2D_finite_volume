      !> @file
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> using the finite volume method
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module fv_operators_class

        use cg_operators_class, only : cg_operators
        use field_class       , only : field
        use parameters_kind   , only : rkind, ikind
        use phy_model_eq_class, only : phy_model_eq
        use td_operators_class, only : td_operators

        implicit none

        private
        public :: fv_operators


        !> @class fv_operators
        !> class encapsulating operators to compute
        !> the time derivatives of the main variables
        !> using the finite volume method
        !>
        !> @param compute_time_dev
        !> compute the time derivatives
        !---------------------------------------------------------------
        type, extends(td_operators) :: fv_operators

          contains

          procedure, nopass :: compute_time_dev

        end type fv_operators


        contains


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> subroutine computing the time derivatives using the
          !> space discretisation operators and the physical model
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p
          !> physical model
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          subroutine compute_time_dev(
     $       field_used,
     $       s,
     $       p_model,
     $       time_dev)

            implicit none


            class(field)                 , intent(in)   :: field_used
            type(cg_operators)           , intent(in)   :: s
            class(phy_model_eq)          , intent(in)   :: p_model
            real(rkind), dimension(:,:,:), intent(inout):: time_dev

            integer(ikind) :: nx
            integer(ikind) :: ny
            integer        :: ne
            integer        :: bc_size
            integer(ikind) :: i
            integer(ikind) :: j
            integer        :: k
            real(rkind), dimension(:,:,:), allocatable :: flux_x
            real(rkind), dimension(:,:,:), allocatable :: flux_y
            
            !<initialize the main tables size
            nx      = size(field_used%nodes,1)
            ny      = size(field_used%nodes,2)
            ne      = p_model%get_eq_nb()
            bc_size = s%get_bc_size()
            
            !<allocate the tables for the intermediate variables
            allocate(flux_x(nx+1,ny  ,ne))
            allocate(flux_y(nx  ,ny+1,ne))

            !<compute the fluxes
            !DEC$ FORCEINLINE RECURSIVE
            call p_model%compute_fluxes(field_used,s,flux_x,flux_y)

            !<compute the time derivatives
            do k=1, ne
               do j=1+bc_size, ny-bc_size
                  do i=1+bc_size, nx-bc_size
                     time_dev(i,j,k)=
     $                    (flux_x(i,j,k)-flux_x(i+1,j,k))/field_used%dx+
     $                    (flux_y(i,j,k)-flux_y(i,j+1,k))/field_used%dy
                  end do
               end do
            end do

            !<deallocate the intermediate variables
            deallocate(flux_x, flux_y)
            
        end subroutine compute_time_dev

      end module fv_operators_class
