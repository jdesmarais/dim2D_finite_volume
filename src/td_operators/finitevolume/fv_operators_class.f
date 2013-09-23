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

        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use parameters_constant, only : earth_gravity_choice
        use parameters_input   , only : nx,ny,ne,gravity_choice
        use parameters_kind    , only : rkind, ikind
        use dim2d_eq_class     , only : dim2d_eq
        use td_operators_class , only : td_operators

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
          function compute_time_dev(field_used,s,p_model)result(time_dev)

            implicit none


            class(field)                    , intent(in) :: field_used
            type(cg_operators)              , intent(in) :: s
            type(dim2d_eq)                  , intent(in) :: p_model
            real(rkind), dimension(nx,ny,ne)             :: time_dev

            integer                            :: bc_size,k
            integer(ikind)                     :: i,j
            real(rkind), dimension(nx+1,ny,ne) :: flux_x
            real(rkind), dimension(nx,ny+1,ne) :: flux_y
            real(rkind), dimension(nx,ny,ne)   :: body_forces

            !<initialize the main tables size
            bc_size = s%get_bc_size()
            
            !<compute the fluxes
            !FORCEINLINE RECURSIVE
            flux_x = p_model%compute_flux_x(field_used,s)

            !FORCEINLINE RECURSIVE
            flux_y = p_model%compute_flux_y(field_used,s)


            !<compute the time derivatives depending on body forces or not
            if(gravity_choice.eq.earth_gravity_choice) then

               !<compute the body forces
               body_forces = p_model%compute_body_forces(field_used,s)

               !<compute the time derivatives
               do k=1, ne
                  do j=1+bc_size, ny-bc_size
                     do i=1+bc_size, nx-bc_size
                        time_dev(i,j,k)=
     $                       (flux_x(i,j,k)-flux_x(i+1,j,k))/field_used%dx+
     $                       (flux_y(i,j,k)-flux_y(i,j+1,k))/field_used%dy+
     $                       body_forces(i,j,k)
                     end do
                  end do
               end do

            else

               !<compute the time derivatives
               do k=1, ne
                  do j=1+bc_size, ny-bc_size
                     do i=1+bc_size, nx-bc_size
                        time_dev(i,j,k)=
     $                       (flux_x(i,j,k)-flux_x(i+1,j,k))/field_used%dx+
     $                       (flux_y(i,j,k)-flux_y(i,j+1,k))/field_used%dy
                     end do
                  end do
               end do

            end if

        end function compute_time_dev

      end module fv_operators_class
