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

        use bc_operators_class , only : bc_operators
        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use parameters_constant, only : earth_gravity_choice,
     $                                  bc_fluxes_choice
        use parameters_input   , only : nx,ny,ne,gravity_choice,
     $                                  bcx_type_choice, bcy_type_choice
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
          !>@param p_model
          !> physical model
          !
          !>@param bc_used
          !> boundary conditions
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          function compute_time_dev(field_used,s,p_model,bc_used)
     $       result(time_dev)

            implicit none


            class(field)                    , intent(in) :: field_used
            type(cg_operators)              , intent(in) :: s
            type(dim2d_eq)                  , intent(in) :: p_model
            type(bc_operators)              , intent(in) :: bc_used
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


            !<if the boundary conditions influence the computation
            !> of the fluxes, then we need to modify the fluxes
            if((bcx_type_choice.eq.bc_fluxes_choice).or.
     $         (bcy_type_choice.eq.bc_fluxes_choice)) then
               call bc_used%apply_bc_on_fluxes(
     $              field_used,s,flux_x,flux_y)
            end if


            !<compute the time derivatives
            !>select if the body forces computation is required
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
