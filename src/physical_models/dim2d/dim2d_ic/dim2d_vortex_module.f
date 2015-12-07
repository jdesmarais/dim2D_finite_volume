      !> @file
      !> module encapsulating subroutines to compute
      !> the velocity field due to vortices
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the velocity field due to vortices
      !
      !> @date
      !> 27_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_vortex_module

        use parameters_kind, only: rkind

        implicit none


        private
        public :: get_vortex_velocity


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function computing the velocity due to a vortex
        !> whose center, spread and rotational are prescribed
        !> \f[
        !> \begin{pmatrix} u(x,y) \\\ v(x,y) \end{pmatrix} = 
        !> \begin{pmatrix} - \displaystyle{\frac{y-y_c}{r}} \\\ \\\ \displaystyle{\frac{x-x_c}{r}} \end{pmatrix} w(r) 
        !> \f]
        !> where 
        !> \f[ r = \sqrt{(x-x_c)^2+(y-y_c)^2} \f]
        !> \f[ w(r) = \begin{cases}
        !>  \omega r & \mbox{if } r<r_c\\\
        !>  \omega \displaystyle{\frac{r_c^2}{r}} & \mbox{if } r \ge r_c
        !> \end{cases}
        !> \f]
        !> \f[ \omega = \textrm{rot} \begin{pmatrix} u \\\ v \end{pmatrix} \f]
        !
        !> @date
        !> 27_09_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> spatial coordinate along the x-axis
        !
        !>@param y
        !> spatial coordinate along the y-axis
        !
        !>@param x_c
        !> spatial coordinate along the x-axis for the vortex center
        !
        !>@param y_c
        !> spatial coordinate along the y-axis for the vortex center
        !
        !>@param r_c
        !> spread of the vortex
        !
        !>@param omega
        !> rotational of the velocity field: strength of the vortex
        !
        !>@return
        !> velocity, \f$ (u,v)^T \f$
        !---------------------------------------------------------------
        function get_vortex_velocity(x,y,x_c,y_c,r_c, omega)
     $       result(velocity)

          implicit none

          real(rkind), intent(in)   :: x
          real(rkind), intent(in)   :: y
          real(rkind), intent(in)   :: x_c
          real(rkind), intent(in)   :: y_c
          real(rkind), intent(in)   :: r_c
          real(rkind), intent(in)   :: omega
          real(rkind), dimension(2) :: velocity

          real(rkind) :: r
          real(rkind) :: v_theta
          real(rkind) :: cos_theta
          real(rkind) :: sin_theta

          !< compute the radial distance from the
          !> vortex center
          r = sqrt((x-x_c)**2+(y-y_c)**2)

          !< compute the orthoradial velocity
          !> of the vortex
          if(r.lt.r_c) then
             v_theta = r*omega
          else
             v_theta = omega*r_c**2/r
          end if

          !< the cos_theta and sin_theta are used
          !> to compute the projection of the
          !> velocity field in cylindrical
          !> coordinates into cartesian coordinates
          !> assuming that the radial component of
          !> the velocity field is null
          cos_theta = (x-x_c)/r
          sin_theta = (y-y_c)/r

          !< projection of the velocity field in
          !> cylindrical into cartesian coordinates
          velocity(1) = - sin_theta*v_theta
          velocity(2) =   cos_theta*v_theta

        end function get_vortex_velocity

      end module dim2d_vortex_module
