      !> @file
      !> class encapsulating subroutines to compute
      !> some simple test equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> some simple test equations
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module pmodel_eq_class
      
        use parameters_bf_layer     , only : interior_pt
        use parameters_constant     , only : scalar
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind
        use pmodel_eq_abstract_class, only : pmodel_eq_abstract
        use sd_operators_class      , only : sd_operators

        implicit none

        private
        public :: pmodel_eq


        !> @class pmodel_eq
        !> class encapsulating operators to compute
        !> simple test equations
        !>
        !> @param get_model_name
        !> get the name of the physical model
        !>
        !> @param get_var_name
        !> get the name of the main variables
        !> (mass)
        !>
        !> @param get_var_longname
        !> get the description of the main variables for the
        !> governing equations of the physical model
        !>
        !> @param get_var_unit
        !> get the units of the main variables
        !>
        !> @param get_var_types
        !> get the type of the main variables
        !> (scalar)
        !>
        !> @param get_eq_nb
        !> get the number of governing equations: 1
        !>
        !> @param apply_initial_conditions
        !> initialize the main variables
        !>
        !> @param compute_fluxes
        !> compute the fluxes along the x- and y-axis
        !>
        !> @param are_openbc_undermined
        !> check whether the open boundary conditions are undermined
        !> at the grid point location
        !---------------------------------------------------------------
        type, extends(pmodel_eq_abstract) :: pmodel_eq
          
          contains

          procedure, nopass :: get_model_name
          procedure, nopass :: get_var_name
          procedure, nopass :: get_var_longname
          procedure, nopass :: get_var_unit
          procedure, nopass :: get_var_type
          procedure, nopass :: get_eq_nb
          procedure, nopass :: apply_ic
          procedure, nopass :: compute_flux_x
          procedure, nopass :: compute_flux_y
          procedure, nopass :: compute_flux_x_nopt
          procedure, nopass :: compute_flux_y_nopt
          procedure, nopass :: compute_flux_x_oneside
          procedure, nopass :: compute_flux_y_oneside
          procedure, nopass :: compute_body_forces
          procedure, nopass :: get_velocity
          procedure, nopass :: are_openbc_undermined

        end type pmodel_eq


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the name of the physical model
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param model_name
        !> character giving the name of the model
        !---------------------------------------------------------------
        function get_model_name() result(model_name)

          implicit none

          character(len=10) :: model_name

          model_name="Simpletest"

        end function get_model_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_name() result(var_pties)

          implicit none

          character(len=10), dimension(ne) :: var_pties

          var_pties(1)="mass"

        end function get_var_name
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the name of the main variables
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable names
        !---------------------------------------------------------------
        function get_var_longname() result(var_pties)

          implicit none

          character(len=33), dimension(ne) :: var_pties

          var_pties(1)="mass density"

        end function get_var_longname


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the units of the main variables
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable units
        !---------------------------------------------------------------
        function get_var_unit() result(var_pties)

          implicit none

          character(len=23), dimension(ne) :: var_pties

          var_pties(1)= "(kg/m3)/(kg/m3)"

        end function get_var_unit


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to get the type of the main variables
        !> (scalar, vector_x, vector_y, scalar)
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param var_name
        !> characters giving the variable type
        !---------------------------------------------------------------
        function get_var_type() result(var_type)

          implicit none

          integer, dimension(ne) :: var_type

          var_type(1)=scalar

        end function get_var_type
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the number of main variables
        !> in the governing equations: 4
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param eq_nb
        !> number of governing equations
        !---------------------------------------------------------------
        function get_eq_nb() result(eq_nb)
          implicit none
          integer :: eq_nb
          eq_nb=1
        end function get_eq_nb
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions to the main
        !> variables of the governing equations
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map

          real(rkind) :: node_s,x_s,y_s

          node_s = nodes(1,1,1)
          x_s = x_map(1)
          y_s = y_map(1)

        end subroutine apply_ic
        
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param sd_operators_used
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-axis
        !---------------------------------------------------------------
        function compute_flux_x(nodes,dx,dy,s) result(flux_x)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          type(sd_operators)                , intent(in)   :: s
          real(rkind), dimension(nx+1,ny,ne)               :: flux_x

          integer     :: i,j
          real(rkind) :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
          do j=bc_size+1, ny-bc_size
             do i=bc_size+1, nx+1-bc_size

c$$$                flux_x(i,j,1) = 10*s%f(nodes,i,j,basic)+
c$$$     $               s%dfdx(nodes,i,j,basic,dx)

                flux_x(i,j,1) = s%f(nodes,i,j,basic)

             end do
          end do

        end function compute_flux_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to apply the initial conditions
        !> to the main variables of the governing
        !> equations
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_y
        !> fluxes along the y-axis
        !---------------------------------------------------------------
        function compute_flux_y(nodes,dx,dy,s) result(flux_y)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          type(sd_operators)                , intent(in)   :: s
          real(rkind), dimension(nx,ny+1,ne)               :: flux_y

          integer     :: i,j
          real(rkind) :: dx_s

          dx_s = dx

          !<fluxes along the x-axis
          do j=bc_size+1, ny+1-bc_size
             do i=bc_size+1, nx-bc_size

c$$$                flux_y(i,j,1) = s%g(nodes,i,j,basic)+
c$$$     $               10*s%dgdy(nodes,i,j,basic,dy)

                flux_y(i,j,1) = s%g(nodes,i,j,basic)

             end do
          end do

        end function compute_flux_y


        subroutine compute_flux_x_nopt(nodes,dx,dy,s,grdpts_id,flux_x)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: s
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: flux_x

          integer(ikind) :: i,j
          real(rkind)    :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
          do j=1+bc_size, size(flux_x,2)-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, size(flux_x,1)-bc_size

                if(grdpts_id(i,j).eq.interior_pt) then

c$$$                   flux_x(i,j,1) = 10*s%f(nodes,i,j,basic)+
c$$$     $               s%dfdx(nodes,i,j,basic,dx)

                   flux_x(i,j,1) = s%f(nodes,i,j,basic)

                end if

             end do
          end do

        end subroutine compute_flux_x_nopt


        subroutine compute_flux_y_nopt(nodes,dx,dy,s,grdpts_id,flux_y)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: s
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: flux_y

          integer(ikind) :: i,j
          real(rkind)    :: dx_s

          dx_s = dx

          !<fluxes along the y-axis
          do j=1+bc_size, size(flux_y,2)-bc_size
             !DEC$ IVDEP
             do i=1+bc_size, size(flux_y,1)-bc_size

                if(grdpts_id(i,j).eq.interior_pt) then

c$$$                   flux_y(i,j,1) = s%g(nodes,i,j,basic)+
c$$$     $                  10*s%dgdy(nodes,i,j,basic,dy)

                   flux_y(i,j,1) = s%g(nodes,i,j,basic)

                end if

             end do
          end do

        end subroutine compute_flux_y_nopt


        function compute_flux_x_oneside(nodes,dx,dy,i,j,s_oneside)
     $     result(flux_x)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          integer(ikind)                    , intent(in)   :: i
          integer(ikind)                    , intent(in)   :: j
          class(sd_operators)               , intent(in)   :: s_oneside
          real(rkind), dimension(ne)                       :: flux_x

          real(rkind) :: dy_s

          dy_s = dy

          !<fluxes along the x-axis
c$$$          flux_x(1) = 10*s_oneside%f(nodes,i,j,basic)+
c$$$     $               s_oneside%dfdx(nodes,i,j,basic,dx)

          flux_x(1) = s_oneside%f(nodes,i,j,basic)

        end function compute_flux_x_oneside


        function compute_flux_y_oneside(nodes,dx,dy,i,j,s_oneside)
     $     result(flux_y)
        
          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)   :: nodes
          real(rkind)                       , intent(in)   :: dx
          real(rkind)                       , intent(in)   :: dy
          integer(ikind)                    , intent(in)   :: i
          integer(ikind)                    , intent(in)   :: j
          class(sd_operators)                , intent(in)  :: s_oneside
          real(rkind), dimension(ne)                       :: flux_y

          real(rkind) :: dx_s

          dx_s = dx


          !<fluxes along the x-axis
c$$$          flux_y(1) = s_oneside%g(nodes,i,j,basic)+
c$$$     $         10*s_oneside%dgdy(nodes,i,j,basic,dy)

          flux_y(1) = s_oneside%g(nodes,i,j,basic)

        end function compute_flux_y_oneside


        function basic(nodes,i,j) result(var)

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = nodes(i,j,1)

        end function basic


        function compute_body_forces(nodes,k) result(body_forces)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          integer                   , intent(in) :: k
          real(rkind)                            :: body_forces

          real(rkind) :: node_s
          integer :: k_s

          body_forces = 0

          node_s = nodes(1)
          k_s = k

        end function compute_body_forces


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> interface to compute the body forces
        !> acting on the cell
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> governing variables at the grid point location
        !
        !>@param velocity
        !> velocity vector at the grid point location
        !--------------------------------------------------------------
        function get_velocity(nodes) result(velocity)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          real(rkind), dimension(2)              :: velocity

          velocity(1) = nodes(1)
          velocity(2) = nodes(1)

        end function get_velocity


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the open boundary conditions
        !> are undermined at the grid point location
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param undermined
        !> check if the open boundary conditions are undermined
        !> at the grid point location
        !--------------------------------------------------------------
        function are_openbc_undermined(nodes) result(undermined)

          implicit none

          real(rkind), dimension(ne), intent(in) :: nodes
          logical                                :: undermined

          real(rkind) :: d_liq, d_vap

          d_liq = 1.1-0.05*(1.1-0.1)
          d_vap = 0.1+0.05*(1.1-0.1)

          if((nodes(1).ge.d_vap).and.(nodes(1).le.d_liq)) then
             undermined = .true.
          else
             undermined = .false.
          end if

          !undermined = .true.

        end function are_openbc_undermined

      end module pmodel_eq_class
