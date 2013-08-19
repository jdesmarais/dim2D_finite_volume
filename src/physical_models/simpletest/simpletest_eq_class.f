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
      module simpletest_eq_class
      
        use field_class        , only : field
        use parameters_constant, only : scalar
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind, rkind
        use phy_model_eq_class , only : phy_model_eq
        use cg_operators_class , only : cg_operators

        implicit none

        private
        public :: simpletest_eq


        !> @class simpletest_eq
        !> class encapsulating operators to compute
        !> simple test equations
        !>
        !> @param get_model_name
        !> get the name of the physcial model
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
        !---------------------------------------------------------------
        type, extends(phy_model_eq) :: simpletest_eq
          
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

        end type simpletest_eq


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

          character(len=32), dimension(ne) :: var_pties

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

          character(len=10), dimension(ne) :: var_pties

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
        subroutine apply_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used


          !<read the input file to know the user choice

          
          !<initialize the field depending on the user choice
          

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
        function compute_flux_x(field_used,s) result(flux_x)
        
          implicit none

          class(field)                      , intent(in) :: field_used
          type(cg_operators)                , intent(in) :: s
          real(rkind), dimension(nx+1,ny,ne)             :: flux_x

          integer :: i,j
          integer :: bc_size


          !<get the size of the boundary layers
          bc_size = s%get_bc_size()


          !<fluxes along the x-axis
          do j=bc_size+1, ny-bc_size
             do i=bc_size+1, nx+1-bc_size

                flux_x(i,j,1) = 10*s%f(field_used,i-1,j,basic)+
     $               s%dfdx(field_used,i-1,j,basic)

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
        function compute_flux_y(field_used,s) result(flux_y)
        
          implicit none

          class(field)                      , intent(in) :: field_used
          type(cg_operators)                , intent(in) :: s
          real(rkind), dimension(nx,ny+1,ne)             :: flux_y

          integer :: i,j
          integer :: bc_size


          !<get the size of the boundary layers
          bc_size = s%get_bc_size()


          !<fluxes along the x-axis
          do j=bc_size+1, ny+1-bc_size
             do i=bc_size+1, nx-bc_size

                flux_y(i,j,1) = s%g(field_used,i,j-1,basic)+
     $               10*s%dgdy(field_used,i,j-1,basic)

             end do
          end do

        end function compute_flux_y


        function basic(field_used,i,j) result(var)

          class(field)       , intent(in) :: field_used
          integer(ikind)     , intent(in) :: i
          integer(ikind)     , intent(in) :: j
          real(rkind)                     :: var

          var = field_used%nodes(i,j,1)

        end function basic

      end module simpletest_eq_class
