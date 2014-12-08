      !> @file
      !> computation of gradients in n1- and n2- directions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the gradients in the n1- and n2- directions
      !> depending on the number of grid points available
      !
      !> @date
      !> 17_11_2014 - Initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module sd_operators_fd_n_module

        use interface_primary, only :
     $       get_primary_var,
     $       gradient_x_proc,
     $       gradient_y_proc

        use parameters_input, only :
     $       ne

        use parameters_kind  , only :
     $       ikind,
     $       rkind

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0


        implicit none

        private
        public ::
     $       gradient_n1_xL0_yL0,
     $       gradient_n1_xL0_yL1,
     $       gradient_n1_xL0_yI,
     $       gradient_n1_xL0_yR1,
     $       gradient_n1_xL0_yR0,
     $       
     $       gradient_n1_xL1_yL0,
     $       gradient_n1_xL1_yL1,
     $       gradient_n1_xL1_yI,
     $       gradient_n1_xL1_yR1,
     $       gradient_n1_xL1_yR0,
     $       
     $       gradient_n1_xI_yL0,
     $       gradient_n1_xI_yL1,
     $       gradient_n1_xI_yI,
     $       gradient_n1_xI_yR1,
     $       gradient_n1_xI_yR0,
     $       
     $       gradient_n1_xR1_yL0,
     $       gradient_n1_xR1_yL1,
     $       gradient_n1_xR1_yI,
     $       gradient_n1_xR1_yR1,
     $       gradient_n1_xR1_yR0,
     $       
     $       gradient_n1_xR0_yL0,
     $       gradient_n1_xR0_yL1,
     $       gradient_n1_xR0_yI,
     $       gradient_n1_xR0_yR1,
     $       gradient_n1_xR0_yR0,
     $       
     $       gradient_n2_xL0_yL0,
     $       gradient_n2_xL0_yL1,
     $       gradient_n2_xL0_yI,
     $       gradient_n2_xL0_yR1,
     $       gradient_n2_xL0_yR0,
     $       
     $       gradient_n2_xL1_yL0,
     $       gradient_n2_xL1_yL1,
     $       gradient_n2_xL1_yI,
     $       gradient_n2_xL1_yR1,
     $       gradient_n2_xL1_yR0,
     $       
     $       gradient_n2_xI_yL0,
     $       gradient_n2_xI_yL1,
     $       gradient_n2_xI_yI,
     $       gradient_n2_xI_yR1,
     $       gradient_n2_xI_yR0,
     $       
     $       gradient_n2_xR1_yL0,
     $       gradient_n2_xR1_yL1,
     $       gradient_n2_xR1_yI,
     $       gradient_n2_xR1_yR1,
     $       gradient_n2_xR1_yR0,
     $       
     $       gradient_n2_xR0_yL0,
     $       gradient_n2_xR0_yL1,
     $       gradient_n2_xR0_yI,
     $       gradient_n2_xR0_yR1,
     $       gradient_n2_xR0_yR0,
     $       
     $       get_gradient_n1,
     $       get_gradient_n2
        
        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL0_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n1_xL0_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL0_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n1_xL0_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL0_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n1_xL0_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL0_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n1_xL0_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL0_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n1_xL0_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL1_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n1_xL1_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL1_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n1_xL1_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL1_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n1_xL1_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL1_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n1_xL1_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xL1_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n1_xL1_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and
        !> no-grid points under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xI_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n1_xI_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and one
        !> grid point under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xI_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n1_xI_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and one
        !> grid point under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xI_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n1_xI_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and no
        !> grid point above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xI_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n1_xI_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and no
        !> grid point above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xI_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n1_xI_yR0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR1_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n1_xR1_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR1_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n1_xR1_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR1_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n1_xR1_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR1_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n1_xR1_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR1_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n1_xR1_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR0_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n1_xR0_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR0_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n1_xR0_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR0_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n1_xR0_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR0_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n1_xR0_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_xR0_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n1(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n1_xR0_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL0_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n2_xL0_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL0_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n2_xL0_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL0_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n2_xL0_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL0_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n2_xL0_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no-grid points on the left side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL0_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n2_xL0_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL1_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n2_xL1_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL1_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n2_xL1_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL1_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n2_xL1_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL1_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n2_xL1_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the left side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xL1_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_L1,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n2_xL1_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and
        !> no-grid points under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xI_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n2_xI_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and one
        !> grid point under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xI_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n2_xI_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and one
        !> grid point under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xI_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n2_xI_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and no
        !> grid point above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xI_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n2_xI_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with two grid points on the left and right sides and no
        !> grid point above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xI_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n2_xI_yR0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR1_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n2_xR1_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR1_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n2_xR1_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR1_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n2_xR1_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR1_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n2_xR1_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with one grid point on the right side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR1_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R1,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n2_xR1_yR0      


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and no-grid points
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR0_yL0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L0,
     $         dx,dy)

        end function gradient_n2_xR0_yL0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and one grid point
        !> under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR0_yL1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L1,
     $         dx,dy)

        end function gradient_n2_xR0_yL1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and two grid points
        !> above and under
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR0_yI(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         dx,dy)

        end function gradient_n2_xR0_yI


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and one grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR0_yR1(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R1,
     $         dx,dy)

        end function gradient_n2_xR0_yR1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !> with no grid point on the right side and no grid point
        !> above
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_xR0_yR0(nodes,i,j,proc,dx,dy)
     $       result(var)

          implicit none
          
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var =  gradient_n2(
     $         nodes,
     $         i,j,
     $         proc,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0,
     $         dx,dy)

        end function gradient_n2_xR0_yR0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left(\frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param grad_x_proc
        !> procedure computing the gradient along the x-direction
        !> at [i,j]
        !
        !>@param grad_y_proc
        !> procedure computing the gradient along the y-direction
        !> at [i,j]
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1(
     $       nodes,
     $       i,j,
     $       var_proc,
     $       grad_x_proc,
     $       grad_y_proc,
     $       dx,dy)
     $       result(var)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          procedure(get_primary_var)                   :: var_proc
          procedure(gradient_x_proc)                   :: grad_x_proc
          procedure(gradient_y_proc)                   :: grad_y_proc
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          real(rkind)                                  :: var

          
          !DEC$ FORCEINLINE RECURSIVE
          var = 0.5d0*SQRT(2.0d0)*(
     $          grad_x_proc(nodes,i,j,var_proc,dx)
     $        - grad_y_proc(nodes,i,j,var_proc,dy))

        end function gradient_n1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left( \frac{\partial u}{\partial x} +
        !> \frac{\partial u}{\partial y} \right)
        !
        !> @date
        !> 17_11_2014 - initial version  - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !
        !>@param i
        !> index along x-axis where the data is evaluated
        !
        !>@param j
        !> index along y-axis where the data is evaluated
        !
        !>@param var_proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param grad_x_proc
        !> procedure computing the gradient along the x-direction
        !> at [i,j]
        !
        !>@param grad_y_proc
        !> procedure computing the gradient along the y-direction
        !> at [i,j]
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2(
     $       nodes,
     $       i,j,
     $       var_proc,
     $       grad_x_proc,
     $       grad_y_proc,
     $       dx,dy)
     $       result(var)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: nodes
          integer(ikind)                  , intent(in) :: i
          integer(ikind)                  , intent(in) :: j
          procedure(get_primary_var)                   :: var_proc
          procedure(gradient_x_proc)                   :: grad_x_proc
          procedure(gradient_y_proc)                   :: grad_y_proc
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          real(rkind)                                  :: var

          
          !DEC$ FORCEINLINE RECURSIVE
          var = 0.5d0*SQRT(2.0d0)*(
     $          grad_x_proc(nodes,i,j,var_proc,dx)
     $        + grad_y_proc(nodes,i,j,var_proc,dy))

        end function gradient_n2


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left( \frac{\partial u}{\partial x} -
        !> \frac{\partial u}{\partial y} \right)
        !
        !> @date
        !> 08_12_2014 - initial version  - J.L. Desmarais
        !
        !>@param var_gradient_x
        !> gradient of the governing variables in the x-direction
        !
        !>@param var_gradient_y
        !> gradient of the governing variables in the y-direction
        !
        !>@param var_gradient_n1
        !> gradient of the governing variables in the n1-direction
        !---------------------------------------------------------------
        function get_gradient_n1(var_gradient_x,var_gradient_y)
     $     result(var_gradient_n1)

          implicit none

          real(rkind), dimension(ne), intent(in)  :: var_gradient_x
          real(rkind), dimension(ne), intent(in)  :: var_gradient_y
          real(rkind), dimension(ne)              :: var_gradient_n1

          integer :: k

          do k=1,ne
             
             var_gradient_n1(k) = 0.5d0*SQRT(2.0d0)*(
     $            var_gradient_x(k)
     $           -var_gradient_y(k))

          end do

        end function get_gradient_n1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\big|_{i,j}=
        !> \frac{\sqrt{2}}{2} \left( \frac{\partial u}{\partial x} +
        !> \frac{\partial u}{\partial y} \right)
        !
        !> @date
        !> 08_12_2014 - initial version  - J.L. Desmarais
        !
        !>@param var_gradient_x
        !> gradient of the governing variables in the x-direction
        !
        !>@param var_gradient_y
        !> gradient of the governing variables in the y-direction
        !
        !>@param var_gradient_n2
        !> gradient of the governing variables in the n2-direction
        !---------------------------------------------------------------
        function get_gradient_n2(var_gradient_x,var_gradient_y)
     $     result(var_gradient_n2)

          implicit none

          real(rkind), dimension(ne), intent(in)  :: var_gradient_x
          real(rkind), dimension(ne), intent(in)  :: var_gradient_y
          real(rkind), dimension(ne)              :: var_gradient_n2

          integer :: k

          do k=1,ne
             
             var_gradient_n2(k) = 0.5d0*SQRT(2.0d0)*(
     $            var_gradient_x(k)
     $           +var_gradient_y(k))

          end do

        end function get_gradient_n2

      end module sd_operators_fd_n_module
