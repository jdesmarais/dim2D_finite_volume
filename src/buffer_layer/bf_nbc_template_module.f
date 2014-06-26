      !> @file
      !> module implementing the initializaion of temporary arrays
      !> when asking the role of grid points at the edge between the
      !> interior and the buffer layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the initializaion of temporary arrays
      !> when asking the role of grid points at the edge between the
      !> interior and the buffer layers
      !> \image html  bf_nbc_template_module.png
      !> \image latex bf_nbc_template_module.eps width=10cm
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module bf_nbc_template_module

        use parameters_bf_layer, only: no_pt,
     $                                 bc_pt,
     $                                 bc_interior_pt,
     $                                 interior_pt

        implicit none

        private
        public :: make_nbc_template_00,
     $            make_nbc_template_01,
     $            make_nbc_template_02,
     $            make_nbc_template_05,
     $            make_nbc_template_06,
     $            make_nbc_template_10,
     $            make_nbc_template_11,
     $            make_nbc_template_12,
     $            make_nbc_template_13,
     $            make_nbc_template_14,
     $            make_nbc_template_15,
     $            make_nbc_template_16,
     $            make_nbc_template_20,
     $            make_nbc_template_21,
     $            make_nbc_template_22,
     $            make_nbc_template_23,
     $            make_nbc_template_24,
     $            make_nbc_template_25,
     $            make_nbc_template_26,
     $            make_nbc_template_31,
     $            make_nbc_template_32,
     $            make_nbc_template_34,
     $            make_nbc_template_35,
     $            make_nbc_template_41,
     $            make_nbc_template_42,
     $            make_nbc_template_43,
     $            make_nbc_template_44,
     $            make_nbc_template_45,
     $            make_nbc_template_50,
     $            make_nbc_template_51,
     $            make_nbc_template_52,
     $            make_nbc_template_53,
     $            make_nbc_template_54,
     $            make_nbc_template_55,
     $            make_nbc_template_56,
     $            make_nbc_template_60,
     $            make_nbc_template_61,
     $            make_nbc_template_62,
     $            make_nbc_template_65,
     $            make_nbc_template_66

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_00() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = no_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = no_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_00


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_10() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = no_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_10


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_20() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = no_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_20

      !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_50() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = no_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_50


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_60() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = no_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = no_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_60


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_01() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = no_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_01


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_11() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_11


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_21() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_21


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_31() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_31


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_41() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_41


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_51() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_51


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_61() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = no_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_61


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_02() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = no_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_02


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_12() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_12


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_22() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_22


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_32() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_32


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_42() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_42


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_52() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_52


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_62() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = no_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_62


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_13() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_13


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_23() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = interior_pt

        end function make_nbc_template_23


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_43() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_43


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_53() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_53


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_14() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_14


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_24() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_24


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_34() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_34


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_44() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_interior_pt

        end function make_nbc_template_44


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_54() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = bc_interior_pt
          nbc_template(2,3) = bc_interior_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_54

      !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_05() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = no_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_05


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_15() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_15


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_25() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_25


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_35() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_35


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_45() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_interior_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_45


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_55() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_interior_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = bc_pt

        end function make_nbc_template_55


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_65() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_interior_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = no_pt
          
          nbc_template(1,3) = bc_pt
          nbc_template(2,3) = bc_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_65


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_06() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = no_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = no_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = no_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_06


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_16() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = no_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_16


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_26() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_interior_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = no_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_26


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_56() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_interior_pt
          nbc_template(3,1) = bc_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = bc_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = no_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_56


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the temporary array giving the role of the grid
        !> points at the edge between the interior domain and the buffer
        !> layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@return nbc_template
        !> temporary array with the role of the grid point
        !--------------------------------------------------------------
        function make_nbc_template_66() result(nbc_template)

          implicit none
          
          integer, dimension(3,3) :: nbc_template

          nbc_template(1,1) = bc_interior_pt
          nbc_template(2,1) = bc_pt
          nbc_template(3,1) = no_pt
          
          nbc_template(1,2) = bc_pt
          nbc_template(2,2) = bc_pt
          nbc_template(3,2) = no_pt
          
          nbc_template(1,3) = no_pt
          nbc_template(2,3) = no_pt
          nbc_template(3,3) = no_pt

        end function make_nbc_template_66

      end module bf_nbc_template_module
