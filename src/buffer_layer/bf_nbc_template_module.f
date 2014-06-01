      module bf_nbc_template_module

        use parameters_bf_layer, only: bc_pt,
     $                                 bc_interior_pt,
     $                                 interior_pt

        implicit none

        private
        public :: make_nbc_template_11,
     $            make_nbc_template_12,
     $            make_nbc_template_13,
     $            make_nbc_template_14,
     $            make_nbc_template_15,
     $            make_nbc_template_21,
     $            make_nbc_template_22,
     $            make_nbc_template_23,
     $            make_nbc_template_24,
     $            make_nbc_template_25,
     $            make_nbc_template_31,
     $            make_nbc_template_32,
     $            make_nbc_template_34,
     $            make_nbc_template_35,
     $            make_nbc_template_41,
     $            make_nbc_template_42,
     $            make_nbc_template_43,
     $            make_nbc_template_44,
     $            make_nbc_template_45,
     $            make_nbc_template_51,
     $            make_nbc_template_52,
     $            make_nbc_template_53,
     $            make_nbc_template_54,
     $            make_nbc_template_55

        contains

        
        function make_nbc_template_11()

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


        function make_nbc_template_21()

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


        function make_nbc_template_31()

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


        function make_nbc_template_41()

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


        function make_nbc_template_51()

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


        function make_nbc_template_12()

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


        function make_nbc_template_22()

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


        function make_nbc_template_32()

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


        function make_nbc_template_42()

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


        function make_nbc_template_52()

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


        function make_nbc_template_13()

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


        function make_nbc_template_23()

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


        function make_nbc_template_43()

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


        function make_nbc_template_53()

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


        function make_nbc_template_14()

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


        function make_nbc_template_24()

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


        function make_nbc_template_34()

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


        function make_nbc_template_44()

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


        function make_nbc_template_54()

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

      
        function make_nbc_template_15()

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


        function make_nbc_template_25()

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


        function make_nbc_template_35()

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


        function make_nbc_template_45()

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


        function make_nbc_template_55()

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

      end module bf_nbc_template_module
