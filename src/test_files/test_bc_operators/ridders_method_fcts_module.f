      module ridders_method_fcts_module

        use ridders_method_module, only :
     $       root_fct_abstract

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public ::
     $       root_fct1,
     $       root_fct2,
     $       root_fct3,
     $       root_fct4
        
        type, abstract, extends(root_fct_abstract) :: root_fct

          integer :: s

        end type root_fct


        type, extends(root_fct) :: root_fct1

          contains

          procedure, pass :: f => fct1

        end type root_fct1

        type, extends(root_fct) :: root_fct2

          contains

          procedure, pass :: f => fct2

        end type root_fct2

        type, extends(root_fct) :: root_fct3

          contains

          procedure, pass :: f => fct3

        end type root_fct3

        type, extends(root_fct) :: root_fct4

          contains

          procedure, pass :: f => fct4

        end type root_fct4


        contains


        function fct1(this,x) result(fx)

          implicit none

          class(root_fct1), intent(in) :: this
          real(rkind)     , intent(in) :: x
          real(rkind)                  :: fx

          integer :: s

          fx = 1.0d0+x

          s=this%s

        end function fct1


        function fct2(this,x) result(fx)

          implicit none

          class(root_fct2), intent(in) :: this
          real(rkind)     , intent(in) :: x
          real(rkind)                  :: fx

          integer :: s

          fx = (x-1.0d0)*(x-2.0d0)

          s = this%s

        end function fct2


        function fct3(this,x) result(fx)

          implicit none

          class(root_fct3), intent(in) :: this
          real(rkind)     , intent(in) :: x
          real(rkind)                  :: fx

          integer :: s
          
          fx = Tanh(x)

          s=this%s

        end function fct3


        function fct4(this,x) result(fx)

          implicit none

          class(root_fct4), intent(in) :: this
          real(rkind)     , intent(in) :: x
          real(rkind)                  :: fx

          integer :: s

          fx = Tanh(x-2.0d0)

          s=this%s

        end function fct4


      end module ridders_method_fcts_module
