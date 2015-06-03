      !> @file
      !> Ridder's method
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> implementation of Ridder's method for the computation of
      !> the root of a non-linear equation
      !
      !> @date
      !> 02_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ridders_method_module

        use check_data_module, only :
     $     is_real_validated

        use parameters_kind, only :
     $     ikind,
     $     rkind

        implicit none


        private
        public ::
     $       root_fct_abstract,
     $       get_root_ridder_method


        !> @class root_fct_abstract
        !> abstract class encapsulating the function
        !> whose root is searched by Ridder's method
        !
        !> @param f
        !> function whose root is searched by Ridder's
        !> method
        !------------------------------------------------------------
        type, abstract :: root_fct_abstract

          contains

          procedure(f_proc), pass, deferred :: f

        end type root_fct_abstract

        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for the function whose root is 
          !> determined by Ridder's method
          !
          !> @date
          !> 02_06_2015 - initial version - J.L. Desmarais
          !
          !>@param x
          !> parameter for the non-linear equation
          !>
          !>@param fx
          !> evaluation of the non-linear equation at x
          !--------------------------------------------------------------
          function f_proc(this,x) result(fx)

            import root_fct_abstract
            import rkind

            class(root_fct_abstract), intent(in) :: this
            real(rkind)             , intent(in) :: x
            real(rkind)                          :: fx

          end function f_proc

        end interface

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> using Ridder's method, return the root of a function
        !> \f$func\f$ know to lie between \f$x1\f$ and \f$x2\f$.
        !> The root, returned as \f$x\f$, will be refined to an
        !> approximate accuracy \f$xacc\f$.
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param func
        !> object containing the non-linear function whose root is
        !> determined. As the non-linear function may need more
        !> arguments than the parameter determining the root, it is
        !> encapsulated inside an object where the other parameters are
        !> defined as attributes
        !
        !>@param x1
        !> lower bracket for the root
        !
        !>@param x2
        !> upper bracket for the root
        !
        !>@param xacc
        !> accuracy
        !
        !>@param x
        !> approximated root
        !--------------------------------------------------------------
        function get_root_ridder_method(func,x1,x2,xacc) result(x)

          implicit none
          
          class(root_fct_abstract), intent(in) :: func
          real(rkind)             , intent(in) :: x1
          real(rkind)             , intent(in) :: x2
          real(rkind)             , intent(in) :: xacc
          real(rkind)                          :: x

          integer    , parameter :: MAXIT=60
          real(rkind), parameter :: UNUSED=-1.11e30

          integer     :: j
          real(rkind) :: fh
          real(rkind) :: fl
          real(rkind) :: fm
          real(rkind) :: fnew
          real(rkind) :: s
          real(rkind) :: xh
          real(rkind) :: xl
          real(rkind) :: xm
          real(rkind) :: xnew

          
          fl = func%f(x1)
          fh = func%f(x2)

          if( ((fl.gt.0.0).and.(fh.lt.0.0)) .or. ((fl.lt.0.0).and.(fh.gt.0.0)) ) then

             xl = x1
             xh = x2
             x  = UNUSED

             do j=1, MAXIT

                xm = 0.5d0*(xl+xh)
                fm = func%f(xm)
                s  = SQRT(fm**2-fl*fh)

                if(is_real_validated(s,0.0d0,.false.)) return

                xnew = xm+(xm-xl)*(sign(1.0d0,fl-fh)*fm/s)

                if(abs(xnew-x).le.xacc) return

                x = xnew
                fnew=func%f(x)

                if(is_real_validated(fnew,0.0d0,.false.)) return

                if(.not.is_real_validated(sign(fm,fnew),fm,.false.)) then
                   xl = xm
                   fl = fm
                   xh = x
                   fh = fnew

                else if (.not.is_real_validated(sign(fl,fnew),fl,.false.)) then
                   xh = x
                   fh = fnew

                else if (.not.is_real_validated(sign(fh,fnew),fh,.false.)) then
                   xl = x
                   fl = fnew

                else
                   print '(''ridder_method_module'')'
                   print '(''get_root_ridder_method'')'
                   print '(''error when brackting the new root estimation'')'
                   stop ''

                end if

                if(abs(xh-xl).le.xacc) return

             end do

             print '(''ridder_method_module'')'
             print '(''get_root_ridder_method'')'
             print '(''exceeded maximum iterations'')'
             stop ''


          else if (is_real_validated(fl,0.0d0,.false.)) then
             x = x1


          else if (is_real_validated(fh,0.0d0,.false.)) then
             x = x2


          else
             print '(''ridder_method_module'')'
             print '(''get_root_ridder_method'')'
             print '(''root must be bracketed'')'
             
             print '(''(x1,f1): '',2F10.4)', x1,fl
             print '(''(x2,fh): '',2F10.4)', x2,fh
             stop ''


          end if

        end function get_root_ridder_method

      end module ridders_method_module
