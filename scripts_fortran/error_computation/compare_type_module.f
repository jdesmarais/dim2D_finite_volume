      module compare_type_module

        use parameters_kind, only :
     $     ikind,
     $     rkind

        implicit none

        private
        public :: 
     $       compare_rkind,
     $       compare_reals,
     $       compare_doubles,
     $       compare_strings


        contains


        !compare real(rkind)
        function compare_rkind(var,cst,detailled) result(same)

          implicit none

          real(rkind)          , intent(in) :: var
          real(rkind)          , intent(in) :: cst
          logical    , optional, intent(in) :: detailled
          logical                           :: same

          integer, parameter :: precision = 1e9
          logical            :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if
          

          if(detailled_op) then
             print *, nint(var*precision)
             print *, nint(cst*precision)
          end if
          
          same=abs(
     $         nint(var*precision)-
     $         nint(cst*precision)).le.1
          
        end function compare_rkind


        !compare reals
        function compare_reals(var,cst,detailled) result(same)

          implicit none

          real(rkind)          , intent(in) :: var
          real(rkind)          , intent(in) :: cst
          logical    , optional, intent(in) :: detailled
          logical                           :: same

          integer, parameter :: precision = 1e6
          logical            :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if
          

          if(detailled_op) then
             print *, nint(var*precision)
             print *, nint(cst*precision)
          end if
          
          same=abs(
     $         nint(var*precision)-
     $         nint(cst*precision)).le.1
          
        end function compare_reals


        !compare doubles
        function compare_doubles(var,cst,detailled) result(same)

          implicit none

          real*8           , intent(in) :: var
          real*8           , intent(in) :: cst
          logical, optional, intent(in) :: detailled
          logical                       :: same

          integer, parameter :: precision = 1e6
          logical            :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if
          

          if(detailled_op) then
             print *, nint(var*precision)
             print *, nint(cst*precision)
          end if
          
          same=abs(
     $         nint(var*precision)-
     $         nint(cst*precision)).le.1
          
        end function compare_doubles


        !compare strings
        function compare_strings(var,cst,detailled) result(same)

          implicit none

          character*(*)              , intent(in) :: var
          character*(*)              , intent(in) :: cst
          logical      , optional    , intent(in) :: detailled
          logical                                 :: same

          logical :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if
          

          if(detailled_op) then
             print *, var
             print *, cst
          end if
          
          same = len(var).eq.len(cst)

          if(same) then
             same = var.eq.cst
          end if
          
        end function compare_strings

      end module compare_type_module
