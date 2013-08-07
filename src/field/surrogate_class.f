      !this stateless type serves only for purposese of extension by
      !other types: it can serve as the child type when that type is
      !inaccessible b/c of Fortran's prohibition against circular
      !references
      !------------------------------------------------------------------
      module surrogate_class

        implicit none
        
        private
        public :: surrogate

        type, abstract ::surrogate

        end type surrogate

      end module surrogate_class
