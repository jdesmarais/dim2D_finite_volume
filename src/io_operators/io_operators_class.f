      !> @file
      !> class encapsulating subroutines to write
      !> data on output files (ex:netcdf files)
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to write
      !> data on output files (ex:netcdf files)
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module io_operators_class

        use field_class       , only : field
        use parameters_input  , only : io_choice
        use parameters_kind   , only : rkind
        use phy_model_eq_class, only : phy_model_eq

        implicit none
        
        private
        public : io_operators


        !> @class io_operators
        !> class encapsulating subroutines to write
        !> data on output files (ex:netcdf files)
        !>
        !> @param write_data
        !> write data on output files
        !---------------------------------------------------------------
        type :: io_operators

          contains
          
          procedure, nopass, non_overridable :: write_data

        end type io_operators


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine choosing the output format according to
        !> the user choice and writing the output data
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param p_model
        !> physical model for the name, description and units
        !> of variables
        !
        !>@param t
        !> reduced simulation time
        !--------------------------------------------------------------
        subroutine write_data(field_used, p_model,bc_size,t)

          implicit none

          class(field)       , intent(in) :: field_used
          class(phy_model_eq), intent(in) :: p_model
          integer            , intent(in) :: bc_size
          real(rkind)        , intent(in) :: t

          type(nf90_operators_wr) :: nf90_writer


          !<select the i/o strategy
          select case(io_choice)
            case(netcdf_choice)
               call nf90_writer%write_data(field_used,p_model,bc_size)

            case default
               print '(''io_operators: write_data'')'
               stop 'io_choice not recognize'
               
          end select

        end subroutine write_data

      end module io_operators_class
