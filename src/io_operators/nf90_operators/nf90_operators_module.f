      !> @file
      !> module containing subroutines to manipulate
      !> netcdf files
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module containing subroutines to manipulate
      !> netcdf files
      !
      !> @date
      !> 14_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module nf90_operators_module

        use netcdf

        use parameters_bf_layer, only :
     $       dct_icr_distance,
     $       search_dcr

        use parameters_constant, only :
     $       institut,
     $       prog_version,
     $       commit,
     $       ref,
     $       convention,
     $       
     $       flow_profile_code,
     $       
     $       ns2d_ic_code,
     $       dim2d_ic_code,
     $       
     $       phase_at_center_code,
     $       
     $       bubble_next_to_wall,
     $       bubble_collapse,
     $       bubble_nucleation,
     $       bubble_spherical_cap,
     $       
     $       bc_code,
     $       obc_type_code,
     $       hedstrom_xy_choice,
     $       poinsot_xy_choice,
     $       yoolodato_xy_choice,
     $       
     $       obc_edge_flux_capillarity,
     $       obc_edge_flux_no_capillarity,
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin,
     $       obc_edge_xy_corner,
     $       obc_edge_xy_flux,
     $       obc_edge_xy_diag_flux,
     $       
     $       surface_type_code,
     $       uniform_surface,
     $       surface_with_heaters,
     $       
     $       wall_source_code,
     $       no_heat_source,
     $       gaussian_heat_source
        
        use parameters_kind, only :
     $       rkind,
     $       ikind

        use parameters_input, only :
     $       
     $       npx,npy,nx,ny,ne,bc_size,
     $       x_min,x_max,y_min,y_max,
     $       t_max,dt,detail_print,
     $       
     $       flow_direction,
     $       flow_velocity,
     $       flow_profile,
     $       T0, phase_at_center,
     $       ic_choice, bc_choice,
     $       
     $       bf_openbc_md_threshold_ac,
     $       bf_openbc_md_threshold,
     $       sigma_P,
     $       obc_type_N, obc_type_S,
     $       obc_type_E, obc_type_W,
     $       obc_eigenqties_strategy,
     $       obc_edge_xy_strategy,
     $       obc_edge_flux_strategy,
     $       obc_edge_overlap_ac,
     $       obc_perturbation_T0_ac,
     $       obc_perturbation_vx0_ac,
     $       obc_perturbation_vy0_ac,
     $       obc_perturbation_T0_amp,
     $       obc_perturbation_vx0_amp,
     $       obc_perturbation_vy0_amp,
     $       io_onefile_per_proc,
     $       debug_adapt_computational_domain,
     $       ic_perturbation_ac,
     $       ic_perturbation_amp,
     $       
     $       ratio_bubble_interface,
     $       
     $       wall_surface_type,
     $       wall_micro_contact_angle,
     $       
     $       wall_heater_center,
     $       wall_heater_length,
     $       wall_heater_variation_angle_length,
     $       wall_heater_micro_contact_angle,
     $       
     $       wall_heat_source_choice,
     $       wall_maximum_heat_flux,
     $       wall_heat_source_center,
     $       wall_heat_source_variance,
     $       
     $       wall_extra_heat_source_choice,
     $       wall_maximum_extra_heat_flux,
     $       wall_extra_heat_source_center,
     $       wall_extra_heat_source_variance

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none

        private
        public :: nf90_handle_err,
     $            nf90_open_file,
     $            nf90_write_header,
     $            nf90_def_var_model,
     $            nf90_put_var_model,
     $            nf90_close_file

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> handle the error encountered while manipulating
        !> netcdf files
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !--------------------------------------------------------------
        subroutine nf90_handle_err(nf90_err)

          implicit none

          integer, intent(in) ::  nf90_err
        
          if (nf90_err.ne.NF90_NOERR) then
             print *, 'Error: ', NF90_STRERROR(nf90_err)
             stop 2
          end if

        end subroutine nf90_handle_err


      
        !open a netcdf file for writing
        subroutine nf90_open_file(filename,ncid)

          implicit none

          character*(*), intent(in)  :: filename
          integer      , intent(out) :: ncid

          integer :: retval

          retval = NF90_CREATE(trim(filename), NF90_NETCDF4, ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_open_file


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function giving the name of the file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying the netcdf file
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> reduced simulation time
        !--------------------------------------------------------------
        subroutine nf90_write_header(ncid, p_model)

          implicit none

          integer        , intent(in) :: ncid
          type(pmodel_eq), intent(in) :: p_model


          character(len=10)     :: title
          character(len=29)     :: history
          integer, dimension(8) :: now_data
          integer               :: retval


          !<retrieve the model name
          title = p_model%get_model_name()

          !<retrieve the current date and time from the system
          call date_and_time(values=now_data)
          write (history,
     $         '(''date'', 1x, i4.4, ''/'', i2.2, ''/'', i2.2, 1x,
     $           ''time'', 1x, i2.2, '':'', i2.2, '':'', i2.2)')
     $         now_data(1), now_data(2), now_data(3),
     $         now_data(5), now_data(6), now_data(7)


          !<write the header
          retval = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = nf90_put_att(ncid,NF90_GLOBAL,'history',history)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'institution',institut)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'source',prog_version)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'prog_commit',commit)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'references',ref)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'convention',convention)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<write the characteristic parameters of the simulation
          retval = nf90_put_att(ncid,NF90_GLOBAL,'x_min',x_min)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'x_max',x_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'y_min',y_min)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'y_max',y_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'t_max',t_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'dt',dt)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,NF90_GLOBAL,'detail_print',detail_print)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !initial conditions
          select case(trim(title))
            case('NS2D')
               retval = nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'initial_conditions',
     $              trim(ns2d_ic_code(ic_choice+1)))
               !DEC$ FORCEINLINE RECURSIVE
               call nf90_handle_err(retval)

            case('DIM2D')
               retval = nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'initial_conditions',
     $              trim(dim2d_ic_code(ic_choice+1)))
               !DEC$ FORCEINLINE RECURSIVE
               call nf90_handle_err(retval)

               call nf90_write_header_ic_dim2d(ncid)

               if(ic_perturbation_ac) then
                  call nf90_handle_err(
     $                 nf90_put_att(
     $                 ncid,
     $                 NF90_GLOBAL,
     $                 'ic_perturbation_ac',
     $                 'activated'))

                  call nf90_handle_err(
     $                 nf90_put_att(
     $                 ncid,
     $                 NF90_GLOBAL,
     $                 'ic_perturbation_amp',
     $                 ic_perturbation_amp))

               end if

          end select


          !boundary conditions
          call nf90_write_header_bc(ncid)          


          !<write the characteristic parameters for the physical model
          call nf90_write_sim_parameters(ncid,p_model)


        end subroutine nf90_write_header


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function writing the i.c. parameters on the
        !> output netcdf file
        !
        !> @date
        !> 05_06_2015 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying the netcdf file
        !--------------------------------------------------------------
        subroutine nf90_write_header_ic_dim2d(ncid)

          implicit none

          integer, intent(in) :: ncid

          integer :: retval


          ! parameters specific to the i.c.
          ! are written

          !1) flow velocity
          !----------------------------------------
          retval = nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'flow_velocity',
     $         flow_velocity)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !2) temperature
          !----------------------------------------
          retval = nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'temperature',
     $         T0)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !3) parameters characteristic to the
          !   initial conditions
          !----------------------------------------
          select case(ic_choice)

          case(bubble_next_to_wall,bubble_collapse,
     $         bubble_nucleation, bubble_spherical_cap)

             if(abs(flow_velocity)>0.0) then
             
               retval = nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'flow_profile',
     $            flow_profile_code(flow_profile+1))
               !DEC$ FORCEINLINE RECURSIVE
               call nf90_handle_err(retval)

             end if

             retval = nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'phase_at_center',
     $         phase_at_center_code(phase_at_center+1))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             if(ic_choice.eq.bubble_collapse) then

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'ratio_bubble_radius_interface',
     $               ratio_bubble_interface)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

             end if

             ! type of surface
             retval = nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'wall_surface_type',
     $            surface_type_code(wall_surface_type+1))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)


             ! type of surface
             retval = nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'wall_micro_contact_angle',
     $            wall_micro_contact_angle)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             
             ! add the parameters specific to the
             ! non-uniform surfaces
             if(wall_surface_type.eq.surface_with_heaters) then
                
                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_heater_center',
     $               wall_heater_center)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_heater_length',
     $               wall_heater_length)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)
                
                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_heater_variation_angle_length',
     $               wall_heater_variation_angle_length)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_heater_micro_contact_angle',
     $               wall_heater_micro_contact_angle)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

             end if


             ! conduction heat source
             if(wall_heat_source_choice.ne.no_heat_source) then
                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_heat_source',
     $               wall_source_code(wall_heat_source_choice+1))
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_maximum_heat_source',
     $               wall_maximum_heat_flux)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                if(wall_heat_source_choice.eq.gaussian_heat_source) then

                   if(wall_surface_type.ne.surface_with_heaters) then

                      retval = nf90_put_att(
     $                     ncid,
     $                     NF90_GLOBAL,
     $                     'wall_heat_source_center',
     $                     wall_heat_source_center)
                      !DEC$ FORCEINLINE RECURSIVE
                      call nf90_handle_err(retval)

                      retval = nf90_put_att(
     $                     ncid,
     $                     NF90_GLOBAL,
     $                     'wall_heat_source_variance',
     $                     wall_heat_source_variance)
                      !DEC$ FORCEINLINE RECURSIVE
                      call nf90_handle_err(retval)

                   end if

                end if

             end if

             
             ! extra heat source
             if(wall_extra_heat_source_choice.ne.no_heat_source) then
                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_extra_heat_source',
     $               wall_source_code(wall_extra_heat_source_choice+1))
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               'wall_maximum_extra_heat_flux',
     $               wall_maximum_extra_heat_flux)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

                if(wall_extra_heat_source_choice.eq.gaussian_heat_source) then

                   if(wall_surface_type.ne.surface_with_heaters) then

                      retval = nf90_put_att(
     $                     ncid,
     $                     NF90_GLOBAL,
     $                     'wall_extra_heat_source_center',
     $                     wall_extra_heat_source_center)
                      !DEC$ FORCEINLINE RECURSIVE
                      call nf90_handle_err(retval)

                      retval = nf90_put_att(
     $                     ncid,
     $                     NF90_GLOBAL,
     $                     'wall_extra_heat_source_variance',
     $                     wall_extra_heat_source_variance)
                      !DEC$ FORCEINLINE RECURSIVE
                      call nf90_handle_err(retval)

                   end if

                end if

             end if

          end select

        end subroutine nf90_write_header_ic_dim2d


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function writing the b.c. parameters on the
        !> output netcdf file
        !
        !> @date
        !> 16_01_2015 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying the netcdf file
        !--------------------------------------------------------------
        subroutine nf90_write_header_bc(ncid)

          implicit none

          integer, intent(in) :: ncid

          integer :: retval

          ! bc_type
          retval = nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'boundary_conditions',
     $         trim(bc_code(bc_choice+1)))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          ! parameters specific to the b.c.
          ! are written
          
          ! open b.c.
          select case(bc_choice)

          case(hedstrom_xy_choice,
     $         poinsot_xy_choice,
     $         yoolodato_xy_choice)


          ! how the eigenquantities are computed at the
          ! borders
          select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_eigenqties_strategy',
     $              'obc_eigenqties_bc'))

            case(obc_eigenqties_lin)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_eigenqties_strategy',
     $              'obc_eigenqties_lin'))

          end select


          ! type of obc_edge_xy_strategy: strategy for anti-corners
          select case(obc_edge_xy_strategy)

            case(obc_edge_xy_corner)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_edge_xy_strategy',
     $              'obc_edge_xy_corner'))

            case(obc_edge_xy_flux)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_edge_xy_strategy',
     $              'obc_edge_xy_flux'))

            case(obc_edge_xy_diag_flux)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_edge_xy_strategy',
     $              'obc_edge_xy_diag_flux'))

          end select


          ! how the fluxes are computed in the
          ! b.c.
          select case(obc_edge_flux_strategy)

            case(obc_edge_flux_capillarity)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_edge_flux_strategy',
     $              'obc_edge_flux_capillarity'))

            case(obc_edge_flux_no_capillarity)

               call nf90_handle_err(
     $              nf90_put_att(
     $              ncid,
     $              NF90_GLOBAL,
     $              'obc_edge_flux_strategy',
     $              'obc_edge_flux_no_capillarity'))

          end select


          ! whether edge overlap anti-corners
          ! for the treatment of the bc_sections
          if(obc_edge_overlap_ac) then
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_edge_overlap_ac',
     $            'activated'))
          else
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_edge_overlap_ac',
     $            'deactivated'))
          end if


          ! activation of the perturbations
          ! in the far field values
          if(obc_perturbation_T0_ac) then
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_T0_ac',
     $            'activated'))

             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_T0_amp',
     $            obc_perturbation_T0_amp))
          end if

          if(obc_perturbation_vx0_ac) then
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_vx0_ac',
     $            'activated'))

             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_vx0_amp',
     $            obc_perturbation_vx0_amp))
          end if

          if(obc_perturbation_vy0_ac) then
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_vy0_ac',
     $            'activated'))

             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'obc_perturbation_vy0_amp',
     $            obc_perturbation_vy0_amp))
          end if


          ! threshold for detection of interface
          if(bf_openbc_md_threshold_ac) then
             retval = nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'openbc_md_threshold_ac',
     $            'activated')
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             retval = nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'openbc_md_threshold',
     $            bf_openbc_md_threshold)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end if


          end select


          ! poinsot and yoo-lodato open b.c.
          select case(bc_choice)

          case(poinsot_xy_choice,
     $         yoolodato_xy_choice)

          retval = nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'openbc_sigma_P',
     $         sigma_P)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          
          retval= nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'openbc_type_N',
     $         trim(obc_type_code(obc_type_N+1)))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval= nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'openbc_type_S',
     $         trim(obc_type_code(obc_type_S+1)))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval= nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'openbc_type_E',
     $         trim(obc_type_code(obc_type_E+1)))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval= nf90_put_att(
     $         ncid,
     $         NF90_GLOBAL,
     $         'openbc_type_W',
     $         trim(obc_type_code(obc_type_W+1)))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          end select


          ! for the adaptation of the computational domain
          if(debug_adapt_computational_domain) then
             
             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'adapt_computational_domain',
     $            'activated'))

             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'dct_icr_distance',
     $            dct_icr_distance))

             call nf90_handle_err(
     $            nf90_put_att(
     $            ncid,
     $            NF90_GLOBAL,
     $            'dct_dcr_distance',
     $            search_dcr))

          end if             

        end subroutine nf90_write_header_bc


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function writing the simulation parameters on the
        !> output netcdf file
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying the netcdf file
        !
        !>@param p_model
        !> physical model
        !--------------------------------------------------------------
        subroutine nf90_write_sim_parameters(ncid,p_model)

          implicit none

          integer        , intent(in) :: ncid
          type(pmodel_eq), intent(in) :: p_model


          character(20), dimension(:), allocatable :: param_name
          real(rkind)  , dimension(:), allocatable :: param_value
          integer                                  :: i, retval

          
          !get the simulation parameters
          call p_model%get_sim_parameters(param_name,param_value)

          if(allocated(param_name)) then

            !save the simulation parameters in the netcdf file
             do i=1, size(param_name)

                retval = nf90_put_att(
     $               ncid,
     $               NF90_GLOBAL,
     $               param_name(i),
     $               param_value(i))
                
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

             end do


             !deallocate unused data
             deallocate(param_name)
             deallocate(param_value)

          end if

        end subroutine nf90_write_sim_parameters


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine defining the variables saved in
        !> the netcdf file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param p_model
        !> physical model
        !
        !>@param coordinate_id
        !> integer table identifying the location of the coordinate
        !> tables in the netcdf file
        !
        !>@param data_id
        !> integer table identifying the location of the main varibles
        !> of the governing equations in the netcdf file
        !--------------------------------------------------------------
        subroutine nf90_def_var_model(
     $     ncid,
     $     p_model,
     $     coordinates_id,
     $     data_id)

          implicit none

          integer               , intent(in)    :: ncid
          type(pmodel_eq)       , intent(in)    :: p_model
          integer, dimension(3) , intent(inout) :: coordinates_id
          integer, dimension(ne), intent(inout) :: data_id


          character*(*), parameter :: T_NAME = 'time'
          character*(*), parameter :: X_NAME = 'x'
          character*(*), parameter :: Y_NAME = 'y'

          character*(*), parameter :: UNITS    = 'units'
          character*(*), parameter :: T_UNITS  = 's/s'
          character*(*), parameter :: X_UNITS  = 'm/m'
          character*(*), parameter :: Y_UNITS  = 'm/m'

          character*(*), parameter :: AXIS     = 'axis'
          character*(*), parameter :: X_AXIS   = 'X'
          character*(*), parameter :: Y_AXIS   = 'Y'
          
          character*(*), parameter :: LONG_NAME  ='long_name'
          character*(*), parameter :: T_LONG_NAME='reduced time'
          character*(*), parameter :: X_LONG_NAME='reduced x-coordinate'
          character*(*), parameter :: Y_LONG_NAME='reduced y-coordinate'

          integer, parameter       :: NDIMS = 3
          integer, dimension(NDIMS):: dimids

          integer     :: NF_MYREAL
          real(rkind) :: NF_FILL_MYREAL

          character(len=10), dimension(ne) :: name_var
          character(len=33), dimension(ne) :: longname_var
          character(len=23), dimension(ne) :: unit_var

          integer :: t_dimid
          integer :: x_dimid
          integer :: y_dimid
      
          integer :: t_varid
          integer :: x_varid
          integer :: y_varid

          integer(ikind) :: nx_domain
          integer(ikind) :: ny_domain

          integer :: k
          integer :: retval


          !<define the type of variables stored
          select case(RKIND)

            case(4)
               NF_MYREAL      = NF90_FLOAT
               NF_FILL_MYREAL = NF90_FILL_FLOAT

            case(8)
               NF_MYREAL      = NF90_DOUBLE
               NF_FILL_MYREAL = NF90_FILL_DOUBLE
               
            case default
               print '(''nf90_operators_wr_class :'')'
               print '(''nf90_def_var_model'')'
               stop 'NF_MYREAL'
          end select


          !<define the dimensions
          if(io_onefile_per_proc) then
             nx_domain = nx
             ny_domain = ny
          else
             nx_domain = npx*(nx-2*bc_size) + 2*bc_size
             ny_domain = npy*(ny-2*bc_size) + 2*bc_size
          end if

          retval = NF90_DEF_DIM(ncid, T_NAME, 1, t_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, X_NAME, nx_domain, x_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, Y_NAME, ny_domain, y_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<define the coordinates variables
          retval = NF90_DEF_VAR(ncid, T_NAME, NF_MYREAL, t_dimid,t_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_VAR(ncid, X_NAME, NF_MYREAL, x_dimid,x_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_DEF_VAR(ncid, Y_NAME, NF_MYREAL, y_dimid,y_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the units to the variables
          retval = NF90_PUT_ATT(ncid, t_varid, UNITS, T_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, UNITS, X_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, UNITS, Y_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the type of axis to the coordinates
          retval = NF90_PUT_ATT(ncid, x_varid, AXIS, X_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, AXIS, Y_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the long name for the description of the coordinates
          retval = NF90_PUT_ATT(ncid, t_varid, LONG_NAME, T_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, LONG_NAME, X_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, LONG_NAME, Y_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<assign the coordinates id
          coordinates_id(1) = t_varid
          coordinates_id(2) = x_varid
          coordinates_id(3) = y_varid


          !<define the extension in x and y directions
          dimids(1) = t_dimid
          dimids(2) = x_dimid
          dimids(3) = y_dimid

          
          !<define the main variables of the governing equations
          !DEC$ FORCEINLINE RECURSIVE
          name_var     = p_model%get_var_name()
          !DEC$ FORCEINLINE RECURSIVE
          longname_var = p_model%get_var_longname()
          !DEC$ FORCEINLINE RECURSIVE
          unit_var     = p_model%get_var_unit()

          do k=1, ne

             !<define the netcdf variables
             retval = NF90_DEF_VAR(ncid, trim(name_var(k)), NF_MYREAL,
     $            dimids, data_id(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign the units to the variables
             retval = NF90_PUT_ATT(ncid, data_id(k), UNITS, unit_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign the long_name to the variables
             retval = NF90_PUT_ATT(ncid, data_id(k), LONG_NAME,
     $            longname_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign missing value to
             !the governing variables
             retval = nf90_put_att(
     $            ncid,
     $            data_id(k),
     $            'missing_data',
     $            NF_FILL_MYREAL)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign fill value to
             !the governing variables
             retval = nf90_put_att(
     $            ncid,
     $            data_id(k),
     $            '_FillValue',
     $            NF_FILL_MYREAL)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end do
          
          !<stop the definition of the variables saved in the file
          retval = NF90_ENDDEF(ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_def_var_model


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> put the main variables in the netcdf file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !> @param ncid
        !> integer identifying the netcdf file opened
        !
        !> @param coordinates_id
        !> table with the id identifying the coordinates
        !
        !> @param data_id
        !> table with the id identifying the data
        !
        !> @param time
        !> reduced simulation time
        !
        !> @param field_used
        !> object encapsulating the main variables
        !
        !> @param start
        !> integer, dimension(2) identifying the sw point
        !> where the partial filling starts
        !
        !> @param count
        !> integer, dimension(2) identifying the profile
        !> of the gridpoints written
        !
        !> @param limit
        !> integer, dimension(4) identifying the borders
        !> of the sub-tile written out of the computational
        !> domain
        !----------------------------------------------------------------
        subroutine nf90_put_var_model(
     $     ncid,
     $     coordinates_id,
     $     data_id,
     $     time, nodes, x_map, y_map,
     $     start, count, limit)

          implicit none
          
          integer                         , intent(in) :: ncid
          integer    , dimension(3)       , intent(in) :: coordinates_id
          integer    , dimension(ne)      , intent(in) :: data_id
          real(rkind)                     , intent(in) :: time
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          integer, dimension(2), optional , intent(in) :: start
          integer, dimension(2), optional , intent(in) :: count
          integer, dimension(4), optional , intent(in) :: limit


          integer    , dimension(3) :: start_op, count_op
          real(rkind), dimension(1) :: time_table
          integer :: t_varid
          integer :: x_varid
          integer :: y_varid
          integer :: k
          integer :: retval


          time_table(1) = time

          t_varid = coordinates_id(1)
          x_varid = coordinates_id(2)
          y_varid = coordinates_id(3)

          
          !< write the time
          retval = NF90_PUT_VAR(ncid, t_varid, time_table)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !< depending if the start and count optional variables
          !> are present, we should write the whole domain or
          !> only a part identified by (start, count) with the
          !> exception of the time
          !-----------------------------------------------------
          if(.not.present(start)) then

             start_op = [1,1,1]
             count_op = [1,nx,ny]

             !< write the x_map coordinates
             retval = NF90_PUT_VAR(ncid, x_varid, x_map)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)


             !< write the y_map coordinates
             retval = NF90_PUT_VAR(ncid, y_varid, y_map)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)
             

             !< write the variables of the governing equations
             do k=1, ne

                retval = NF90_PUT_VAR(
     $               ncid,
     $               data_id(k),
     $               nodes(:,:,k),
     $               START=start_op,
     $               COUNT=count_op)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

             end do

         else

            !< write the x_map coordinates
            retval = NF90_PUT_VAR(
     $           ncid,
     $           x_varid,
     $           x_map(limit(1):limit(2)),
     $           START=[start(1)],
     $           COUNT=[count(1)])
            !DEC$ FORCEINLINE RECURSIVE
            call nf90_handle_err(retval)


            !< write the y_map coordinates
            retval = NF90_PUT_VAR(
     $           ncid,
     $           y_varid,
     $           y_map(limit(3):limit(4)),
     $           START=[start(2)],
     $           COUNT=[count(2)])
            !DEC$ FORCEINLINE RECURSIVE
            call nf90_handle_err(retval)

            
            !< write the variables of the governing equations
            do k=1, ne

               retval = NF90_PUT_VAR(
     $              ncid,
     $              data_id(k),
     $              nodes(limit(1):limit(2),limit(3):limit(4),k),
     $              START=[1, start(1), start(2)], 
     $              COUNT=[1, count(1), count(2)])
               !DEC$ FORCEINLINE RECURSIVE
               call nf90_handle_err(retval)

            end do

          end if

        end subroutine nf90_put_var_model


        !close netcdf file
        subroutine nf90_close_file(ncid)
        
          implicit none

          integer, intent(in) :: ncid

          integer :: retval

          retval = NF90_CLOSE(ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_close_file

      end module nf90_operators_module
