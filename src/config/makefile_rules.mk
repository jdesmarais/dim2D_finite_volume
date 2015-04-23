help:
	@(echo 'make all      : create executable for all test files')
	@(echo 'make clean    : remove .o files')
	@(echo 'make cleanall : remove .o and executable files')
	@(echo 'make clear    : remove .err .out and .nc')
	@(echo 'make test     : display makefile options')

	@(echo 'to change the cluster configuration, change the env. var. AUGEANSTABLES_CLUSTER')
	@(echo 'to change the compiler, change the env. var. AUGEANSTABLES_COMPILER')
	@(echo 'to activate profiling, change the env. var. AUGEANSTABLES_PROFILE')
	@(echo 'to activate analyzing, change the env. var. AUGEANSTABLES_ANALYSE')
	@(echo 'to activate tracing, change the env. var. AUGEANSTABLES_TRACE')
	@(echo 'to activate report, change the env. var. AUGEANSTABLES_OPTREPORT')
	@(echo 'all these variables can also be set in the corresponding module file')
	@(echo '')

doc:
	$(config_dir)/generate_doc.sh
	@(echo 'documentation generated')
	@(echo 'opening documentation...')
	@(firefox ./doc/html/files.html &)

dir:
	@($(config_dir)/display_dir.sh 'bc_dir:         ' $(bc_dir)          )
	@($(config_dir)/display_dir.sh 'field_dir:      ' $(field_dir)	     )
	@($(config_dir)/display_dir.sh 'io_dir:         ' $(io_dir)	     )
	@($(config_dir)/display_dir.sh 'mpi_dir:        ' $(mpi_dir)	     )
	@($(config_dir)/display_dir.sh 'param_dir:      ' $(param_dir)       )
	@($(config_dir)/display_dir.sh 'phy_eq_dir:     ' $(phy_eq_dir)      )
	@($(config_dir)/display_dir.sh 'sd_dir:         ' $(sd_dir)	     )
	@($(config_dir)/display_dir.sh 'ti_dir:         ' $(ti_dir)	     )
	@($(config_dir)/display_dir.sh 'td_dir:         ' $(td_dir)          )
	@($(config_dir)/display_dir.sh 'test_dir:       ' $(test_dir)        )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'config_dir:     ' $(config_dir)	     )
	@(echo '')
	@($(config_dir)/display_dir.sh 'cg_dir:         ' $(cg_dir) 	     )
	@($(config_dir)/display_dir.sh 'mt_dir:         ' $(mt_dir)          )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'rbc_dir:        ' $(rbc_dir)	     )
	@($(config_dir)/display_dir.sh 'pbc_dir:        ' $(pbc_dir)	     )
	@($(config_dir)/display_dir.sh 'wbc_dir:        ' $(wbc_dir)	     )
	@($(config_dir)/display_dir.sh 'wrbc_dir:       ' $(wrbc_dir)        )
	@($(config_dir)/display_dir.sh 'obc_dir:        ' $(obc_dir)         )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'nf90_dir:       ' $(nf90_dir)        )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'mpi_bc_dir:     ' $(mpi_bc_dir)      )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'simpletest_dir: ' $(simpletest_dir)  )
	@($(config_dir)/display_dir.sh 'wave1d_dir:     ' $(wave1d_dir)      )
	@($(config_dir)/display_dir.sh 'wave2d_dir:     ' $(wave2d_dir)      )
	@($(config_dir)/display_dir.sh 'ns2d_dir:       ' $(ns2d_dir)        )
	@($(config_dir)/display_dir.sh 'dim2d_dir:      ' $(dim2d_dir)       )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'rk3tvd_dir:     ' $(rk3tvd_dir)      )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'fv_dir:         ' $(fv_dir)	     )
	@(echo '' )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'ns2d_ic:        ' $(ns2d_ic)         )
	@($(config_dir)/display_dir.sh 'dim2d_ic:       ' $(dim2d_ic)        )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'ns2d_sic:       ' $(ns2d_sic)       )
	@($(config_dir)/display_dir.sh 'ns2d_pic:       ' $(ns2d_pic)       )
	@($(config_dir)/display_dir.sh 'ns2d_vic:       ' $(ns2d_vic)       )
	@($(config_dir)/display_dir.sh 'ns2d_sxic:      ' $(ns2d_sxic)      )
	@($(config_dir)/display_dir.sh 'ns2d_syic:      ' $(ns2d_syic)      )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'hobc_dir:       ' $(hobc_dir)        )
	@($(config_dir)/display_dir.sh 'hcobc_dir:      ' $(hcobc_dir)       )
	@($(config_dir)/display_dir.sh 'hrobc_dir:      ' $(hrobc_dir)       )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'pobc_dir:       ' $(pobc_dir)        )
	@($(config_dir)/display_dir.sh 'pobc_ns2d_dir:  ' $(pobc_ns2d_dir)   )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'yobc_dir:       ' $(yobc_dir)        )
	@($(config_dir)/display_dir.sh 'yobc_ns2d_dir:  ' $(yobc_ns2d_dir)   )
	@(echo '' )
	@($(config_dir)/display_dir.sh 'bf_layer_dir:   ' $(bf_layer_dir)    )
	@($(config_dir)/display_dir.sh 'nbf_layer_dir:  ' $(nbf_layer_dir)   )
	@($(config_dir)/display_dir.sh 'dbf_layer_dir:  ' $(dbf_layer_dir)   )
	@($(config_dir)/display_dir.sh 'sbf_layer_dir:  ' $(sbf_layer_dir)   )
	@($(config_dir)/display_dir.sh 'iobf_layer_dir: ' $(iobf_layer_dir)  )
	@($(config_dir)/display_dir.sh 'cbf_layer_dir:  ' $(cbf_layer_dir)   )

cdir:
	@(echo 'sd_cdir:   '   	 $(sd_cdir)          )
	@(echo 'pm_cdir:   '	 $(pm_cdir)          )
	@(echo 'bc_cdir:   '	 $(bc_cdir)          )
	@(echo 'ic_cdir:   '     $(ic_cdir)          )
	@(echo 'pobc_cdir: '	 $(pobc_cdir)        )
	@(echo 'yobc_cdir: '	 $(yobc_cdir)        )

test:
	@(echo '')

	@(echo 'source code folder: ' $(src_dir))
	@(echo 'cluster:            ' $(cluster))
	@(echo '')

	@(echo 'profile:   ' $(profile))
	@(echo 'analyse:   ' $(analyse))
	@(echo 'trace:     ' $(trace))
	@(echo 'report:    ' $(report))
	@(echo 'scalasca:  ' $(scalasca))
	@(echo '')

	@(echo 'PREP:    ' $(PREP))
	@(echo 'FC:      ' $(FC))
	@(echo 'FFLAGS:  ' $(FFLAGS))
	@(echo 'LDFLAGS: ' $(LDFLAGS))
	@(echo 'INCLUDE: ' $(INCLUDE))
	@(echo 'LIBS:    ' $(LIBS))
	@(echo '')

	@(echo 'sd_choice: ' $(sd_choice))
	@(echo 'pm_choice: ' $(pm_choice))
	@(echo 'ic_choice: ' $(ic_choice))
	@(echo 'bc_choice: ' $(bc_choice))
	@(echo 'td_choice: ' $(td_choice))
	@(echo 'ti_choice: ' $(ti_choice))
	@(echo 'io_choice: ' $(io_choice))

	@(echo $(test_dir))

tests:	test_field\
	test_cg_operators\
	test_dim2d_prim\
	test_dim2d_fluxes\
	test_dim2d_eq\
	test_fv_operators\
	test_bc_periodic\
	test_bc_reflection\
	test_wall_xy_module\
	test_bc_wall\
	test_rk3tvd\
	test_nf90_operators\
	test_dim2d_ic\
	test_rk3tvd_dim2d

tests_par: test_field_par\
	test_mpi_mg_bc

sims:	sim_dim2d

#objects
#$(eval file := $(shell basename $@))
#@(printf "compilation: %-60s\r" $(file))
%.o:	%.f
	@(basename $@)
	@($(PREP) $(FC) $(FFLAGS) -c $< $(INCLUDE))

#executable
%:	%.o
	@(basename $@)
	@($(PREP) $(FC) $(LDFLAGS) -o $@ $^ $(LIBS))

#main commmand
all:	$(CMD)

#cleaning previous simulations
clear:
	rm *.err *.out *.nc

#cleaning
clean:	
	rm -f *.o *.mod core *~

#complete cleaning
cleanall:
	rm -f $(subst .exe,_mono.exe,$(CMD)) $(CMD) *.o *.mod core *~


