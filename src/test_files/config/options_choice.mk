#analyse--
#scalasca instrumentation
ifeq ($(strip $(scalasca)),true)
	PREP=skin
endif


#IFORT options
ifeq ($(FC),ifort)
#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_IFORT)
	LDFLAGS+=$(LDFLAGS_PROFILE_IFORT)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_IFORT)
	LDFLAGS+=$(LDFLAGS_ANALYSE_IFORT)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_IFORT)
	LDFLAGS+=$(LDFLAGS_OP_IFORT)

#OPTIMIZATION REPORT
ifeq ($(strip $(report)),true)
	FFLAGS+=-opt-report3
	FFLAGS+=-opt-report-phase hpo_vectorization
	FFLAGS+=-opt-report-file intel_report.txt 
endif


endif
endif
INCLUDE +=$(INCLUDE_IFORT)
LIBS	+=$(LIBS_IFORT)
endif

#GFORTRAN options
ifeq ($(FC),gfortran)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_PROFILE_GFORTRAN)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_ANALYSE_GFORTRAN)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_OP_GFORTRAN)

endif
endif

INCLUDE +=$(INCLUDE_GFORTRAN)
LIBS	+=$(LIBS_GFORTRAN)
endif

#mpiifort options
ifeq ($(FC),mpiifort)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_IFORT)
	LDFLAGS+=$(LDFLAGS_PROFILE_IFORT)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_IFORT)
	LDFLAGS+=$(LDFLAGS_ANALYSE_IFORT)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_IFORT)
	LDFLAGS+=$(LDFLAGS_OP_IFORT)

#OPTIMIZATION REPORT
ifeq ($(strip $(report)),true)
	FFLAGS+=$(FFLAGS_OPT_REPORT_IFORT)
	LDFLAGS+=$(LDFLAGS_OPT_REPORT_IFORT)
endif

endif
endif

INCLUDE +=$(INCLUDE_IFORT)
LIBS	+=$(LIBS_IFORT)
endif


#mpif90 options
ifeq ($(FC),mpif90)

#PROFILE
ifeq ($(strip $(profile)),true)
	FFLAGS +=$(FFLAGS_PROFILE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_PROFILE_GFORTRAN)
else	

#ANALYSE
ifeq ($(strip $(analyse)),true)
	FFLAGS +=$(FFLAGS_ANALYSE_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_ANALYSE_GFORTRAN)
else

#OPTIMIZED
	FFLAGS +=$(FFLAGS_OP_GFORTRAN)
	LDFLAGS+=$(LDFLAGS_OP_GFORTRAN)

endif
endif

INCLUDE +=$(INCLUDE_GFORTRAN)
INCLUDE +=$(INCLUDE_GCC_PAR)
LIBS	+=$(LIBS_GFORTRAN)
LIBS	+=$(LIBS_GCC_PAR)
endif