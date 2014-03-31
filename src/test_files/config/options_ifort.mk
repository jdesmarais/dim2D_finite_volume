#FFLAGS_OP_IFORT		 = -O2 -p -g -ipo-c -132 -inline-factor=300 -warn all
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=ipo_inl -inline-factor=500
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=ipo_align -inline-factor=500 #check if the data are aligned
FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -inline-factor=500 -warn noalignment #check if there are data that are unaligned
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -vec-report3 -inline-factor=500
#FFLAGS_OP_IFORT		 = -O2 -132 -ipo-c #-D	ALIGNED #-inline-factor=400
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -vec-report3 -guide -inline-factor=400
#FFLAGS_OP_IFORT		 = -O2 -ipo-c -132 -opt-report3 -opt-report-phase=hlo_unroll -inline-factor=400

FFLAGS_PROFILE_IFORT	 = -O0 -g -132 -warn all -check all -r8 -traceback	#profiling options 
FFLAGS_OPT_REPORT_IFORT  = -opt-report3 -opt-report-phase hpo_vectorization -opt-report-file report_intel.txt

#LDFLAGS_OP_IFORT     	 = -O2 -p -g -ipo -132 -inline-factor=300 -warn all
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=ipo_inl -inline-factor=500 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=ipo_align -inline-factor=500 #check if the data are aligned
LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -inline-factor=500 -warn noalignment #check if there are data that are unaligned
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -vec-report3 -inline-factor=300 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -132 -ipo #-DALIGNED #-inline-factor=400 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -vec-report3 -guide -inline-factor=400 #-warn all -check all #GOOD options
#LDFLAGS_OP_IFORT     	 = -O2 -ipo -132 -opt-report3 -opt-report-phase=hlo_unroll -inline-factor=400 #-warn all -check all #GOOD options
LDFLAGS_PROFILE_IFORT	 = -O0 -g -i-dynamic -132 -r8			#profiling options
LDFLAGS_OPT_REPORT_IFORT = -opt-report3 -opt-report-phase hpo_vectorization -opt-report-file report_intel.txt
