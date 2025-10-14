# -*- mode: makefile -*-
# Makefile for PHITS, last revised 2021/04/13
# If you want to use many cores for compiling PHITS, you have to type
# make -j xx,
# where xx is number of cores

### Machine Dependent variables, please set your environment
ENVFLAGS = LinIfort
# Linux Intel Fortran (ifort):                             LinIfort
# Linux Intel Fortran (ifx):                               LinIFX
# Linux gfortran:                                          LinGfort
# Linux NVIDIA compiler (former PGI Fortran):              LinNVIDIA
# Windows gfortran (4.8 or 7 series):                      WinGfort
# Windows Cygwin gfortran (4.8 or 7 series):               WinCygGfort
# Mac Intel fortran:                                       MacIfort
# Mac PGI fortran:                                         MacPGIfort
# Mac gfortran (4.71 or later):                            MacGfort
# SGI ICE X (JAEA supercomputer since Nov. 2015):          SGIICEX
# HPE SGI8600 (JAEA supercomputer since Nov. 2020):        SGI8600
# FX1 (you have to use "gmake" instead of "make"):         FX1
### If you succeeded in compiling PHITS in new environment,
### please send your makefile to phits-office@jaea.go.jp

### If you want to use MPI, delete # in the next line
# USEMPI = true
### If you want to use OpenMP, delete # in the next line.
# USEOMP = true

### If you set USEOMP = true,
### You have to set an environmental variable OMP_NUM_THREADS before executing PHITS.
### If you do not know the appropriate value of OMP_NUM_THREADS for your machine,
### please ask your machine administrator, or turn off USEOMP.

### If you want to make dependence (make.depend), you have to type
### make depend DEVEL=1

DEBUG ?= false

### This is for future option. Currently, PHITS cannot use GPU.
# USEGPU = true

### If you want to use HDF5, set true in the following line and
### set correct HDF5 installed path where lib/ and inlude/ directories exist
### to HDF5PATH.
USEHDF5 ?= false
ifeq ($(USEHDF5),true)
 HDF5PATH = ../utility/HDF5
 ifeq ("$(wildcard $(HDF5PATH)/include/*.h)","")
  $(error Error! header files are not found in $(HDF5PATH)/include)
 else
  INCLUDES = -I$(HDF5PATH)/include
  ifeq ("$(wildcard $(HDF5PATH)/include/*.mod)","")
   ifeq ("$(wildcard $(HDF5PATH)/mod/*.mod)","")
    ifeq ("$(wildcard $(HDF5PATH)/mod/static/*.mod)","")
     $(error Error! module file are not found in $(HDF5PATH)/include)
    else
     INCLUDES += -I$(HDF5PATH)/mod/static
    endif
   else
    INCLUDES += -I$(HDF5PATH)/mod/
   endif
  endif
 endif
 ifeq ("$(wildcard $(HDF5PATH)/lib/*.a)","")
  ifeq ("$(wildcard $(HDF5PATH)/lib/*.so)","")
   $(error Error! library files are not found in $(HDF5PATH)/lib)
  endif
 endif
 OPTFLAGS = -Dhdf5
 LDLIBS += -L$(HDF5PATH)/lib
 LDLIBS += -Bstatic -lhdf5_fortran -lhdf5 -lz -Wl,-rpath,$(HDF5PATH)/lib
 # LDLIBS += -Bstatic -lhdf5_fortran -lhdf5 -lz -lsz -laec  # for Linux IFX
endif

##### target name
ifeq ($(ENVFLAGS),)
 $(error Error! ENVFLAGS is not defined)
else
 TARGET   = ../phits_$(ENVFLAGS)
 ifeq ($(USEOMP),true)
  TARGET := $(TARGET)_OMP
 endif
 ifeq ($(USEMPI),true)
  # $(warning You may need to set the fortran compiler that is invoked from your mpifort command)
  TARGET := $(TARGET)_MPI
 endif
 ifeq ($(DEBUG),true)
  TARGET := $(TARGET)_DBG
 endif
 ifeq ($(USEGPU),true)
  TARGET := $(TARGET)_GPU
 endif
endif

ifeq ($(ENVFLAGS),WinGfort)
 RM = del
else
 RM = rm
endif


### SRCS8
ifeq ($(USEMPI),true)
 ifeq ($(ENVFLAGS),WinGfort)
  SRCS8 = mpi-win.f unix90.f
 else
  SRCS8 = mpi-lin.f unix90.f
 endif
else
 SRCS8 = mpi-non.f unix90.f
endif

### Linux Intel Fortran (ifort)
ifeq ($(ENVFLAGS),LinIfort)
 IntelFortran = true
 SRCS8  += mdp-uni90.f
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -fp-stack-check -check uninit -check bounds
 else
  FFLAGS = -O3 -fpconstant
  ### FFLAGS  = -O3 -fpconstant -mcmodel=large -shared-intel ### in the case of large mdas
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -qopenmp  ### if -qopenmp does not work, try -openmp
  LD = $(FC) -static-intel
 endif
 ifeq ($(USEMPI),true)
  FC = mpiifort
  LDLIBS += -lmpi ### -lmpi sometimes -lmpich, or not necessary
 else
  FC = ifort
 endif
endif

### Linux Intel Fortran (ifx)
ifeq ($(ENVFLAGS),LinIFX)
  IntelFortran = true
  LDLIBS += -static-intel
  SRCS8  += mdp-uni90.f
  ifeq ($(DEBUG),true)
#   FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -check uninit,bounds -mcmodel=medium -Wl,--no-relax ### in the case of checking uninit. Note that many errors occur.
   FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -check bounds
  else
   FFLAGS = -O3 -fpconstant
  endif
  ifeq ($(USEOMP),true)
   FFLAGS += -fiopenmp  ### if -fiopenmp does not work, try -qopenmp
   LD = $(FC) -static-intel
   LDLIBS += -qopenmp-link=static
  endif
  ifeq ($(USEMPI),true)
   FC = mpiifx
   LDLIBS += -lmpi ### -lmpi sometimes -lmpich, or not necessary
  else
   FC = ifx
  endif
endif

### Linux PGI Fortran
ifeq ($(ENVFLAGS),LinNVIDIA)
 FC = pgfortran
 SRCS8  += mdp-uni90.f
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -g
 else
  FFLAGS = -O2 -r8 -fast -mcmodel=medium -Mfixed -Mpreprocess -Mbackslash
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -mp
 endif
 ifeq ($(USEMPI),true)
  FC = mpif90
  # FFLAGS += -Mmpi=mpich # for compiling w/ pgfortran
 endif
 ifeq ($(USEGPU),true)
    FFLAGS += -acc -ta=tesla:cc70,rdc -Mcuda,cc70 -Mcuda=cuda11.0,ptxinfo
 endif

endif

### Linux, gfortran (4.8, 7.0 or later)
ifeq ($(ENVFLAGS),LinGfort)
 OtherGfortran = true
 FC        = gfortran
 SRCS8    += mdp-uni90.f
 FFLAGS    = -fdefault-double-8 -fdefault-real-8 -fdollar-ok -std=legacy
 ifeq ($(DEBUG),true)
  FFLAGS += -O0 -g
 endif
 ifeq ($(USEOMP),true)
  FFLAGS  += -fopenmp
 endif
 ifeq ($(USEMPI),true)
  FC       = mpif90
  #LDLIBS  += -lmpi_usempi -lmpi_mpifh -lmpi
 endif
endif

### SGI ICE X (JAEA supercomputer since Nov. 2015)
ifeq ($(ENVFLAGS),SGIICEX)
 IntelFortran = true
 SRCS8  += mdp-uni90.f
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -fp-stack-check
 else
  FFLAGS = -O3 -fpconstant
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -qopenmp
 endif
 ifeq ($(USEMPI),true)
  FC = mpif90
  FFLAGS += -lmpi
 else
  FC = ifort
 endif
endif

### SGI 8600 (JAEA supercomputer since Dec. 2020)
ifeq ($(ENVFLAGS),SGI8600)
   IntelFortran = true
   SRCS8  += mdp-uni90.f
   ifeq ($(DEBUG),true)
       FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -fp-stack-check
   else
       FFLAGS = -O3 -fpconstant
   endif
   ifeq ($(USEOMP),true)
       FFLAGS += -qopenmp
   endif
   ifeq ($(USEMPI),true)
       FC = mpif90
       FFLAGS += -lmpi
   else
       FC = ifort
   endif
endif

### FX1, OpenMP uncheck, you have to use "gmake" instead of "make"
ifeq ($(ENVFLAGS),FX1)
 FC       = frt
 SRCS8   += mdp-uni.f
 FFLAGS   = -fw -X9 -Am $(OMPFLAGS) \
            -xzdistfunc,Gau -xtime_b -xshellspar -xfnf \
            -xsmmmod.xi,smmmod.eps
 FFLAGS  += -Kfast,ocl,preex
 # OPTFLAGS = -Kfast,ocl,preex
 ifeq ($(USEOMP),true)
  FFLAGS += -KOMP
 endif
 ifeq ($(USEMPI),true)
  FC      = mpifrt
 endif
endif

### Mac Intel fortran
ifeq ($(ENVFLAGS),MacIfort)
 VER := $(shell sw_vers -productVersion | sed -r 's/^([0-9]*)\..*/\1/')
 ifeq ($(shell [ $(VER) -ge 14 ] && echo true), true)
  FFLAGS1 = -Wl,-ld_classic
 else
  FFLAGS1 =
 endif
 FC       = ifort
 IntelFortran = true
 SRCS8  += mdp-uni90.f
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -fpconstant -g -traceback -fpe0 -fp-stack-check -check uninit -check bounds $(FFLAGS1)
 else
  FFLAGS = -O3 -fpconstant $(FFLAGS1)
 endif
 ifeq ($(USEOMP),true)
   FFLAGS  += -qopenmp
    LD = $(FC) -qopenmp-link=static
 endif
 ifeq ($(USEMPI),true)
  FC       = mpif90
  FFLAGS  += -I/usr/local/lib -I/usr/local/include $(FFLAGS1)
  LDLIBS  += -lmpi_f90 -lmpi_f77 -lmpi -lm
  LDFLAGS += -L/usr/local/lib  -Wl,-commons,use_dylibs
 endif
endif

### Mac PGI Fortran
### PGI Fortran cannot accept CR+LF return code
### You have to change all return codes in input files, including xsdir.jnd
ifeq ($(ENVFLAGS),MacPGIfort)
 FC = pgfortran
 SRCS8  += mdp-uni90.f
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -g
 else
  FFLAGS = -O2 -Mipa=fast
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -mp
 endif
 ifeq ($(USEMPI),true)
  FC = mpif90
  # FFLAGS += -Mmpi=mpich # for compiling w/ pgfortran
 endif
endif

### Mac gfortran
ifeq ($(ENVFLAGS),MacGfort)
 OtherGfortran = true
 LDLIBS += -static-libgfortran -static-libgcc -static-libquadmath
 SRCS8    += mdp-uni90.f
 FC        = gfortran
 ifeq ($(DEBUG),true)
  FFLAGS = -O0 -g -fdefault-double-8 -fdefault-real-8 -fdollar-ok -std=legacy -ffpe-summary=none
 else
  FFLAGS = -O3 -fdefault-double-8 -fdefault-real-8 -fdollar-ok -std=legacy -ffpe-summary=none
 endif
 ifeq ($(USEMPI),true)
  FC       = mpifort
  FFLAGS  +=   $(shell mpifort --showme:compile)
  LDFLAGS  =   $(shell mpifort --showme:link)
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -fopenmp
 endif
endif

### Windows gfortran (4.8, 7.0 or later, MPI unchecked)
### You have to install "gfortran". make command is generally "mingw32-make"
### Download site for GCC including gfortran for Windows
### Please see /phits/document/Install_gfortran_Windows-**.pdf
ifeq ($(ENVFLAGS),WinGfort)
 OtherGfortran = true
 SRCS8   += mdp-win.f
 FC       = gfortran
 FFLAGS   = -fdefault-double-8 -fdefault-real-8 -fdollar-ok -std=legacy
 ifeq ($(USEMPI),true)
   FFLAGS  += -I"C:\Program FIles\MPICH2\include"
   LDFLAGS  = -I"C:\Program FIles\MPICH2\include"
   LDLIBS  += -L"C:\Program FIles\MPICH2\lib" -lfmpich2 -lfmpich2g -lmpi
 endif
 ifeq ($(USEOMP),true)
  FFLAGS += -fopenmp
 endif
endif

### Cygwin gfortran (4.8, 7.0 or later, MPI unchecked, Cygwin-32bit is not supported, use Cygwin-64bit.)
ifeq ($(ENVFLAGS),WinCygGfort)
 OtherGfortran = true
 SRCS8     += mdp-win.f
 FC         = gfortran
 FFLAGS     = -O0 -fdefault-double-8 -fdefault-real-8 -fdollar-ok
 ifeq ($(USEOMP),true)
  FFLAGS += -fopenmp
 endif
 ifeq ($(USEMPI),true)
  # FFLAGS  += -I/usr/local/lib -I/usr/local/include
  # LDFLAGS  = -L/usr/local/lib -Wl,-commons,use_dylibs
  # LDLIBS   =  -lmpi_f90 -lmpi_f77 -lmpi -lm
 else
  LDFLAGS  = -Wl,--stack,8388608
 endif
endif

#=======================================================================
# modules
#=======================================================================
SRCS0 = moddas_character.f moddas_fragdata.f moddas_material.f \
	moddas_mesh.f moddas_multiplier.f moddas_region.f \
	moddas_region_mtrg.f moddas_repeated_collisions.f moddas_source.f \
	moddas_tally.f moddas_variance_reduction.f moddas_ggs.f \
	moddas_bends.f \
	moddas.f \
	mod_ompparallel.f mod_counter.f mod_rsouin.f partmod.f talmod.f \
	mmbankmod.f membankmod.f ggbankmod.f \
	ggmbankmod.f eventtalmod.f restalmod.f \
	ngsdatamod.f smm.f gammod.f ggmarraymod.f ggm05mod.f \
	egs5mod.f \
	sumtallymod.f anatallymod.f tetramod.f qmd_cood2_mod.f \
	sangelmod.f usrtalmod.f  neutrinomod.f ndata2mod.f magmod.f \
	fragdatamod.f elcmod.f dedx_file.f latdatamod.f levdenmod.f \
	tallsm4.f fissionmod.f lafdatamod.f cosmicmod.f \
	elemdatamod.f itsmode.f etsart.f tsedmod.f liboutmod.f\
	charvarmod.f cvaloutmod.f \
	udm_parameter.f udm_utility.f \

#=======================================================================
# modules Fortran 90
#=======================================================================
SRCS0_90 = 	$(wildcard udm_part_*.f90) $(wildcard udm_int_*.f90)

#=======================================================================
# for machine dependent or user defined source and analysis
#=======================================================================
SRCS1 = \
	usrsors.f  usrmgf1.f  usrmgf3.f  anal-002.f \
	usrdfn1.f  usrdfn2.f  usrdfn3.f usrresol.f  usranatal.f \
	talls_4dtrack.f

#=======================================================================
# user defined model
#=======================================================================
SRCS1_90 = udm_manager.f90 decay_utility.f90

#=======================================================================
# for param.inc
#=======================================================================
SRCS2 = \
	analyz.f   celimp.f   dataup.f   getflt.f   magtrs.f   \
	nreac.f    ovly12.f   ovly13.f   partrs.f   range.f    \
	read00.f   read01.f   read02.f   sors.f     talls00.f  \
	talls01.f  talls02.f  talls03.f  talls04.f  talls05.f  \
	talls06.f  talls07.f  tallsm1.f  tallsm2.f  tallsm3.f  \
	update.f   wrnt12.f   wrnt13.f   read03.f   marscg.f   \
	ggs00.f    ggs01.f    ggs02.f    ggs03.f    wrnt10.f   \
	geocntl.f  ggm01.f    ggm02.f    ggm03.f    ggm04.f    \
	ggm05.f    ggm06.f    ggm07.f    ggm08.f    ovly14.f   \
	ovly15.f   usrelst1.f usrelst2.f usrmgt1.f  usrmgt2.f  \
	talls08.f  statistic.f itrminmax.f ggs04.f  a-angel.f  \
	talls09.f  errcheck.f errwrite.f ggm09.f  multiplier.f

#=======================================================================
# for new package
#=======================================================================
SRCS3 = \
	main.f     dklos.f   ncasc.f     nelst.f    nevap.f    \
	sdml.f     gem.f     gemset.f    utl01.f    \
	utl02.f    jbook.f   masdis.f    atima01.f  atima02.f  \
	atima03.f  fismul.f  kurotama0.f usrtally.f incelf.f   \
	incelfin.f inclp.f   inclin.f    dwbain.f   dwbaD.f    \
	dwbaE.f    photnucl.f  aama.f    muoncap.f ompdump.f   \
	muonint.f  risource.f dpa.f \
	sangel.f   gamma2mumu.f scinucl.f gamma1.f

#=======================================================================
# for old package
#=======================================================================
SRCS4 = \
	bert.f     bertin.f   bert-bl0.f bert-bl1.f bert-bl2.f \
	utlnmtc.f  gamlib.f   erupin.f   erup.f     fissn.f    \
	isobert.f  isodat.f   randmc.f   energy.f   ndata01.f  \
	mars00.f   mars01.f   mars02.f   mars03.f   mars04.f

#=======================================================================
# for JAM
#=======================================================================
SRCS5 = \
	jamin.f    jam.f      jamdat.f   jamcoll.f  jamdec.f   \
	jamcross.f jampdf.f   jamsoft.f  jamhij.f   jamhard.f  \
	jambuu.f   jamana.f   pyjet.f    pythia.f   pysigh.f

#=======================================================================
# for JQMD
#=======================================================================
SRCS6 = \
	qmd00.f    qmdcoll.f  qmddflt.f  qmdgrnd.f  qmdinit.f   \
	qmdmfld.f                                               \
	qmd00_qmd-nr.f     qmdcoll_qmd-nr.f    qmdgrnd_qmd-nr.f \
	qmdmfld_qmd-nr.f qmdscm_qmd-nr.f

#=======================================================================
# for ANGEL except for utl03.f, a-func.f, a-utl00.f
#=======================================================================
SRCS7 = \
	utl03.f    a-func.f   a-utl00.f     \
	a-main0.f  a-main1.f  a-hsect.f  a-line.f   a-wtext.f

#=======================================================================
# for restart function and sumtally and anatally options
#=======================================================================
SRCS9 = \
	restart.f resutl.f \
	resttrack.f restcross.f restheat.f restdeposit.f restyield.f \
	restlet.f   restsed.f   restdpa.f  restproduct.f resttime.f  \
	restdeposit2.f reststar.f restpoint.f restadjnt.f restwwg.f  \
	restvlm.f \
	sumtally.f \
	sumtaltrack.f sumtaldeposit.f \
	sumtalcross.f sumtaldeposit2.f sumtaldpa.f sumtalheat.f \
        sumtallet.f   sumtalproduct.f  sumtalsed.f sumtalstar.f \
	sumtaltime.f  sumtalyield.f sumtalpoint.f sumtaladjnt.f \
	sumtalwwg.f   sumtalvlm.f \
        anatally.f anainp.f anatal_rearrange.f \
	anataltrack.f anataldeposit.f anatalyield.f anatalproduct.f \
	anatalcross.f anataldpa.f anatalstar.f anataltime.f \
	anatallet.f anatalsed.f anatalpoint.f

#=======================================================================
# for EGS5
#=======================================================================
SRCS10 = \
	egs5.f  egs5init.f      egs5pfpl.f      egs5pcoll.f     \
	egs5efpl.f      egs5ecoll.f     egs5ede.f       egs5annih.f     \
	egs2phits.f     egs5esteps0.f	egs5edxde.f	egs5init2.f	\
	egs5elast.f     egs5ededx.f   egs_elmgfd.f  etsmode.f etsmode02.f \
    kurbuc.f

#=======================================================================
# for BITMAP
#=======================================================================
SRCS11 = \
	bitmap_main.f bitmap_a-main0.f bitmap_a_angel.f bitmap_util.f  \
	vtk_main.f vtk_gshow.f vtk_util.f

#=======================================================================
# for f90
#=======================================================================
SRCSf90 = \
	cosmicray.f90

#=======================================================================
# for hdf5
#=======================================================================
SRCShdf5 = tetrahdf5.F90

OBJS = $(SRCS0:.f=.o) \
       $(SRCS0_90:.f90=.o) \
       $(SRCS1:.f=.o) \
       $(SRCS1_90:.f90=.o) \
       $(SRCS2:.f=.o) \
       $(SRCS3:.f=.o) \
       $(SRCS4:.f=.o) \
       $(SRCS5:.f=.o) \
       $(SRCS6:.f=.o) \
       $(SRCS7:.f=.o) \
       $(SRCS8:.f=.o) \
       $(SRCS9:.f=.o) \
       $(SRCS10:.f=.o) \
       $(SRCS11:.f=.o) \
       $(SRCSf90:.f90=.o) \
       $(SRCShdf5:.F90=.o)

LD = $(LINK.f)
ifeq ($(ENVFLAGS),MacIfort)
  LD = $(FC) -static-intel $(FFLAGS)
  ifeq ($(USEOMP),true)
    LD = $(FC) -qopenmp-link=static  $(FFLAGS)
  endif
endif
# LD = $(FC)

ifeq ($(ENVFLAGS),WinGfort)
 ifeq ($(USEMPI),true)
$(TARGET) : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)
 else
$(TARGET) : $(OBJS)
	$(LD) -o $@ $^
 endif
else
$(TARGET) : $(OBJS)
	$(LD) -o $@ $^ $(LDLIBS)
endif

# According to the description [10.2 Catalogue of Built-In Rules]
# for 'Linking a single object file' in the GNU make manual,
# "... when there is no object file whose  name derives from the
# executable file name, you must write an explicit recipe for linking."
# though there're some built-in rules...

# .f.o:
#	$(FC) -c $(FFLAGS) $(INCLUDES) $(OPTFLAGS) $*.f $(LIBS)
# For 'Compiling Fortran and Ratfor programs' there's no need
# to write an explicit rule, because the following rules suffice.
##### The GNU make built-in rules
#%.o: %.f
# #  recipe to execute (built-in):
#         $(COMPILE.f) $(OUTPUT_OPTION) $<
# # default
# COMPILE.f = $(FC) $(FFLAGS) $(TARGET_ARCH) -c

### default rule for a module (nothing done, just a dependency)
# c.f. http://lagrange.mechse.illinois.edu/f90_mod_deps/
#      http://qiita.com/chirimen/items/af8479d8c0c7cf545b1e
%.mod : %.f %o %.f90 %.F90
	@:

### clean-up and recreate the file which describes dependencies
#     for all fortran source files w/ extension '.f'.
#     Note the use of 'makedepfort' with option '-m', which removes
#     the dependency on 'mpif.h'.
DEPFILE = make.depend
include $(DEPFILE)
ifdef DEVEL
 depend $(DEPFILE):
	@$(RM) -f $(DEPFILE)
	@makedepfort *.f *.f90 *.F90> $(DEPFILE)
endif

##### Special (explicit) compile rules for some object files.
gammod.o gamma1.o :
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f

ifeq ($(DEBUG),true)
 usrtally.o:
 ifeq ($(ENVFLAGS),LinIfort)
	$(FC) -c $(FFLAGS) -check bounds -gen-interfaces -warn interfaces $*.f
 endif
 ifeq ($(ENVFLAGS),LinNVIDIA)
	$(FC) -c $(FFLAGS) -Mbounds $*.f
 endif
endif

ifeq ($(OtherGfortran),true)
  read00.o read03.o risource.o:
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f
endif
ifeq ($(IntelFortran),true)
  restdeposit.o:
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f
  sumtaltrack.o sumtaldeposit.o sumtalcross.o sumtaldeposit2.o \
  sumtaldpa.o sumtalheat.o sumtallet.o sumtalproduct.o sumtalsed.o \
  sumtalstar.o sumtaltime.o sumtalyield.o sumtalpoint.o:
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f
  ifeq ($(USEOMP),true)
    aama.o:
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f
  endif
endif

tetrahdf5.o:
	$(FC) -c $(FFLAGS) $(INCLUDES) $(OPTFLAGS) $*.F90

### f90
ifeq ($(ENVFLAGS),LinNVIDIA)
%.o: %.f90
	$(FC) -c $(filter-out -Mfixed,$(filter-out -O%,$(FFLAGS))) -O3 $*.f90
%.o: %.F90
	$(FC) -c $(filter-out -Mfixed,$(filter-out -O%,$(FFLAGS))) -O3 $*.F90
else
%.o: %.f90
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O3 $*.f90
%.o: %.F90
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O3 $*.F90
endif

### Mac Intel fortran
ifeq ($(ENVFLAGS),MacIfort)
atima01.o talls06.o:
	$(FC) -c $(filter-out -O%,$(FFLAGS)) -O0 $*.f
 endif

clean:
	@$(RM) -f $(OBJS) *.mod
ifdef DEVEL
	@$(RM) -f $(DEPFILE)
endif

.PHONY: depend
