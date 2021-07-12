# Makefile

# MACHINE = debug
MACHINE = mac
# MACHINE = rocks
# MACHINE = solidmatters
# MACHINE = ubuntu
# MACHINE = ubuntu-debug
# MACHINE = millipede-intel
# MACHINE = millipede


ifeq ($(MACHINE), debug)
CC = gcc
CC_FLAGS = -O0 -g -w
LAPACK = -framework veclib
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), mac)
CC = gcc
CC_FLAGS = -O3
LAPACK = -llapack
#LAPACK = -framework veclib
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), obsolete)
CC = gcc
CC_FLAGS = -O3
LAPACK = -framework Accelerate
FFTW = -lfftw3 -lm -L/opt/local/lib
FFTWI = -I/opt/local/include
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), rocks)
CC = gcc
CC_FLAGS = -O3 -lm -L/share/apps/intel/composer_xe_2015/mkl/lib/intel64 -I/share/apps/intel/composer_xe_2015/mkl/include -lpthread
LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFTW = -lfftw3
endif

ifeq ($(MACHINE), solidmatters)
CC = icc
FFTWI = -I/software/apps/fftw/intel/3.3.3/include
FFTW = -lfftw3 -L/software/apps/fftw/intel/3.3.3/lib
CC_FLAGS = -B -O3 -lm
LAPACK = -llapack
LAPACKL = -L/usr/lib64/atlas
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), ubuntu)
CC = gcc
CC_FLAGS = -O3 -L/usr/lib/x86_64-linux-gnu -lm
LAPACK = -llapack
FFTW = -lfftw3
endif

ifeq ($(MACHINE), ubuntu-debug)
CC = gcc
CC_FLAGS =  -O0 -g -w -L/usr/lib/x86_64-linux-gnu -lm
LAPACK = -llapack
FFTW = -lfftw3
endif

ifeq ($(MACHINE), millipede-intel)
CC = icc
FFTWI = -I/cm/shared/apps/fftw/gcc/64/3.2.2/include
FFTW = -lfftw3 -L/cm/shared/apps/fftw/gcc/64/3.2.2/lib
CC_FLAGS = -B -O3 -lm
MKLDIR = /cm/shared/apps/intel-cluster-runtime/2.1/mkl/10.2.0.013/lib/em64t
LAPACK = -L$(MKLDIR) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath=$(MKLDIR)
CPPFLAGS = $(FFTWI)
endif

ifeq ($(MACHINE), millipede)
CC = gcc
FFTWI = -I/cm/shared/apps/fftw/gcc/64/3.2.2/include
FFTW = -lfftw3 -L/cm/shared/apps/fftw/gcc/64/3.2.2/lib
CC_FLAGS = -B -O3 -lm
LAPACK = -llapack -lblas -lgcc
LAPACKL = -L/cm/shared/apps/lapack/gcc/64/3.2.1 -L/cm/shared/apps/blas/gcc/1/lib64
CPPFLAGS = $(FFTWI)
endif


all: VibroNISE Spectra

VibroNISE:			VibroNISE.o SubsMod.o ParmsHandleMod.o VOverlapsMod.o BasisMod.o ToolsMod.o RandomMod.o
	$(CC) -o VibroNISE	VibroNISE.o SubsMod.o ParmsHandleMod.o VOverlapsMod.o BasisMod.o ToolsMod.o RandomMod.o $(LAPACK) $(CC_FLAGS)

Spectra:			Spectra.o ParmsHandleMod.o ToolsMod.o
	$(CC) -o Spectra	Spectra.o ParmsHandleMod.o ToolsMod.o $(LAPACK) $(FFTW) $(CC_FLAGS)

VibroNISE.o:			VibroNISE.c SubsMod.h ParmsHandleMod.h VOverlapsMod.h BasisMod.h ToolsMod.h RandomMod.h ParmsMod.h GlobalsMod.h

Spectra.o:			Spectra.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

SubsMod.o:			SubsMod.c SubsMod.h BasisMod.h ToolsMod.h RandomMod.h ParmsMod.h GlobalsMod.h

ParmsHandleMod.o:		ParmsHandleMod.c ParmsHandleMod.h ToolsMod.h ParmsMod.h GlobalsMod.h

VOverlapsMod.o:			VOverlapsMod.c VOverlapsMod.h GlobalsMod.h

BasisMod.o:			BasisMod.c BasisMod.h ToolsMod.h GlobalsMod.h

ToolsMod.o:			ToolsMod.c ToolsMod.h GlobalsMod.h

RandomMod.o:			RandomMod.c RandomMod.h

clean:
	rm *.o VibroNISE Spectra
