########## User-definable stuff ##########
#
###Compiler and compilation options
COMP = mpicc
OPTIONS = -Wall -O3
#
### Behavioural flags
#Use double precision integer (enable in general)
DEFINEFLAGS += -D_LONGIDS
DEFINEFLAGS += -D_DEBUG
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
GSL_INC = -I/home/alonso/include
GSL_LIB = -L/home/alonso/lib
#FFTW
FFTW_INC = 
FFTW_LIB = 
#
########## End of user-definable ##########

#OPTIONS += -fopenmp

#ifeq ($(strip $(USE_SINGLE_PRECISION)),yes)
#DEFINEFLAGS += -D_SPREC
#LIB_FFTW = -lfftw3f_omp -lfftw3f
#else
#LIB_FFTW = -lfftw3_omp -lfftw3
#endif
LIB_FFTW = -lfftw3f_mpi -lfftw3f

OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(GSL_INC) $(FFTW_INC)
LIB_ALL = $(GSL_LIB) $(FFTW_LIB) -lgsl -lgslcblas $(LIB_FFTW) -lm

COMMONO = src/common.o
GRIDO = src/grid_tools.o
IOO = src/io.o
MAIN = src/main.c
OFILES = $(COMMONO) $(GRIDO) $(IOO)

EXE = DensTools

default : $(EXE)

%.o : %.c
	$(COMP) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(EXE) : $(OFILES)
	$(COMP) $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXE) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner : 
	rm -f *~ src/*.o src/*~ $(EXE)
