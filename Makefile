# C source code
CSRC	= two_photos_psf_estim.c two_photos_psf_estim_main.c homography.c image.c ls.c io_pgm.c
# C++ source code
CXXSRC	= ./HomoOrsa/libOrsaHomography/lib_orsa_homography.cpp ./HomoOrsa/libOrsa/conditioning.cpp \
	./HomoOrsa/libOrsa/homography_model.cpp ./HomoOrsa/libOrsa/orsa_model.cpp  \
	./HomoOrsa/extras/sift/demo_lib_sift.cpp	 ./HomoOrsa/extras/sift/domain.cpp ./HomoOrsa/extras/sift/filter.cpp \
	./HomoOrsa/extras/sift/library.cpp ./HomoOrsa/extras/sift/splines.cpp ./HomoOrsa/extras/libNumerics/matrix.cpp \
	./HomoOrsa/extras/libNumerics/numerics.cpp ./HomoOrsa/extras/libNumerics/vector.cpp


# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= two_photos_psf_estim

default	: $(BIN)

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	=  -lm -lfftw3f -lblas -llapack  -lstdc++ -lgomp

# include dir for fftw, lapack, blas
INCDIR=-I/opt/local/include -I./HomoOrsa -I./HomoOrsa/extras
LIBDIR=-L/opt/local/lib

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif


# partial compilation of C source code
%.o: %.c %.h
	$(CC) $(INCDIR) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp 
	$(CXX) $(INCDIR) -c -o $@  $< $(CXXFLAGS)


# link all the opject code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) $(LIBDIR) -o $@ $(OBJ) $(LDFLAGS)

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN)


tar:
	tar -zc -f two_photos_psf_estim.tar.gz  --exclude .svn --exclude two_photos_psf_estim.tar.gz -C../ two_photos_psf_estim/ 
