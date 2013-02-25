# C source code
CSRC	= two_photos_psf_estim.c two_photos_psf_estim_main.c image.c ls.c io_pgm.c
# C++ source code

#Third-Party Software (ORSA-homography)
CXXSRC	= ./third_party/lib_orsa_homography.cpp \
	  ./third_party/OrsaHomography_20120515/src/libOrsa/conditioning.cpp \
	  ./third_party/OrsaHomography_20120515/src/libOrsa/homography_model.cpp \
	  ./third_party/OrsaHomography_20120515/src/libOrsa/orsa_model.cpp  \
	  ./third_party/OrsaHomography_20120515/src/extras/sift/demo_lib_sift.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/sift/domain.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/sift/filter.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/sift/library.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/sift/splines.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/libNumerics/matrix.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/libNumerics/numerics.cpp \
	  ./third_party/OrsaHomography_20120515/src/extras/libNumerics/vector.cpp

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
COPT = -g


# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra -Wno-write-strings -ansi  -Werror
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi -Werror
# link flags
LDFLAGS	=  -lm -lfftw3f -lblas -llapack  -lstdc++ -lgomp

# include dir for fftw, lapack, blas
INCDIR=-I/usr/local/include \
       -I./third_party/OrsaHomography_20120515/src/demo \
       -I./third_party/OrsaHomography_20120515/src/extras \
       -I./third_party/OrsaHomography_20120515/src/extras/libMatch \
       -I./third_party/OrsaHomography_20120515/src

LIBDIR=-L/usr/local/lib

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

