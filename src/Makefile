IDIR =../include
IDIR2=../third-party/eigen
IDIR3=../third-party/libigl
INCLUDES  = -I./ -I${IDIR} -I${IDIR2} -I${IDIR3}
LDIR =../lib
LIB = -lFloatTetwild -lgeogram -lfmt -lgmp -fopenmp -lpthread

SRC = $(wildcard *.cpp)

default: build
	echo "Start Build"

# CXX=g++
CXX=mpic++

EXE = main

CXXFLAGS = -O3 -std=c++11 ${INCLUDES}

#The value of the make variable VPATH specifies a list of directories that make should search. 
vpath %.cpp $(sort $(dir $(SRC)))

OBJ = $(notdir $(SRC:.cpp=.o))  

$(warning $(OBJ) $(EXE) $(sort $(dir $(SRC))))

build: $(EXE)

$(EXE): $(OBJ) 
	${CXX} $(OBJ) -L$(LDIR) $(LIB) -o $(EXE)

clean:
	rm -f *.o *.cuda *.host *.a *.h *.tmp *.x ${EXE}

# Compilation rules

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $(notdir $@)

test: $(EXE)
	mpirun -np 1 ./$(EXE)
