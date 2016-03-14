
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
BIN		 = ./bin/
NSSRC		 = ./NS/

#-Llibpath -llibname
#libname is the lib filename without prefix "lib" .e.g 
#libmylib.a --> -lmylib
METIS_LIB	 = -lmetis

OBJSFILE		 :=$(BIN)main.o 
OBJSFILE		 +=$(BIN)MPIStruct.o
OBJSFILE		 +=$(BIN)RootStruct.o
# ORIGINAL CYCAS 2 OBJECTS
OBJSFILE		 +=$(BIN)navier.o 
OBJSFILE		 +=$(BIN)readparam.o
OBJSFILE		 +=$(BIN)tools.o
OBJSFILE		 +=$(BIN)scalar.o
OBJSFILE		 +=$(BIN)bc.o
OBJSFILE		 +=$(BIN)pressure.o
OBJSFILE		 +=$(BIN)gradient.o
OBJSFILE		 +=$(BIN)velocity.o
OBJSFILE		 +=$(BIN)geo.o
OBJSFILE 		 +=$(BIN)io.o

MYLINKER = mpicxx ${PCC_LINKER_FLAGS} ${CPPFLAGS}


CLEANFILES       = rhs.vtk solution.vtk


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

play: 
	mpirun -np 4 ./cycas2

cl:
	${RM} $(BIN)*.o

cycas2:  $(OBJSFILE) chkopts
	-${MYLINKER} -o $@ $(OBJSFILE)  ${PETSC_LIB} $(METIS_LIB)

$(BIN)main.o: main.cpp 
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)MPIStruct.o: MPIStructure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)RootStruct.o: RootStructure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)navier.o: $(NSSRC)navier.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)readparam.o: $(NSSRC)readparamfile.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)tools.o: $(NSSRC)tools.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)scalar.o: $(NSSRC)navier_scalar.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)bc.o: $(NSSRC)navier_bc.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)pressure.o: $(NSSRC)navier_pressure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)gradient.o: $(NSSRC)navier_gradient.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)velocity.o: $(NSSRC)navier_velocity.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)geo.o: $(NSSRC)geometry.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)io.o: $(NSSRC)navierIO.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ 

include ${PETSC_DIR}/lib/petsc/conf/test
