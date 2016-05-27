
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         = -g
FPPFLAGS         =
BIN		 = ./bin/
NSSRC		 = ./NS/

#-Llibpath -llibname
#libname is the lib filename without prefix "lib" .e.g 
#libmylib.a --> -lmylib
METIS_LIB	 = -lmetis

INCLUDE = -I./NS -I./

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

test:
	mpirun -np 2 ./cycas2  
	
play: 
	mpirun -np 2 ./cycas2  -log_summary

cl:
	${RM} $(BIN)*.o

cycas2:  $(OBJSFILE) chkopts
	-${MYLINKER} -o $@ $(OBJSFILE)  ${PETSC_LIB} $(METIS_LIB) $(CPPFLAGS)

$(BIN)main.o: $(NSSRC)main.cpp 
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)MPIStruct.o: $(NSSRC)MPIStructure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)RootStruct.o: $(NSSRC)RootStructure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)navier.o: $(NSSRC)navier.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)readparam.o: $(NSSRC)readparamfile.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)tools.o: $(NSSRC)tools.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)scalar.o: $(NSSRC)navier_scalar.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)bc.o: $(NSSRC)navier_bc.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)pressure.o: $(NSSRC)navier_pressure.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)gradient.o: $(NSSRC)navier_gradient.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)velocity.o: $(NSSRC)navier_velocity.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)geo.o: $(NSSRC)geometry.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

$(BIN)io.o: $(NSSRC)navierIO.cpp
	$(PETSC_CXXCOMPILE) -o $@ $^ $(CPPFLAGS) $(INCLUDE)

include ${PETSC_DIR}/lib/petsc/conf/test
