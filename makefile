
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
BIN		 = ./bin/
NSSRC		 = ./NS/

OBJSFILE		 :=$(BIN)main.o 
OBJSFILE		 +=$(BIN)dataProcess.o 
OBJSFILE		 +=$(BIN)navier.o 
OBJSFILE		 +=$(BIN)readparam.o
OBJSFILE		 +=$(BIN)tools.o
OBJSFILE		 +=$(BIN)scalar.o
OBJSFILE		 +=$(BIN)bc.o
OBJSFILE		 +=$(BIN)pressure.o
OBJSFILE		 +=$(BIN)gradient.o
OBJSFILE		 +=$(BIN)velocity.o

MYLINKER = mpicxx ${PCC_LINKER_FLAGS} ${CPPFLAGS}


CLEANFILES       = rhs.vtk solution.vtk


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

play: 
	mpiexec -np 4 ./cycas2

cl:
	${RM} $(BIN)*.o

cycas2:  $(OBJSFILE) chkopts
	-${MYLINKER} -o $@ $(OBJSFILE)  ${PETSC_LIB}

$(BIN)main.o: main.cpp 
	$(PETSC_CXXCOMPILE) -o $@ $^ 

$(BIN)dataProcess.o: dataProcess.cpp
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


include ${PETSC_DIR}/lib/petsc/conf/test
