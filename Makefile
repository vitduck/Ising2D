F90 = ifort 
FFLAGS = -O2 -i4 -r8 -warn -CB -traceback -nogen-interfaces\
		 -vec-report1 -mcmodel medium -shared-intel

SRC = ising.f90 main.f90
OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

EXE = mc.x

# suffix rules 
%.o: %.f90
	$(F90) $(FFLAGS) -c $<

# make rules 
mc: $(OBJ) 
	$(F90) $(FFLAGS) $(OBJ) -o ${EXE}

clean: 
	-rm -f ${OBJ} ${MOD}
