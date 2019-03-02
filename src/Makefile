# GNU Makefile
.SUFFIXES: .o .f90
FC = gfortran
FCFLAGS = -O3
#FCFLAGS = -fbounds-check
PROGRAM = ../mkXtal
OBJ = ../../stringParseMods/lineParse.o ../../stringParseMods/domMod/domtype.o matrix.o tessel.o mkXtal.o

${PROGRAM}:${OBJ}
	${FC} -o ${PROGRAM} ${OBJ} ${FCFLAGS}

# to clean up stuff
clean:
	rm -f *.o *.mod *.MOD
veryclean: clean
	rm -f *~ $(PROGRAM)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<


