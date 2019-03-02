#!/bin/bash

F90="gfortran"
EXE="../mkXtal"
if [ -z $1 ]; then
   OP="-O3"
else
   OP="-fbounds-check"
fi
echo "compiling with options: " ${OP}

# source file list
libFils="lineParse.o domtype.o matrix.o "
appFils="tessel.o mkXtal.o"

# compilation to object files
for line in ${libFils}; do
   ${F90} ${OP} -c ${line%.*}.f90
done
for line in ${appFils}; do
   ${F90} ${OP} -c ${line%.*}.f90
done

# Linking
${F90} ${libFils} ${appFils} -o ${EXE}


