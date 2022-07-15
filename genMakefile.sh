#!/bin/bash
LMKL=" -lmkl_intel_lp64 -lmkl_core -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack95_lp64"

CFLAGS=" -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict"
LFLAGS=" -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict "$LMKL

CFLAGS_OPT=" -c -traceback -O3 -ipo"
LFLAGS_OPT=" -traceback -O3 -ipo "$LMKL

FoBiS.py build -s ./src -m Makefile -compiler intel -fc "mpiifort" -lflags "$LFLAGS" -cflags "$CFLAGS"

FoBiS.py build -s ./src -m Makefile.opt -compiler intel -fc "mpiifort" -lflags "$LFLAGS_OPT" -cflags "$CFLAGS_OPT" --obj_dir ./obj_opt --mod_dir ./mod_opt

#build all fortran programs from src
#please note that this is a dirty trick, we need a better solution
#maybe to gig in FoBiS.py options/code
echo "all: \$(addprefix \$(DEXE),\$(EXES))" >> Makefile
