#=======================================================================
#  FIRST-CUT MAKEFILE FOR ATHENA CODE
#
#-------------------  object files  ------------------------------------
OBJ =   add1_2name.o   \
	binary_dump.o \
	data_output.o \
	eigensystems_InConsVars.o \
	eigensystems_InPrimVars.o \
	init_dt.o \
	init_grid_block.o \
	integrate_1step.o \
	lr_states.o \
	main.o \
	printd.o \
	roe_fluxes.o \
	set_bval_arrays.o \
	set_ghost_zones.o \
	shkset.o \
	twoibw.o \
	shu_osher.o \
	linear_wave.o
#-------------------  macro definitions  -------------------------------
EXEDIR = ../bin/
# macros for sun or SGI
CC  = gcc
OPTS = -O3 -Wall
LDR = gcc -O3
LIB = -lm

#--------------------  implicit rules  ---------------------------------
.c.o:
	${CC} ${OPTS} -c $<
#---------------------  targets  ---------------------------------------
help:
	@echo Type 'make compile' to generate ATHENA executable
	@echo Type 'make clean'   to remove *.o files
#-----------------------------------------------------------------------
clean:
	rm -f *.o Makedepend
	touch Makedepend
#-----------------------------------------------------------------------
compile: ${OBJ}
	${LDR} -o ${EXEDIR}athena ${OBJ} ${LIB}
#-----------------------------------------------------------------------
depend:
	rm -f Makedepend
	${CC} -M ${OPTS} ${OBJ} > Makedepend

include Makedepend
