#
# Makefile for mdlovofit
#
# This make file will try to compile mdlovofit with gfortran
# For doing this, just type
#
#          make 
#
# If you want to compile with some specific fortran compiler, you 
# must change the line below to the path of your fortran compiler. 
#
# Version 13.060
#
FORTRAN = gfortran
#FORTRAN = ifort
#
# Change "AUTO" to the fortran command you want.
#
# Change the flags of the compilation if you want:
#
FLAGS = -O3 
ifeq ($(MAKECMDGOALS),devel)
FLAGS = -Wall -fcheck=bounds -g -fbacktrace -ffpe-trap=zero,overflow,underflow
endif
#
# Source files:
SRC=./src
# Object files:
OBJ=./obj
# Executable files
BIN=./bin

all : $(BIN)/mdlovofit $(BIN)/user_field.tcl warning

devel : $(BIN)/mdlovofit $(BIN)/user_field.tcl warning 

#
# mdlovofit
#
$(BIN)/mdlovofit : $(OBJ)/common.o $(OBJ)/lib.o $(OBJ)/mdlovofit.o $(BIN)/user_field.tcl
	$(FORTRAN) $(FLAGS) -o $(BIN)/mdlovofit $(OBJ)/common.o $(OBJ)/mdlovofit.o $(OBJ)/lib.o

$(OBJ)/lib.o : $(SRC)/lib.f90 
	$(FORTRAN) $(FLAGS) -std=legacy -c -o $(OBJ)/lib.o $(SRC)/lib.f90

$(OBJ)/mdlovofit.o : $(SRC)/mdlovofit.f90
	$(FORTRAN) $(FLAGS) -c -o $(OBJ)/mdlovofit.o $(SRC)/mdlovofit.f90

$(OBJ)/common.o : $(SRC)/common.f90 
	$(FORTRAN) $(FLAGS) -c -o $(OBJ)/common.o $(SRC)/common.f90

$(BIN)/user_field.tcl : $(SRC)/user_field.tcl
	cp -f $(SRC)/user_field.tcl $(BIN)/user_field.tcl

warning : 

	@echo " ------------------------------------------------------ " 
	@echo " Compilation finished successfully. "
	@echo " IMPORTANT: Add the executable directory to your path. "
	@echo " Current executable files directory: `pwd`/bin "
	@echo " ------------------------------------------------------ " 

clean :
	rm -f $(OBJ)/* 
	rm -f ./*.mod


