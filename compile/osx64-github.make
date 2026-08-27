#------------------------------------------------------------------------------
# This one is for the GitHub Actions CI.
# Compilers information
#     FC: Fortran Compiler
#     CC: C Compiler
#     LD: Linker
#     FFLAGS : Fortran flags
#     CFLAGS : C Flags
#     LDFLAGS: Linker flags
#     COMPILE_LIB_NAME : Name of the compiled library file
#     INSTALL_LIB_NAME : Name of the installed library file
#------------------------------------------------------------------------------

FC := $(shell bash -c "compgen -c gfortran | sort -V | tail -n1")
CC = gcc
LD := $(FC)

FFLAGS=-fpic -fno-second-underscore -std=legacy -ffixed-line-length-none -O2 -march=native
CFLAGS=-fpic -O2 -march=native
LDFLAGS=-shared

COMPILE_LIB_NAME=libirbem.$(OS).$(ENV).so
INSTALL_LIB_NAME=libirbem.so
