# Makefile for compiling Compose
# ==============================
#
#
#
# Options:
# --------
#
#
# BOUNDCHECK=(0|1)
#
# Note: Switch on bound checking for Fortran part.
#
# HDF5=(0|1)
#
# Enable writing in HDF5-format for the data tables



NAME = compose
EXEC = $(NAME)

# default flag settings (BOUNDCHECK=0, HDF5 = 0)

BOUNDCHECK = 0
HDF5 = 0
ifeq ($(BOUNDCHECK),1)
   FC_FLAGS_BOUNDCHECK = #-fbounds-check
endif

ifeq ($(HDF5),1)
   HDF5_LIB = -lhdf5
   HDF5_C = -Dhdf5
   FC = h5pfc
   CC = h5pcc
   LINK = $(FC)
else
   HDF5_LIB =
   HDF5_C =
   FC = gfortran
   CC = gcc
   LINK = $(FC)
endif

DEBUG = -g -DDEBUG -debug -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow
FC_FLAGS = -O3  -cpp -pg -fopenmp #$(DEBUG) #$(HDF5_C)
CC_FLAGS = -g -O3 -cpp

LD_LIB = $(HDF5_LIB)


ifeq ($(HDF5),1)
  SRC_F = composemodules.f90 \
          out_to_json.f90 \
          hdf5compose.f90 \
	  get_tables.f90 \
          compose.f90
  SRC_C = hdf5writecompose.c \
	  hdf5readcompose.c
else
  SRC_F = composemodules.f90 \
          out_to_json.f90 \
	  get_tables.f90 \
          compose.f90
  SRC_C =
endif


OBJ_C := $(SRC_C:.c=.o)
OBJ_F := $(SRC_F:.f90=.o)

$(NAME):  $(OBJ_F) $(OBJ_C)
	rm -f $(EXEC);
	@echo building compose;
	$(LINK)  $(FC_FLAGS) -o $(EXEC) $(OBJ_F) $(OBJ_C) $(LD_LIB)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $< -o $@


%.o : %.c
	$(CC) -c $(CC_FLAGS) $< -o $@


clean:
	rm -f *.o *.mod compose
