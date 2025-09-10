# Compiler setting
BUILDDIR = build
CC = gcc
MPI_CC = mpicc

# Base flag
CFLAGS  = -Wall -Wextra -I./include

# Optimization flags
OPT_CFLAGS = -Ofast -flto  -march=native
OMPFLAG = -fopenmp

# Source directory
SRC_DIR = src

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Source files
SRC_SERIAL       = $(SRC_DIR)/stencil_serial.c
SRC_SERIAL_NO_OMP = $(SRC_DIR)/stencil_serial_nomp.c
SRC_PARALLEL     = $(SRC_DIR)/stencil_parallel.c
SRC_PARALLEL_2   = $(SRC_DIR)/stencil_parallel_2.c

# Object files
OBJ_SERIAL        = $(BUILDDIR)/stencil_serial.o
OBJ_SERIAL_NO_OMP = $(BUILDDIR)/stencil_serial_nomp.o
OBJ_PARALLEL      = $(BUILDDIR)/stencil_parallel.o
OBJ_PARALLEL_NO_OMP = $(BUILDDIR)/stencil_parallel_nomp.o
OBJ_PARALLEL_2    = $(BUILDDIR)/stencil_parallel_2.o

# Executables
EXEC_SERIAL      = stencil_serial
EXEC_SERIAL_NO_OMP = stencil_serial_nomp
EXEC_PARALLEL    = stencil_parallel
EXEC_PARALLEL_NO_OMP = stencil_parallel_nomp
EXEC_PARALLEL_2  = stencil_parallel_2

# Default target
all: $(EXEC_SERIAL) $(EXEC_SERIAL_NO_OMP) $(EXEC_PARALLEL) $(EXEC_PARALLEL_NO_OMP) $(EXEC_PARALLEL_2)

# =================== SERIAL WITH OPENMP ===================
$(EXEC_SERIAL): $(OBJ_SERIAL)
	$(CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -o $@ $^

$(OBJ_SERIAL): $(SRC_SERIAL) | $(BUILDDIR)
	$(CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -c $< -o $@

# =================== SERIAL WITHOUT OPENMP ===================
$(EXEC_SERIAL_NO_OMP): $(OBJ_SERIAL_NO_OMP)
	$(CC) $(CFLAGS) $(OPT_CFLAGS) -o $@ $^

$(OBJ_SERIAL_NO_OMP): $(SRC_SERIAL_NO_OMP) | $(BUILDDIR)
	$(CC) $(CFLAGS) $(OPT_CFLAGS) -c $< -o $@

# =================== PARALLEL WITH MPI+OPENMP ===================
$(EXEC_PARALLEL): $(OBJ_PARALLEL)
	$(MPI_CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -o $@ $^

$(OBJ_PARALLEL): $(SRC_PARALLEL) | $(BUILDDIR)
	$(MPI_CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -c $< -o $@

# =================== PARALLEL WITH MPI ONLY ===================
# Compiles the *same source* as the hybrid version, but without the -fopenmp flag.
# The compiler will simply ignore the #pragma omp directives.
$(EXEC_PARALLEL_NO_OMP): $(OBJ_PARALLEL_NO_OMP)
	$(MPI_CC) $(CFLAGS) $(OPT_CFLAGS) -o $@ $^

$(OBJ_PARALLEL_NO_OMP): $(SRC_PARALLEL) | $(BUILDDIR)
	$(MPI_CC) $(CFLAGS) $(OPT_CFLAGS) -c $< -o $@

# =================== PARALLEL WITH MPI+OPENMP (2) ===================
$(EXEC_PARALLEL_2): $(OBJ_PARALLEL_2)
	$(MPI_CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -o $@ $^

$(OBJ_PARALLEL_2): $(SRC_PARALLEL_2) | $(BUILDDIR)
	$(MPI_CC) $(CFLAGS) $(OMPFLAG) $(OPT_CFLAGS) -c $< -o $@

# Clean target
clean:
	rm -rf $(BUILDDIR) $(EXEC_SERIAL) $(EXEC_SERIAL_NO_OMP) $(EXEC_PARALLEL) $(EXEC_PARALLEL_NO_OMP) $(EXEC_PARALLEL_2)

.PHONY: all clean