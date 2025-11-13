# Makefile for MODTRAN6 MPI+HDF5 LUT Generator
#
# Build command: make
# Clean: make clean

# Compiler and flags
MPICXX = /opt/mpich-3.3.2/bin/mpicxx
MODTRAN_HOME = /home/judy/MODTRAN6.0

# Include and library paths
INCLUDES = -I$(MODTRAN_HOME)/developer/include
LIBDIRS = -L$(MODTRAN_HOME)/bin/linux
RPATH = -Wl,-rpath='$$ORIGIN',-rpath-link='$(MODTRAN_HOME)/bin/linux',-lmod6rt

# HDF5 flags (using pkg-config)
HDF5_FLAGS = $(shell pkg-config --cflags --libs hdf5)

# Additional libraries
LIBS = -lcurl

# Output executable
TARGET = lut_mpi_hdf5

# Source file
SRC = lut_mpi_hdf5.cc

# Build rule
all: $(TARGET)

$(TARGET): $(SRC)
	@echo "Compiling $(TARGET)..."
	$(MPICXX) $(SRC) -o $(TARGET) \
		$(INCLUDES) \
		$(LIBDIRS) \
		$(RPATH) \
		$(HDF5_FLAGS) \
		$(LIBS)
	@echo "Build complete: $(TARGET)"
	@echo ""
	@echo "To run:"
	@echo "  export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:$(MODTRAN_HOME)/bin/linux:\$$LD_LIBRARY_PATH"
	@echo "  /opt/mpich-3.3.2/bin/mpirun -np 4 ./$(TARGET) template.json output.h5 0.0 1.0 11 0.0 4.0 151"

clean:
	@echo "Cleaning..."
	rm -f $(TARGET)
	@echo "Clean complete"

.PHONY: all clean
