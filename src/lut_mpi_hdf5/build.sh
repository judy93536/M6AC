#!/bin/bash
#
# Build script for MODTRAN6 MPI+HDF5 LUT Generator
#
# Usage: ./build.sh

set -e  # Exit on error

MPICXX=/opt/mpich-3.3.2/bin/mpicxx
MODTRAN_HOME=/home/judy/MODTRAN6.0

echo "======================================================================"
echo "Building MODTRAN6 MPI+HDF5 LUT Generator"
echo "======================================================================"

# Check if MPICH exists
if [ ! -f "$MPICXX" ]; then
    echo "ERROR: MPICH compiler not found at $MPICXX"
    echo "Please install MPICH 3.3.2 or update MPICXX path in this script"
    exit 1
fi

# Check if MODTRAN6 exists
if [ ! -d "$MODTRAN_HOME" ]; then
    echo "ERROR: MODTRAN6 not found at $MODTRAN_HOME"
    echo "Please update MODTRAN_HOME path in this script"
    exit 1
fi

echo "Compiler: $MPICXX"
echo "MODTRAN:  $MODTRAN_HOME"
echo ""

# Compile
echo "Compiling lut_mpi_hdf5.cc..."

$MPICXX lut_mpi_hdf5.cc -o lut_mpi_hdf5 \
    -I$MODTRAN_HOME/developer/include/ \
    -L$MODTRAN_HOME/bin/linux \
    -Wl,-rpath='$ORIGIN',-rpath-link="$MODTRAN_HOME/bin/linux",-lmod6rt \
    $(pkg-config --cflags --libs hdf5) \
    -lcurl

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================================"
    echo "Build successful!"
    echo "======================================================================"
    echo "Executable: ./lut_mpi_hdf5"
    echo ""
    echo "To run:"
    echo "  export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:$MODTRAN_HOME/bin/linux:\$LD_LIBRARY_PATH"
    echo "  /opt/mpich-3.3.2/bin/mpirun -np 4 ./lut_mpi_hdf5 \\"
    echo "      template.json output.h5 \\"
    echo "      0.0 1.0 11 \\"  # albedo_min albedo_max n_albedo
    echo "      0.0 4.0 151"     # h2o_min h2o_max n_h2o
    echo ""
    echo "See README_LUT_MPI_HDF5.md for detailed usage"
else
    echo ""
    echo "======================================================================"
    echo "Build FAILED"
    echo "======================================================================"
    exit 1
fi
