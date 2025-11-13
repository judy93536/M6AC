# M6AC Source Code

This directory contains compiled/external source code for M6AC.

## Contents

### lut_mpi_hdf5/

MODTRAN6 MPI+HDF5 LUT generator (C++)

**Purpose:** Fast parallel LUT generation for atmospheric correction

**Performance:**
- 11×151 = 1,661 MODTRAN runs in ~9 minutes (4 processes)
- Outputs single HDF5 file with complete LUT grid

**Build:**
```bash
cd lut_mpi_hdf5
make
# or
./build.sh
```

**Usage:**
```bash
export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:/home/judy/MODTRAN6.0/bin/linux:$LD_LIBRARY_PATH

/opt/mpich-3.3.2/bin/mpirun -np 32 ./lut_mpi_hdf5 \
    ../../data/ang20190624t230039_rdn_v2u1/ang20190624t230039_modtran_template.json \
    ../../data/ang20190624t230039_rdn_v2u1/luts/ang20190624t230039_lut.h5 \
    0.0 1.0 11 \
    0.0 4.0 151
```

See `lut_mpi_hdf5/README_LUT_MPI_HDF5.md` for detailed documentation.

## Requirements

- MODTRAN6.0+ with C++ API
- MPICH 3.3.2 (installed at `/opt/mpich-3.3.2/`)
- HDF5 C library
- pkg-config
- libcurl

## Git Policy

**Tracked:**
- ✓ Source code (`.cc`, `.cpp`, `.h`)
- ✓ Build scripts (`Makefile`, `build.sh`)
- ✓ Documentation (`.md`)

**Not Tracked:**
- ✗ Compiled binaries (added to `.gitignore`)
- ✗ Object files (`.o`)
- ✗ MODTRAN output files (`.tp6`, `.chn`, etc.)
