# MODTRAN MPI+HDF5 LUT Generation System

## Overview

This system generates Look-Up Tables (LUTs) for atmospheric correction using MODTRAN 6 radiative transfer calculations. It uses MPI (Message Passing Interface) for parallel execution and stores results in HDF5 format for efficient data access.

**Key Features:**
- MPI-based parallelization (avoids FORTRAN thread-safety issues)
- HDF5 output format with complete spectral data
- Dynamic task distribution for load balancing
- ~1 second per MODTRAN run (tested with 4 processes)
- Supports arbitrary parameter grids (albedo × H2O)

**Performance:**
- Test: 6 runs in 6.14 seconds (4 MPI processes)
- Production LUT (11×151 = 1661 runs):
  - 4 processes: ~9.4 minutes
  - 8 processes: ~4.0 minutes
  - 16 processes: ~1.9 minutes

## System Requirements

### Software Dependencies
- MODTRAN 6.0.2r3 with C++ API
- MPICH 3.3.2 (or compatible MPI implementation)
- HDF5 with C library
- nlohmann/json (header-only, included in MODTRAN developer kit)
- libcurl

### Environment
- Linux (tested on Ubuntu with kernel 6.8.0-85)
- C++11 compatible compiler (g++)
- Sufficient disk space for output (~64KB per test, scales with grid size)

## Installation

### 1. Install MPICH 3.3.2

```bash
cd /tmp
wget https://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
tar -xzf mpich-3.3.2.tar.gz
cd mpich-3.3.2

./configure --prefix=/opt/mpich-3.3.2 \
    --enable-shared \
    --disable-static \
    --with-device=ch3:nemesis \
    FFLAGS="-fallow-argument-mismatch"

make -j 8
sudo make install
```

**Why MPICH 3.3.2?**
- System MPI may conflict with MATLAB's bundled MPI libraries
- MPICH provides a clean, isolated MPI environment
- Compatible with MODTRAN 6 FORTRAN core

### 2. Verify Installation

```bash
/opt/mpich-3.3.2/bin/mpirun --version
# Should show: HYDRA build details, MPICH version 3.3.2
```

## Compilation

### Build Command

```bash
cd /home/judy/MODTRAN6.0/developer/examples

/opt/mpich-3.3.2/bin/mpicxx lut_mpi_hdf5.cc -o lut_mpi_hdf5 \
    -I../include/ \
    -L../../bin/linux \
    -Wl,-rpath='$ORIGIN',-rpath-link='../../bin/linux',-lmod6rt \
    $(pkg-config --cflags --libs hdf5) \
    -lcurl
```

### Compilation Notes

- **mpicxx**: MPICH C++ compiler wrapper (handles MPI includes/libs automatically)
- **-I../include/**: MODTRAN C++ API headers
- **-L../../bin/linux**: MODTRAN shared library location
- **-Wl,-rpath='$ORIGIN'**: Runtime library search path (relative to executable)
- **pkg-config --cflags --libs hdf5**: HDF5 compiler/linker flags
- **-lcurl**: Required by MODTRAN library

### Expected Output

```
lut_mpi_hdf5 executable created successfully
```

## Configuration

### JSON Template Structure

The MODTRAN configuration template must be prepared with correct geometry and atmospheric parameters. Example:

```json
{
    "MODTRAN": [
        {
            "MODTRANINPUT": {
                "NAME": "ang20190624t230039_lut",
                "CASE": 0,
                "DESCRIPTION": "AVIRIS-NG LUT generation",
                "RTOPTIONS": {
                    "MODTRN": "RT_CORRK_FAST",
                    "IEMSCT": "RT_SOLAR_AND_THERMAL",
                    "IMULT": "RT_ISAACS_SCALED_AT_OBS"
                },
                "ATMOSPHERE": {
                    "CO2MX": 410.0,
                    "H2OSTR": 1.8,
                    "H2OUNIT": "g"
                },
                "GEOMETRY": {
                    "H1ALT": 4.830596380703965,
                    "H2ALT": 0.30000,
                    "LENN": 0,
                    "IDAY": 175,
                    "IPARM": 11,
                    "PARM1": 34.19660186767578,
                    "PARM2": 118.17120361328125,
                    "GMTIME": 23.010799407958984
                },
                "SURFACE": {
                    "SURFTYPE": "REFL_CONSTANT",
                    "SURREF": 0.5,
                    "GNDALT": 0.3
                },
                "SPECTRAL": {
                    "V1": 370.0,
                    "V2": 2500.0,
                    "DV": 5.0,
                    "FWHM": 6.0,
                    "FILTNM": "/raid/ATMCOR/M6AC/data/filters/aviris_ng.flt"
                },
                "FILEOPTIONS": {
                    "NOPRNT": 2,
                    "FLROOT": "/path/to/modtran_output/basename"
                }
            }
        }
    ]
}
```

### Critical Geometry Parameters

**IMPORTANT**: These values prevent MODTRAN geometry warnings:

- `H2ALT: 0.30000` - Ground altitude (km) - **MUST be less than H1ALT**
- `LENN: 0` - Earth radius model (automatically set by MODTRAN when H2ALT < H1ALT)
- `H1ALT: 4.830596380703965` - Sensor altitude (km)

If H2ALT >= H1ALT, MODTRAN will issue warnings and auto-correct values.

### Template Location

Default template for ang20190624t230039:
```
/raid/ATMCOR/M6AC/data/ang20190624t230039/modtran_template.json
```

## Running the Code

### Command Syntax

```bash
export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:/home/judy/MODTRAN6.0/bin/linux:$LD_LIBRARY_PATH

/opt/mpich-3.3.2/bin/mpirun -np <NUM_PROCESSES> ./lut_mpi_hdf5 \
    <TEMPLATE_JSON> \
    <OUTPUT_H5> \
    <ALBEDO_MIN> <ALBEDO_MAX> <N_ALBEDO> \
    <H2O_MIN> <H2O_MAX> <N_H2O>
```

### Example: Test Run (2×3 grid)

```bash
cd /home/judy/MODTRAN6.0/developer/examples

export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:/home/judy/MODTRAN6.0/bin/linux:$LD_LIBRARY_PATH

time /opt/mpich-3.3.2/bin/mpirun -np 4 ./lut_mpi_hdf5 \
    /raid/ATMCOR/M6AC/data/ang20190624t230039/modtran_template.json \
    /raid/ATMCOR/M6AC/tests/MOD6_VAL/test_lut_mpi_final.h5 \
    0.0 1.0 2 \
    0.0 4.0 3
```

**Expected output:**
```
========================================
MPI + HDF5 LUT Generation
========================================

Template:  /raid/ATMCOR/M6AC/data/ang20190624t230039/modtran_template.json
Output:    /raid/ATMCOR/M6AC/tests/MOD6_VAL/test_lut_mpi_final.h5
Processes: 4 (3 workers)

Parameter grid: 2 albedos × 3 H2O = 6 runs
...
Task 0 completed by Worker 1 (1/6 = 16.7%)
...
SUCCESS!

real    0m6.140s
user    0m23.269s
sys     0m1.116s
```

### Example: Production Run (11×151 grid)

```bash
cd /home/judy/MODTRAN6.0/developer/examples

export LD_LIBRARY_PATH=/opt/mpich-3.3.2/lib:/home/judy/MODTRAN6.0/bin/linux:$LD_LIBRARY_PATH

time /opt/mpich-3.3.2/bin/mpirun -np 16 ./lut_mpi_hdf5 \
    /raid/ATMCOR/M6AC/data/ang20190624t230039/modtran_template.json \
    /raid/ATMCOR/M6AC/data/ang20190624t230039/lut/ang20190624t230039_modtran_lut.h5 \
    0.0 1.0 11 \
    0.0 4.0 151
```

**Estimated runtime:** ~1.9 minutes (1661 runs with 16 processes)

## Output Format

### HDF5 Structure

The output HDF5 file contains:

**1D Arrays (Parameter Grids):**
- `albedo_grid` - [n_albedo] - Surface albedo values
- `h2o_grid` - [n_h2o] - Water vapor column values (g/cm²)

**3D Arrays (Spectral Data):** Shape: [n_albedo, n_h2o, n_channels]
- `total_radiance` - Total at-sensor radiance (W/(cm²·sr·nm))
- `transmission` - Total atmospheric transmission
- `channel_radiance` - Channel-specific radiance
- `grnd_rflt` - Ground-reflected radiance
- `mult_scat` - Multiple scattering component
- `sing_scat` - Single scattering component

### Wavelength Grid

AVIRIS-NG: 425 channels from 370-2500 nm
```python
wavelengths = np.linspace(370, 2500, 425)
```

### Accessing Data with Python

```python
import h5py
import numpy as np

# Open file
f = h5py.File('output.h5', 'r')

# Read parameter grids
albedo_grid = f['albedo_grid'][:]  # Shape: (n_albedo,)
h2o_grid = f['h2o_grid'][:]        # Shape: (n_h2o,)

# Read spectral data
total_radiance = f['total_radiance'][:]  # Shape: (n_albedo, n_h2o, 425)
transmission = f['transmission'][:]      # Shape: (n_albedo, n_h2o, 425)

# Access specific case
i_albedo = 5  # Index
j_h2o = 100   # Index
spectrum = total_radiance[i_albedo, j_h2o, :]  # Shape: (425,)

f.close()
```

## Validation

### Plot Full Spectrum

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File('test_lut_mpi_final.h5', 'r')

albedo_grid = f['albedo_grid'][:]
h2o_grid = f['h2o_grid'][:]
total_radiance = f['total_radiance'][:]

wavelengths = np.linspace(370, 2500, 425)

fig, axes = plt.subplots(len(albedo_grid), len(h2o_grid), figsize=(15, 10))
for i, alb in enumerate(albedo_grid):
    for j, h2o in enumerate(h2o_grid):
        ax = axes[i, j]
        rad = total_radiance[i, j, :]
        ax.plot(wavelengths, rad * 1e6)  # Convert to µW/(cm²·sr·nm)
        ax.set_title(f'Albedo={alb:.1f}, H2O={h2o:.1f} g/cm²')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Radiance (µW/(cm²·sr·nm))')
        ax.set_xlim(370, 2500)  # Full AVIRIS-NG range
        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('lut_validation.png', dpi=150)
f.close()
```

### Expected Radiance Range

- **Minimum:** ~2.3×10⁻¹¹ W/(cm²·sr·nm) (deep absorption bands)
- **Maximum:** ~2.3×10⁻⁵ W/(cm²·sr·nm) (peak solar spectrum)
- **Mean:** ~5.4×10⁻⁶ W/(cm²·sr·nm)

### Spectral Features to Verify

- **UV/Blue peak:** 370-500 nm (~20-22 µW/(cm²·sr·nm))
- **Water vapor absorption:** 940, 1130, 1380, 1880 nm
- **Oxygen absorption:** 760 nm (A-band)
- **Albedo effect:** Higher albedo → higher radiance across all wavelengths
- **H2O effect:** Stronger absorption with higher water vapor

## Architecture Details

### Why MPI Instead of OpenMP?

**Problem:** MODTRAN's FORTRAN core uses module-level allocatable arrays that are not thread-safe.

**Attempted with OpenMP:**
```
forrtl: severe (151): allocatable array is already allocated
```

**Solution:** MPI provides separate process memory spaces, avoiding FORTRAN thread conflicts.

### Master-Worker Pattern

- **Master (rank 0):** Distributes tasks, collects results, writes HDF5
- **Workers (rank 1-N):** Execute MODTRAN runs, send spectral data to master
- **Dynamic scheduling:** Workers request new tasks upon completion (load balancing)

### Per-Process JSON Modification

Each worker:
1. Receives task parameters (albedo, H2O)
2. Creates local copy of JSON template
3. Updates `SURREF` (albedo) and `H2OSTR` (water vapor)
4. Sets `NOFILE: "FC_NOFILES"` to disable file I/O conflicts
5. Executes MODTRAN via C++ API
6. Extracts channel radiance data
7. Sends results to master via MPI

## Troubleshooting

### Problem: "cannot find -lmod6rt"

**Solution:** Verify MODTRAN library path
```bash
ls -l /home/judy/MODTRAN6.0/bin/linux/libmod6rt.so
export LD_LIBRARY_PATH=/home/judy/MODTRAN6.0/bin/linux:$LD_LIBRARY_PATH
```

### Problem: "H2ALT has been reset" warnings

**Solution:** Update template JSON with corrected geometry
```json
"GEOMETRY": {
    "H2ALT": 0.30000,
    "LENN": 0
}
```
Ensure H2ALT < H1ALT (ground altitude less than sensor altitude).

### Problem: Execution hangs

**Causes:**
1. File I/O conflicts between parallel processes
2. Missing `NOFILE: "FC_NOFILES"` in worker JSON configs

**Solution:** Code automatically sets this flag in `lut_mpi_hdf5.cc:213`

### Problem: MPI version conflicts

**Symptoms:**
```
symbol lookup error: /usr/lib/x86_64-linux-gnu/libmpi.so.40
```

**Solution:** Use isolated MPICH installation
```bash
/opt/mpich-3.3.2/bin/mpirun -np 4 ./lut_mpi_hdf5 ...
```

## File Locations

### Source Code
- `lut_mpi_hdf5.cc`: Main MPI+HDF5 LUT generator
- Location: `/home/judy/MODTRAN6.0/developer/examples/`

### Configuration
- Template: `/raid/ATMCOR/M6AC/data/ang20190624t230039/modtran_template.json`
- Filter: `/raid/ATMCOR/M6AC/data/filters/aviris_ng.flt`

### Output
- Test HDF5: `/raid/ATMCOR/M6AC/tests/MOD6_VAL/test_lut_mpi_final.h5`
- Production: `/raid/ATMCOR/M6AC/data/ang20190624t230039/lut/`

### Libraries
- MPICH: `/opt/mpich-3.3.2/`
- MODTRAN: `/home/judy/MODTRAN6.0/bin/linux/`

## References

- **MODTRAN 6 Documentation:** Available in `/home/judy/MODTRAN6.0/documentation/`
- **nlohmann/json:** https://github.com/nlohmann/json
- **HDF5:** https://www.hdfgroup.org/solutions/hdf5/
- **MPICH:** https://www.mpich.org/

## Version History

- **2025-11-09:** Initial production system
  - MPI-based parallel execution
  - HDF5 output format
  - Validated with AVIRIS-NG ang20190624t230039
  - Performance: ~1 sec/run (4 processes)

## Contact

For questions about this implementation, refer to:
- MODTRAN Support: modtran@spectral.com
- AFRL: jeannette.van_den_bosch@us.af.mil
