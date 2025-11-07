# M6AC - MODTRAN6-based Atmospheric Correction for AVIRIS-NG

M6AC is a physics-based atmospheric correction method for AVIRIS-NG hyperspectral imagery, combining proven retrieval algorithms with MODTRAN6's native spectral response function handling.

## Features

- **CIBR Water Vapor Retrieval**: Continuum Interpolated Band Ratio method for per-pixel water column estimation
- **FLAASH-style Aerosol Estimation**: Dark dense vegetation (DDV) method for visibility/aerosol optical depth
- **MODTRAN6 Integration**: Native SRF convolution eliminates manual band alignment
- **LUT-based Correction**: Fast operational processing via pre-computed lookup tables
- **Bad Band Masking**: Proper handling of atmospheric absorption features

## Background

M6AC builds on the Modified MODTRAN4 Atmospheric Correction (M4AC) method published in:

> Northrop, J., et al. (2021). "AVIRIS-NG Reflectance..." [Your paper citation]

The key innovation of M4AC was manual convolution of MODTRAN4 outputs with AVIRIS-NG spectral response functions. MODTRAN6's native filter file support (`FILTNM` parameter) now handles this automatically, simplifying the workflow while maintaining scientific accuracy.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/M6AC.git
cd M6AC

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install M6AC in development mode
pip install -e .
```

## Requirements

- Python 3.8+
- MODTRAN6.0+ (licensed separately from Spectral Sciences Inc.)
- AVIRIS-NG spectral response function file
- NumPy, SciPy, Spectral Python, matplotlib

## Quick Start

```python
from m6ac import M6AC
from m6ac.retrieval import water_vapor_cibr, aerosol_ddv

# Initialize atmospheric correction
atcor = M6AC(
    modtran_path='/home/judy/MODTRAN6.0',
    sensor_filter='/raid/AVIRIS_NG/data/tools/aviris_ng.flt'
)

# Load AVIRIS-NG radiance
radiance, wavelengths = atcor.load_aviris_ng(
    '/path/to/ang_rdn_img'
)

# Retrieve atmospheric parameters
h2o_map = water_vapor_cibr(radiance, wavelengths)
vis_map = aerosol_ddv(radiance, wavelengths, h2o_map)

# Apply atmospheric correction using pre-computed LUT
reflectance = atcor.correct(radiance, h2o_map, vis_map)

# Save output
atcor.save_reflectance(reflectance, 'output_rfl')
```

## LUT Generation

Generate scene-specific MODTRAN6 lookup tables:

```python
from m6ac.lut import generate_lut

# Create LUT for specific scene geometry
lut = generate_lut(
    output_path='./luts/scene_lut.npy',
    solar_zenith=41.73,
    sensor_zenith=8.92,
    relative_azimuth=100.54,
    elevation=0.3,  # km
    h2o_range=(0.0, 4.0, 40),  # min, max, steps (g/cm²)
    vis_range=(5.0, 100.0, 20),  # min, max, steps (km)
    albedo_range=(0.0, 1.0, 11)  # min, max, steps
)
```

## Workflow

1. **Parameter Retrieval**: Estimate atmospheric water vapor and aerosol optical depth from radiance
2. **LUT Generation**: Pre-compute MODTRAN6 radiative transfer for parameter space
3. **Atmospheric Correction**: Interpolate LUT and apply correction pixel-by-pixel
4. **Quality Control**: Mask bad bands and apply uncertainty estimates

## Project Structure

```
M6AC/
├── m6ac/                  # Main package
│   ├── __init__.py
│   ├── core.py           # M6AC main class
│   ├── retrieval.py      # Water vapor & aerosol retrieval
│   ├── lut.py            # LUT generation & interpolation
│   ├── modtran.py        # MODTRAN6 interface
│   └── utils.py          # Utilities & I/O
├── examples/             # Example scripts
├── tests/                # Unit tests
├── docs/                 # Documentation
├── setup.py              # Package setup
├── requirements.txt      # Dependencies
└── README.md            # This file
```

## Comparison to ISOFIT

| Feature | M6AC | ISOFIT |
|---------|------|--------|
| Physics basis | MODTRAN6 LUTs | sRTMnet emulator |
| Speed | Fast (LUT lookup) | Slow (13+ hours full OE) |
| Water band masking | ✅ Correct | ❌ Artifacts |
| Operational | ✅ Production-ready | ⚠️ Research tool |
| Dependencies | MODTRAN6 license | None |

## License

[Your chosen license - e.g., MIT, BSD, GPL]

## Citation

If you use M6AC in your research, please cite:

```bibtex
@article{northrop2021aviris,
  title={AVIRIS-NG Reflectance...},
  author={Northrop, J. and ...},
  journal={...},
  year={2021}
}
```

## Contact

[Your contact information]

## Acknowledgments

- MODTRAN6 by Spectral Sciences Inc.
- AVIRIS-NG data from NASA JPL
- CIBR method from Green et al. (1989)
- FLAASH DDV method from Kaufman et al.
