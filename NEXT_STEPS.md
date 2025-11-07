# M6AC Development - Next Steps

## What's Been Created

✅ Project structure at `/raid/ATMCOR/M6AC/`
✅ Python virtual environment (`venv/`)
✅ Git repository initialized
✅ Package structure (`m6ac/` directory)
✅ README.md with documentation
✅ setup.py for installation
✅ .gitignore for Python projects

## To Use This Project

### 1. Activate the Virtual Environment

```bash
cd /raid/ATMCOR/M6AC
source venv/bin/activate  # Your prompt will show (venv)
```

### 2. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
pip install -e .  # Install M6AC in development mode
```

### 3. Development Workflow

The conda environment you're currently in won't interfere. When you want to work on M6AC:
- Activate the venv: `source venv/bin/activate`
- Make changes to code
- Test: `python -m pytest tests/`
- Deactivate when done: `deactivate`

## Modules To Implement

Priority order for porting your existing M4AC code:

### 1. `m6ac/utils.py` (NEXT - utilities needed by everything)
- `load_aviris_ng()` - Load ENVI format
- `save_reflectance()` - Write output
- `mask_bad_bands()` - Zero out water absorption bands

### 2. `m6ac/retrieval.py` (HIGH PRIORITY - your proven algorithms)
Port from `/raid/Python_MODTRAN/MODTRAN/pyWork/spectralUtils.py`:
- `water_vapor_cibr()` - CIBR method (0.94µm and 1.14µm bands)
- `aerosol_ddv()` - Dark dense vegetation visibility estimation
- Helper functions for band selection

### 3. `m6ac/modtran.py` (HIGH PRIORITY - MODTRAN6 interface)
Port from `/raid/Python_MODTRAN/MODTRAN6/local-packages/RT_MOD6.py`:
- `ModtranRunner` class - Execute MODTRAN6 from Python
- `parse_chn()` - Read channel radiance files
- `parse_tp7()` - Read tape7 transmission files
- JSON template management

### 4. `m6ac/lut.py` (MEDIUM PRIORITY - LUT generation & interpolation)
Port from your existing LUT generation scripts:
- `LUTGenerator` - Create MODTRAN6 LUTs
  - Generate parameter grid (H2O, visibility, albedo)
  - Create JSON input for each LUT point
  - Run MODTRAN6 in batch
  - Parse outputs into numpy array
- `LUTInterpolator` - Fast pixel-by-pixel correction
  - Load pre-computed LUT
  - Multi-dimensional interpolation (scipy.interpolate)

### 5. `m6ac/cli.py` (LOW PRIORITY - command-line tools)
- `m6ac-correct` - Process a scene
- `m6ac-lut` - Generate LUTs

### 6. `examples/` (MEDIUM PRIORITY - usage examples)
- `process_caltech_scene.py` - Complete workflow
- `generate_scene_lut.py` - LUT generation
- `compare_to_m4ac.py` - Validation

### 7. `tests/` (MEDIUM PRIORITY - unit tests)
- `test_retrieval.py` - CIBR, DDV tests
- `test_lut.py` - Interpolation tests
- `test_utils.py` - I/O tests

## Key Differences from M4AC

| Task | M4AC (MODTRAN4) | M6AC (MODTRAN6) |
|------|-----------------|-----------------|
| SRF convolution | Manual (your code) | Native `FILTNM` parameter |
| Input | TAPE5 text files | JSON |
| Execution | Batch scripts | Python API |
| Output parsing | Parse TAPE7 | Parse .chn/.tp7 |

## Testing Against M4AC Results

Your M4AC output for ang20190624t230039:
- Location: `/raid/AVIRIS_NG/imagery/ang20190624t230039_rdn_v2u1/M4AC/`
- File: `ang20190624t230039_rfl_v2u1_img`
- Stats: Mean rfl = 0.1265, NDVI = 0.428

Use this as validation when M6AC is working.

## Git Workflow

```bash
# Check status
git status

# Add changes
git add m6ac/retrieval.py

# Commit
git commit -m "Implement CIBR water vapor retrieval"

# Create GitHub repo (on github.com), then:
git remote add origin https://github.com/yourusername/M6AC.git
git branch -M main  # Rename master to main (optional)
git push -u origin main
```

## Resources

- Your M4AC code: `/raid/Python_MODTRAN/MODTRAN/pyWork/`
- MODTRAN6 work: `/raid/Python_MODTRAN/MODTRAN6/`
- IDL code: `/raid/AVIRIS_NG/idl/work/`
- AVIRIS-NG filter: `/raid/AVIRIS_NG/data/tools/aviris_ng.flt`
- MODTRAN6 install: `/home/judy/MODTRAN6.0/`

## Questions for Implementation

1. Do your existing MODTRAN6 LUTs need regeneration, or can we use them?
2. Preferred H2O/vis LUT grid resolution?
3. Which MODTRAN6 Python API do you prefer (developer/ vs. custom)?
4. License for GitHub (MIT, BSD, GPL)?
