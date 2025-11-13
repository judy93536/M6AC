# M6AC Workflow Step 1: Scene Metadata Extraction - COMPLETE ✓

## Summary

Successfully implemented and tested a one-call convenience function for extracting all AVIRIS-NG collection parameters needed for atmospheric correction.

---

## What Was Implemented

### New Function: `extract_scene_metadata()`

**Location:** `/raid/ATMCOR/M6AC/m6ac/utils.py` (lines 328-471)

**Purpose:** Extract all collection parameters from AVIRIS-NG `*_obs_ort` and `*_loc_ort` files in a single function call.

**Signature:**
```python
def extract_scene_metadata(
    scene_id: str,
    base_path: Optional[Union[str, Path]] = None,
    border: int = 30
) -> Dict[str, any]
```

**Usage:**
```python
from m6ac.utils import extract_scene_metadata

# Simple usage - assumes standard path structure
metadata = extract_scene_metadata('ang20190624t230039_rdn_v2u1')

# Custom base path
metadata = extract_scene_metadata(
    'ang20190624t230039_rdn_v2u1',
    base_path='/custom/path/to/scene'
)
```

---

## Extracted Parameters

The function returns a complete metadata dictionary with:

### Date & Time (from filename)
- `scene_id`: Scene identifier
- `year`, `month`, `day`, `hour`, `minute`, `second`
- `day_of_year`: 1-366
- `utc_time`: Decimal hours

### Solar Geometry (from *_obs_ort)
- `solar_zenith`: Solar zenith angle (degrees, scene center)
- `solar_azimuth`: Solar azimuth angle (degrees, scene center)
- `earth_sun_dist`: Earth-Sun distance (AU)

### Location (from *_loc_ort)
- `latitude`: Scene center latitude (degrees)
- `longitude`: Scene center longitude (degrees)
- `ground_altitude`: Ground elevation (km, scene center)

### Sensor
- `sensor_altitude`: Sensor altitude (km, scene center)

### File Paths (for reference)
- `radiance_path`: Path to radiance image
- `obs_path`: Path to observation geometry
- `loc_path`: Path to location file
- `base_path`: Scene base directory

---

## Test Results

**Scene:** ang20190624t230039_rdn_v2u1

```
COLLECTION PARAMETERS:
  Date/Time: 2019-06-24 23:00:39 UTC
  Day of Year: 175
  UTC Time: 23.0108 hours

  Solar Zenith: 41.6730°
  Solar Azimuth: 267.6376°
  Earth-Sun Distance: 1.016419 AU

  Sensor Altitude: 5.1544 km
  Ground Altitude: 0.3208 km
  Latitude: 33.851173°
  Longitude: -116.977646°
```

**Status:** ✅ All parameters extracted successfully

---

## Files Created/Modified

### Modified
1. `/raid/ATMCOR/M6AC/m6ac/utils.py`
   - Added `extract_scene_metadata()` function
   - Converts numpy types to Python native types for JSON compatibility

2. `/raid/ATMCOR/M6AC/m6ac/__init__.py`
   - Exported `extract_scene_metadata` for easy import

### Created
1. `/raid/ATMCOR/M6AC/examples/workflow_step1_extract_metadata.py`
   - Example script demonstrating Step 1 of workflow
   - Shows how to use extracted metadata for LUT generation

2. `/raid/ATMCOR/M6AC/data/ang20190624t230039_rdn_v2u1/scene_metadata.json`
   - JSON file with extracted parameters
   - Ready for downstream processing

---

## Integration with Workflow

### Current Workflow (Original Code)
```python
# Multiple manual steps
hdrFile = '/raid/AVIRIS/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1_img.hdr'
imgFile = '/raid/AVIRIS/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1_img'
L, wave = ht.getAVNG(imgFile, (30,30))
# ... manual parameter specification
```

### M6AC Workflow (Productionized)
```python
# Single function call extracts everything
from m6ac.utils import extract_scene_metadata

metadata = extract_scene_metadata('ang20190624t230039_rdn_v2u1')

# Parameters ready for LUT generation
print(f"Solar zenith: {metadata['solar_zenith']:.4f}°")
print(f"Sensor altitude: {metadata['sensor_altitude']:.4f} km")
```

---

## Next Workflow Step: LUT Generation

With metadata extracted, the next step is to generate MODTRAN LUTs:

```python
from m6ac.lut import LUTGenerator

# Initialize generator
generator = LUTGenerator(
    sensor_filter_path='/raid/AVIRIS_NG/data/tools/aviris_ng.flt'
)

# Generate LUT grid using extracted metadata
luts = generator.generate_lut_grid(
    sensor_altitude_km=metadata['sensor_altitude'],
    ground_altitude_km=metadata['ground_altitude'],
    solar_zenith_deg=metadata['solar_zenith'],
    latitude=metadata['latitude'],
    longitude=metadata['longitude'],
    day_of_year=metadata['day_of_year'],
    utc_time=metadata['utc_time'],
    h2o_values=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],  # g/cm²
    vis_values=[5, 10, 20, 35, 50, 80, 120],              # km
    surface_albedo=0.3
)
# → Generates 8×7 = 56 HDF5 LUT files
```

---

## Benefits

1. **Single function call** - No need to manually call 4 separate functions
2. **Automatic path handling** - Assumes standard AVIRIS-NG directory structure
3. **JSON-compatible** - Returns native Python types for easy serialization
4. **Scene-center averaging** - Uses middle 10% of image for representative parameters
5. **Border handling** - Excludes edge pixels when computing scene-center values
6. **Complete documentation** - Comprehensive docstring with examples

---

## Testing

**Test script:** `examples/workflow_step1_extract_metadata.py`

**Run:**
```bash
cd /raid/ATMCOR/M6AC
python3 examples/workflow_step1_extract_metadata.py
```

**Expected output:**
- ✓ Metadata extraction from .hdr and binary files
- ✓ Display of all collection parameters
- ✓ JSON file saved to `data/{scene_id}/scene_metadata.json`
- ✓ Next steps guidance

---

## Code Quality

- ✅ Comprehensive docstring with examples
- ✅ Type hints for all parameters
- ✅ Handles both 2D and 3D numpy arrays
- ✅ Converts numpy types to Python types
- ✅ Exported in `__init__.py` for easy import
- ✅ Tested with real AVIRIS-NG data
- ✅ JSON serialization verified

---

## Production Ready

This function is **production-ready** and can be used as the first step in processing all 20+ AVIRIS-NG scenes.

**Workflow automation:**
```python
scene_ids = [
    'ang20190624t230039_rdn_v2u1',
    'ang20190625t193845_rdn_v2u1',
    # ... add all scenes
]

for scene_id in scene_ids:
    metadata = extract_scene_metadata(scene_id)
    # Save metadata for each scene
    # Proceed to LUT generation...
```

---

## Date Completed

**2025-11-13**

**Next:** Workflow Step 2 - LUT Grid Generation
