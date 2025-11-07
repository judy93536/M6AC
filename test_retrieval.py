#!/usr/bin/env python3
"""
Test script for m6ac/retrieval.py

Tests CIBR water vapor and DDV aerosol retrieval algorithms
using the Caltech AVIRIS-NG scene (ang20190624t230039)
"""
import numpy as np
import sys
import importlib.util

# Direct import to avoid __init__.py issues
spec_utils = importlib.util.spec_from_file_location("utils", "/raid/ATMCOR/M6AC/m6ac/utils.py")
utils = importlib.util.module_from_spec(spec_utils)
spec_utils.loader.exec_module(utils)

spec_retrieval = importlib.util.spec_from_file_location("retrieval", "/raid/ATMCOR/M6AC/m6ac/retrieval.py")
retrieval = importlib.util.module_from_spec(spec_retrieval)
spec_retrieval.loader.exec_module(retrieval)

print("="*70)
print("M6AC Retrieval Algorithm Test")
print("="*70)

# Load AVIRIS-NG scene
scene_path = "/raid/AVIRIS_NG/imagery/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1_img"
obs_path = "/raid/AVIRIS_NG/imagery/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1_obs_ort"
loc_path = "/raid/AVIRIS_NG/imagery/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1_loc_ort"

print("\n1. Loading AVIRIS-NG data...")
radiance, wavelengths = utils.load_aviris_ng(scene_path, border=30)
geometry = utils.load_obs_geometry(obs_path, border=30)
location = utils.load_location(loc_path, border=30)

print(f"   Radiance shape: {radiance.shape}")
print(f"   Wavelengths: {len(wavelengths)} bands, {wavelengths[0]:.2f}-{wavelengths[-1]:.2f} nm")

# Get scene parameters
scene_params = utils.get_scene_center_params(geometry, location)
solar_zenith = scene_params['solar_zenith']
print(f"   Solar zenith: {solar_zenith:.2f}°")

# Convert radiance to apparent reflectance for retrieval algorithms
# Using simple conversion - will use proper MODTRAN solar irradiance later
PI = np.pi
mu_s = np.cos(np.deg2rad(solar_zenith))

# Simple TOA solar irradiance approximation (W/m²/sr/nm)
# Scaled to give reasonable reflectance values
solar_irrad_approx = 1500.0 / (wavelengths / 550.0)**2

# Convert radiance (µW/cm²/nm/sr) to apparent reflectance
# radiance * pi / (solar_irrad * mu_s)
# Note: 1 µW/cm² = 0.01 W/m², so multiply by 0.01
apparent_rfl = np.zeros_like(radiance, dtype=np.float32)
for i in range(len(wavelengths)):
    apparent_rfl[:, :, i] = radiance[:, :, i] * 0.01 * PI / (solar_irrad_approx[i] * mu_s)

# Clip to reasonable range
apparent_rfl = np.clip(apparent_rfl, 0, 2.0)

print(f"   Apparent reflectance range: {np.nanmin(apparent_rfl):.4f} - {np.nanmax(apparent_rfl):.4f}")
print(f"   Apparent reflectance mean: {np.nanmean(apparent_rfl):.4f}")

print(f"\n2. Testing CIBR water vapor retrieval (original coefficients)...")
h2o_ratio = retrieval.water_vapor_cibr(apparent_rfl, wavelengths, method='1.14um', use_refined_coeffs=False)
print(f"   H2O ratio shape: {h2o_ratio.shape}")
print(f"   H2O ratio range: {np.nanmin(h2o_ratio):.4f} - {np.nanmax(h2o_ratio):.4f}")
print(f"   H2O ratio mean: {np.nanmean(h2o_ratio):.4f}")
print(f"   (Lower ratio = more water vapor absorption)")

print(f"\n3. Testing CIBR water vapor retrieval (refined coefficients)...")
h2o_ratio_refined = retrieval.water_vapor_cibr(apparent_rfl, wavelengths, method='1.14um', use_refined_coeffs=True)
print(f"   H2O ratio range: {np.nanmin(h2o_ratio_refined):.4f} - {np.nanmax(h2o_ratio_refined):.4f}")
print(f"   H2O ratio mean: {np.nanmean(h2o_ratio_refined):.4f}")
print(f"   Difference from original: {np.nanmean(h2o_ratio_refined - h2o_ratio):.6f}")

print(f"\n4. Testing NDVI calculation...")
ndvi = retrieval.compute_ndvi(apparent_rfl, wavelengths)
print(f"   NDVI shape: {ndvi.shape}")
print(f"   NDVI range: {np.nanmin(ndvi):.4f} - {np.nanmax(ndvi):.4f}")
print(f"   NDVI mean: {np.nanmean(ndvi):.4f}")

# Count vegetation pixels at different thresholds
veg_08 = np.sum(ndvi > 0.8)
veg_09 = np.sum(ndvi > 0.9)
veg_092 = np.sum(ndvi > 0.92)
total_pixels = ndvi.size
print(f"   Vegetation pixels (NDVI > 0.8): {veg_08} ({100.0*veg_08/total_pixels:.1f}%)")
print(f"   Vegetation pixels (NDVI > 0.9): {veg_09} ({100.0*veg_09/total_pixels:.1f}%)")
print(f"   Vegetation pixels (NDVI > 0.92): {veg_092} ({100.0*veg_092/total_pixels:.1f}%)")

print(f"\n5. Testing DDV aerosol/visibility retrieval...")
try:
    visibility_map, ddv_mask = retrieval.aerosol_ddv(apparent_rfl, wavelengths,
                                                      ndvi_threshold=0.8,
                                                      ratio_range=(1.99, 2.01),
                                                      dilation_iterations=6)

    print(f"   Visibility map shape: {visibility_map.shape}")
    print(f"   DDV mask shape: {ddv_mask.shape}")

    ddv_pixels = np.sum(ddv_mask)
    ddv_percent = 100.0 * ddv_pixels / ddv_mask.size
    print(f"   DDV pixels found: {ddv_pixels} ({ddv_percent:.2f}%)")

    if ddv_pixels > 0:
        valid_vis = visibility_map[ddv_mask]
        print(f"   Visibility range: {np.min(valid_vis):.2f} - {np.max(valid_vis):.2f} km")
        print(f"   Visibility mean: {np.mean(valid_vis):.2f} km")
        print(f"   Visibility median: {np.median(valid_vis):.2f} km")
    else:
        print(f"   WARNING: No DDV pixels found!")

except Exception as e:
    print(f"   ERROR: {e}")
    import traceback
    traceback.print_exc()

print(f"\n6. Testing scene visibility estimation...")
try:
    scene_vis = retrieval.estimate_scene_visibility(apparent_rfl, wavelengths, default_visibility=35.0)
    print(f"   Scene visibility: {scene_vis:.2f} km")
    print(f"\n   NOTE: M4AC result for this scene was ~29.6 km")

except Exception as e:
    print(f"   ERROR: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("Test completed!")
print("="*70)
