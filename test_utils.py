#!/usr/bin/env python3
"""Test script for m6ac.utils module"""
import numpy as np

# Direct file import to avoid __init__.py issues
import importlib.util
spec = importlib.util.spec_from_file_location("utils", "/raid/ATMCOR/M6AC/m6ac/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)

# Test with our ang20190624t230039 scene
scene_path = '/raid/AVIRIS_NG/imagery/ang20190624t230039_rdn_v2u1/ang20190624t230039_rdn_v2u1'

print("="*80)
print("Testing M6AC utils.py")
print("="*80)

# Test 1: Parse filename
print("\n1. Parsing filename...")
metadata = utils.parse_aviris_filename('ang20190624t230039')
for key, val in metadata.items():
    print(f"   {key}: {val}")

# Test 2: Load radiance
print("\n2. Loading radiance...")
radiance, wavelengths = utils.load_aviris_ng(scene_path + '_img')
print(f"   Radiance shape: {radiance.shape}")
print(f"   Wavelengths: {len(wavelengths)} bands, {wavelengths[0]:.2f}-{wavelengths[-1]:.2f} nm")
print(f"   Radiance range: {float(radiance.min()):.2f} - {float(radiance.max()):.2f}")

# Test 3: Load observation geometry
print("\n3. Loading observation geometry...")
geometry = utils.load_obs_geometry(scene_path + '_obs_ort')
print(f"   Solar zenith (center): {float(geometry['solar_zenith'][527, 318]):.2f}°")
print(f"   Solar azimuth (center): {float(geometry['solar_azimuth'][527, 318]):.2f}°")
print(f"   Solar elevation (center): {float(geometry['solar_elevation'][527, 318]):.2f}°")

# Test 4: Load location
print("\n4. Loading location...")
location = utils.load_location(scene_path + '_loc_ort')
print(f"   Latitude (center): {float(location['latitude'][527, 318]):.4f}°")
print(f"   Longitude (center): {float(location['longitude'][527, 318]):.4f}°")
print(f"   Elevation (center): {float(location['elevation'][527, 318]):.1f} m")

# Test 5: Get scene center parameters
print("\n5. Computing scene-center parameters for MODTRAN...")
params = utils.get_scene_center_params(geometry, location)
for key, val in params.items():
    print(f"   {key}: {val:.4f}")

print("\n" + "="*80)
print("All tests passed!")
print("="*80)
