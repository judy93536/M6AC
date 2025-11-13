#!/usr/bin/env python3
"""
M6AC Workflow - Step 1: Extract Scene Metadata

This script demonstrates the first step in the M6AC atmospheric correction
workflow: extracting all collection parameters from AVIRIS-NG files.

Author: M6AC Development Team
Date: 2025-11-13
"""
import sys
sys.path.insert(0, '/raid/ATMCOR/M6AC/m6ac')

from utils import extract_scene_metadata
import json
from pathlib import Path


def main():
    # Scene to process
    scene_id = 'ang20190624t230039_rdn_v2u1'

    print("=" * 70)
    print("M6AC WORKFLOW - STEP 1: EXTRACT SCENE METADATA")
    print("=" * 70)
    print(f"\nScene ID: {scene_id}\n")

    # Extract all collection parameters in one call
    print("Extracting metadata from AVIRIS-NG files...")
    metadata = extract_scene_metadata(scene_id)
    print("✓ Metadata extraction complete\n")

    # Display extracted parameters
    print("COLLECTION PARAMETERS:")
    print(f"  Date/Time: {metadata['year']}-{metadata['month']:02d}-{metadata['day']:02d} "
          f"{metadata['hour']:02d}:{metadata['minute']:02d}:{metadata['second']:02d} UTC")
    print(f"  Day of Year: {metadata['day_of_year']}")
    print(f"  UTC Time: {metadata['utc_time']:.4f} hours")
    print()
    print(f"  Solar Zenith: {metadata['solar_zenith']:.4f}°")
    print(f"  Solar Azimuth: {metadata['solar_azimuth']:.4f}°")
    print(f"  Earth-Sun Distance: {metadata['earth_sun_dist']:.6f} AU")
    print()
    print(f"  Sensor Altitude: {metadata['sensor_altitude']:.4f} km")
    print(f"  Ground Altitude: {metadata['ground_altitude']:.4f} km")
    print(f"  Latitude: {metadata['latitude']:.6f}°")
    print(f"  Longitude: {metadata['longitude']:.6f}°")

    # Save to JSON
    output_dir = Path(f'/raid/ATMCOR/M6AC/data/{scene_id}')
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / 'scene_metadata.json'
    with open(output_file, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"\n✓ Metadata saved to: {output_file}")

    print("\n" + "=" * 70)
    print("NEXT STEPS:")
    print("=" * 70)
    print("1. Generate LUT grid using these parameters")
    print("2. Load radiance image and create validity mask")
    print("3. Retrieve atmospheric parameters (H2O, visibility)")
    print("4. Apply atmospheric correction")
    print()
    print("Example:")
    print("  from m6ac.lut import LUTGenerator")
    print("  generator = LUTGenerator(")
    print("      sensor_filter_path='/raid/AVIRIS_NG/data/tools/aviris_ng.flt'")
    print("  )")
    print("  luts = generator.generate_lut_grid(")
    print(f"      sensor_altitude_km={metadata['sensor_altitude']:.4f},")
    print(f"      ground_altitude_km={metadata['ground_altitude']:.4f},")
    print(f"      solar_zenith_deg={metadata['solar_zenith']:.4f},")
    print(f"      latitude={metadata['latitude']:.4f},")
    print(f"      longitude={metadata['longitude']:.4f},")
    print(f"      day_of_year={metadata['day_of_year']},")
    print(f"      utc_time={metadata['utc_time']:.4f},")
    print("      h2o_values=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],")
    print("      vis_values=[5, 10, 20, 35, 50, 80, 120],")
    print("      surface_albedo=0.3")
    print("  )")
    print()


if __name__ == '__main__':
    main()
