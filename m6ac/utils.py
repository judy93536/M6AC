"""
Utility functions for loading AVIRIS-NG data and metadata extraction
"""
import numpy as np
import spectral.io.envi as envi
from pathlib import Path
from typing import Tuple, Optional, Union, Dict
from datetime import datetime


def parse_aviris_filename(filename: str) -> Dict[str, any]:
    """
    Parse AVIRIS-NG filename to extract date/time metadata

    Filename format: angYYYYMMDDtHHMMSS
    Example: ang20190624t230039

    Parameters
    ----------
    filename : str
        AVIRIS-NG scene identifier

    Returns
    -------
    metadata : dict
        Dictionary containing year, month, day, hour, minute, second, day_of_year
    """
    # Remove path and extension if present
    scene_id = Path(filename).stem

    # Extract date and time components
    # Format: angYYYYMMDDtHHMMSS
    if not scene_id.startswith('ang'):
        raise ValueError(f"Invalid AVIRIS-NG filename: {filename}")

    date_time = scene_id[3:]  # Remove 'ang' prefix
    year = int(date_time[0:4])
    month = int(date_time[4:6])
    day = int(date_time[6:8])
    hour = int(date_time[9:11])
    minute = int(date_time[11:13])
    second = int(date_time[13:15])

    # Calculate day of year
    dt = datetime(year, month, day)
    day_of_year = dt.timetuple().tm_yday

    # UTC time in decimal hours
    utc_time = hour + minute / 60.0 + second / 3600.0

    return {
        'year': year,
        'month': month,
        'day': day,
        'hour': hour,
        'minute': minute,
        'second': second,
        'day_of_year': day_of_year,
        'utc_time': utc_time
    }


def load_aviris_ng(
    radiance_path: Union[str, Path],
    border: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load AVIRIS-NG radiance image

    Parameters
    ----------
    radiance_path : str or Path
        Path to AVIRIS-NG radiance image (without extension)
    border : int, optional
        Border pixels to remove from each edge (default: 0)

    Returns
    -------
    radiance : ndarray
        Radiance image (ny, nx, nbands) in µW/(cm² nm sr)
    wavelengths : ndarray
        Wavelength vector (nbands,) in nm
    """
    radiance_path = str(radiance_path)

    # Load ENVI file
    img = envi.open(radiance_path + '.hdr', radiance_path)
    radiance = img.load()

    # Extract wavelengths
    wavelengths = np.array([float(w) for w in img.metadata['wavelength']])

    # Remove border pixels if requested
    if border > 0:
        radiance = radiance[border:-border, border:-border, :]

    return radiance, wavelengths


def load_obs_geometry(
    obs_path: Union[str, Path],
    border: int = 0
) -> Dict[str, np.ndarray]:
    """
    Load AVIRIS-NG observation geometry file (_obs_ort)

    Parameters
    ----------
    obs_path : str or Path
        Path to observation geometry file (without extension)
    border : int, optional
        Border pixels to remove (default: 0)

    Returns
    -------
    geometry : dict
        Dictionary containing:
        - path_length: Sensor-to-ground distance (m)
        - solar_azimuth: Solar azimuth angle (deg, 0-360 CW from N)
        - solar_zenith: Solar zenith angle (deg, 0-90 from zenith)
        - solar_elevation: Solar elevation angle (deg, 90-zenith)
        - utc_time: UTC time (decimal hours)
        - earth_sun_dist: Earth-sun distance (AU)
    """
    obs_path = str(obs_path)

    # Load observation file
    obs = envi.open(obs_path + '.hdr', obs_path)
    obs_data = obs.load()

    if border > 0:
        obs_data = obs_data[border:-border, border:-border, :]

    # Extract bands (0-indexed)
    # Band 1: Path length (m)
    # Band 4: To-sun azimuth (0-360° CW from N)
    # Band 5: To-sun zenith (0-90° from zenith)
    # Band 10: UTC Time (decimal hours)
    # Band 11: Earth-sun distance (AU)

    geometry = {
        'path_length': obs_data[:, :, 0],
        'solar_azimuth': obs_data[:, :, 3],
        'solar_zenith': obs_data[:, :, 4],
        'solar_elevation': 90.0 - obs_data[:, :, 4],
        'utc_time': obs_data[:, :, 9],
        'earth_sun_dist': obs_data[:, :, 10]
    }

    return geometry


def load_location(
    loc_path: Union[str, Path],
    border: int = 0
) -> Dict[str, np.ndarray]:
    """
    Load AVIRIS-NG location file (_loc_ort)

    Parameters
    ----------
    loc_path : str or Path
        Path to location file (without extension)
    border : int, optional
        Border pixels to remove (default: 0)

    Returns
    -------
    location : dict
        Dictionary containing:
        - longitude: Longitude in WGS-84 (deg)
        - latitude: Latitude in WGS-84 (deg)
        - elevation: Elevation (m)
    """
    loc_path = str(loc_path)

    # Load location file
    loc = envi.open(loc_path + '.hdr', loc_path)
    loc_data = loc.load()

    if border > 0:
        loc_data = loc_data[border:-border, border:-border, :]

    # Extract bands (0-indexed)
    # Band 1: Longitude (WGS-84 deg)
    # Band 2: Latitude (WGS-84 deg)
    # Band 3: Elevation (m)

    location = {
        'longitude': loc_data[:, :, 0],
        'latitude': loc_data[:, :, 1],
        'elevation': loc_data[:, :, 2]
    }

    return location


def get_scene_center_params(
    geometry: Dict[str, np.ndarray],
    location: Dict[str, np.ndarray]
) -> Dict[str, float]:
    """
    Extract scene-center parameters for MODTRAN configuration

    Parameters
    ----------
    geometry : dict
        Observation geometry dictionary from load_obs_geometry()
    location : dict
        Location dictionary from load_location()

    Returns
    -------
    params : dict
        Scene-center parameters:
        - sensor_altitude: Sensor altitude (km)
        - ground_altitude: Ground elevation (km)
        - latitude: Center latitude (deg)
        - longitude: Center longitude (deg)
        - solar_zenith: Center solar zenith (deg)
        - solar_azimuth: Center solar azimuth (deg)
    """
    # Get shape - handle both 2D and 3D arrays from spectral library
    shape = geometry['path_length'].shape
    if len(shape) == 3:
        ny, nx, _ = shape
    else:
        ny, nx = shape

    # Use center region (middle 10% of image)
    cy_start, cy_end = int(ny * 0.45), int(ny * 0.55)
    cx_start, cx_end = int(nx * 0.45), int(nx * 0.55)

    center_slice = (slice(cy_start, cy_end), slice(cx_start, cx_end))

    # Calculate sensor altitude from path length + ground elevation
    path_km = np.nanmean(geometry['path_length'][center_slice]) / 1000.0
    elev_km = np.nanmean(location['elevation'][center_slice]) / 1000.0
    sensor_alt_km = path_km + elev_km

    params = {
        'sensor_altitude': sensor_alt_km,
        'ground_altitude': elev_km,
        'latitude': np.nanmean(location['latitude'][center_slice]),
        'longitude': np.nanmean(location['longitude'][center_slice]),
        'solar_zenith': np.nanmean(geometry['solar_zenith'][center_slice]),
        'solar_azimuth': np.nanmean(geometry['solar_azimuth'][center_slice]),
    }

    return params


def mask_bad_bands(
    reflectance: np.ndarray,
    wavelengths: np.ndarray,
    bad_ranges: list = [(1350, 1450), (1800, 1950), (2450, 2600)]
) -> np.ndarray:
    """
    Mask bad bands (water absorption regions) by setting to zero

    Parameters
    ----------
    reflectance : ndarray
        Reflectance image (ny, nx, nbands)
    wavelengths : ndarray
        Wavelength vector (nbands,)
    bad_ranges : list of tuples, optional
        List of (min_nm, max_nm) ranges to mask
        Default: [(1350, 1450), (1800, 1950), (2450, 2600)]

    Returns
    -------
    reflectance_masked : ndarray
        Reflectance with bad bands set to zero
    """
    reflectance_masked = reflectance.copy()

    for wl_min, wl_max in bad_ranges:
        mask = (wavelengths >= wl_min) & (wavelengths <= wl_max)
        reflectance_masked[:, :, mask] = 0.0

    return reflectance_masked


def save_reflectance(
    reflectance: np.ndarray,
    output_path: Union[str, Path],
    template_path: Union[str, Path],
    wavelengths: np.ndarray
):
    """
    Save reflectance image in ENVI format

    Parameters
    ----------
    reflectance : ndarray
        Reflectance image (ny, nx, nbands)
    output_path : str or Path
        Output file path (without extension)
    template_path : str or Path
        Template radiance file for copying metadata
    wavelengths : ndarray
        Wavelength vector (nbands,)
    """
    output_path = str(output_path)
    template_path = str(template_path)

    # Load template metadata
    template = envi.open(template_path + '.hdr', template_path)
    metadata = template.metadata.copy()

    # Update metadata for reflectance
    metadata['description'] = 'M6AC Surface Reflectance'
    metadata['data type'] = '4'  # 32-bit float
    metadata['wavelength'] = [str(w) for w in wavelengths]
    metadata['data ignore value'] = '0'

    # Save reflectance
    envi.save_image(
        output_path + '.hdr',
        reflectance,
        metadata=metadata,
        force=True,
        interleave='bil'
    )


def extract_scene_metadata(
    scene_id: str,
    base_path: Optional[Union[str, Path]] = None,
    border: int = 30
) -> Dict[str, any]:
    """
    Extract all collection parameters from AVIRIS-NG files in one call

    This convenience function combines filename parsing, observation geometry,
    location data, and scene-center parameter extraction into a single call.

    Parameters
    ----------
    scene_id : str
        AVIRIS-NG scene identifier (e.g., 'ang20190624t230039_rdn_v2u1')
    base_path : str or Path, optional
        Base path to scene directory. If None, assumes standard structure:
        /raid/AVIRIS_NG/imagery/{scene_id}/
    border : int, optional
        Border pixels to remove when computing scene-center parameters (default: 30)

    Returns
    -------
    metadata : dict
        Complete metadata dictionary containing:

        **From filename:**
        - scene_id: str - Full scene identifier
        - year, month, day, hour, minute, second: int
        - day_of_year: int (1-366)
        - utc_time: float (decimal hours)

        **From observation geometry (_obs_ort):**
        - solar_zenith: float (degrees, scene center)
        - solar_azimuth: float (degrees, scene center)
        - earth_sun_dist: float (AU)

        **From location (_loc_ort):**
        - latitude: float (degrees, scene center)
        - longitude: float (degrees, scene center)
        - ground_altitude: float (km, scene center)

        **Computed:**
        - sensor_altitude: float (km, scene center)

        **File paths:**
        - radiance_path: str - Path to radiance image (without extension)
        - obs_path: str - Path to observation geometry (without extension)
        - loc_path: str - Path to location file (without extension)

    Examples
    --------
    >>> # Extract metadata for a scene
    >>> metadata = extract_scene_metadata('ang20190624t230039_rdn_v2u1')
    >>> print(f"Solar zenith: {metadata['solar_zenith']:.2f} degrees")
    >>> print(f"Sensor altitude: {metadata['sensor_altitude']:.3f} km")

    >>> # Use custom base path
    >>> metadata = extract_scene_metadata(
    ...     'ang20190624t230039_rdn_v2u1',
    ...     base_path='/custom/path/to/scene'
    ... )

    >>> # Use directly with MODTRAN LUT generation
    >>> from m6ac.lut import LUTGenerator
    >>> generator = LUTGenerator(sensor_filter_path='aviris_ng.flt')
    >>> metadata = extract_scene_metadata('ang20190624t230039_rdn_v2u1')
    >>> lut = generator.generate_lut(
    ...     sensor_altitude_km=metadata['sensor_altitude'],
    ...     ground_altitude_km=metadata['ground_altitude'],
    ...     solar_zenith_deg=metadata['solar_zenith'],
    ...     latitude=metadata['latitude'],
    ...     longitude=metadata['longitude'],
    ...     day_of_year=metadata['day_of_year'],
    ...     utc_time=metadata['utc_time'],
    ...     h2o_column_gcm2=2.0,
    ...     visibility_km=35.0
    ... )
    """
    # Remove any file extensions from scene_id
    scene_id = Path(scene_id).stem

    # Determine base path
    if base_path is None:
        base_path = Path('/raid/AVIRIS_NG/imagery') / scene_id
    else:
        base_path = Path(base_path)

    # Construct file paths
    radiance_path = base_path / f'{scene_id}_img'
    obs_path = base_path / f'{scene_id}_obs_ort'
    loc_path = base_path / f'{scene_id}_loc_ort'

    # 1. Parse filename for date/time
    filename_params = parse_aviris_filename(scene_id)

    # 2. Load observation geometry
    geometry = load_obs_geometry(obs_path, border=border)

    # 3. Load location data
    location = load_location(loc_path, border=border)

    # 4. Extract scene-center parameters
    scene_params = get_scene_center_params(geometry, location)

    # 5. Combine all parameters into unified dictionary
    # Convert numpy types to native Python types for JSON serialization
    metadata = {
        # Scene identifier
        'scene_id': scene_id,

        # Date and time (from filename)
        'year': filename_params['year'],
        'month': filename_params['month'],
        'day': filename_params['day'],
        'hour': filename_params['hour'],
        'minute': filename_params['minute'],
        'second': filename_params['second'],
        'day_of_year': filename_params['day_of_year'],
        'utc_time': filename_params['utc_time'],

        # Solar geometry (scene center)
        'solar_zenith': float(scene_params['solar_zenith']),
        'solar_azimuth': float(scene_params['solar_azimuth']),

        # Location (scene center)
        'latitude': float(scene_params['latitude']),
        'longitude': float(scene_params['longitude']),

        # Altitude (scene center)
        'sensor_altitude': float(scene_params['sensor_altitude']),
        'ground_altitude': float(scene_params['ground_altitude']),

        # Earth-Sun distance (scene average)
        'earth_sun_dist': float(geometry['earth_sun_dist'].mean()),

        # File paths for reference
        'radiance_path': str(radiance_path),
        'obs_path': str(obs_path),
        'loc_path': str(loc_path),
        'base_path': str(base_path),
    }

    return metadata
