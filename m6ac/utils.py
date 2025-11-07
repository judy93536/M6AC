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
