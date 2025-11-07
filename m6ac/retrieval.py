"""
Atmospheric parameter retrieval algorithms (CIBR and DDV methods)

Based on proven M4AC algorithms from Northrop (2021)
"""
import numpy as np
from scipy import ndimage
from typing import Tuple, Optional


def find_band_indices(wavelengths: np.ndarray, wl_min: float, wl_max: float) -> np.ndarray:
    """
    Find indices of bands within wavelength range

    Parameters
    ----------
    wavelengths : ndarray
        Wavelength array (nm)
    wl_min : float
        Minimum wavelength (nm)
    wl_max : float
        Maximum wavelength (nm)

    Returns
    -------
    indices : ndarray
        Band indices within range
    """
    return np.where((wavelengths >= wl_min) & (wavelengths <= wl_max))[0]


def water_vapor_cibr(
    radiance: np.ndarray,
    wavelengths: np.ndarray,
    method: str = '1.14um',
    use_refined_coeffs: bool = False
) -> np.ndarray:
    """
    Water vapor retrieval using Continuum Interpolated Band Ratio (CIBR)

    Based on Green et al. (1989) method used in FLAASH and M4AC.

    Parameters
    ----------
    radiance : ndarray
        At-sensor radiance or apparent reflectance (ny, nx, nbands)
    wavelengths : ndarray
        Wavelength vector (nbands,) in nm
    method : str, optional
        Which absorption feature to use:
        - '1.14um': 1140 nm feature (default, more sensitive)
        - '0.94um': 940 nm feature (alternative)
    use_refined_coeffs : bool, optional
        Use refined coefficients from M4AC Caltech scene calibration.
        Original Green et al. (1989): a=0.571, b=0.429
        Refined (Caltech 2019): a=0.615388, b=0.384612
        Default: False (use original coefficients)
        Note: Refined coeffs only tested on ang20190624t230039

    Returns
    -------
    h2o_ratio : ndarray
        Water vapor absorption ratio (ny, nx)
        Higher values indicate more water vapor

    References
    ----------
    Green, R.O., Carrere, V., Conel, J.E., 1989. Measurement of Atmospheric
    Water Vapor Using the Airborne Visible/Infrared Imaging Spectrometer.
    Image Processing '89, Sparks, Nevada, 23 May 1989.

    Notes
    -----
    To convert ratio to water vapor column (g/cm²), use calibration curve
    from MODTRAN LUT (see M4AC paper Fig. 3).
    """
    if method == '1.14um':
        # 1.14 µm method (higher sensitivity)
        # Absorption channel: ~1130 nm
        idx_abs = find_band_indices(wavelengths, 1116, 1153)

        # Reference channels
        idx_ref1 = find_band_indices(wavelengths, 1057, 1077)  # 1.065 µm
        idx_ref2 = find_band_indices(wavelengths, 1211, 1231)  # 1.225 µm

        # Weighting coefficients for continuum interpolation
        if use_refined_coeffs:
            # Refined from M4AC IDL code (flaash_map.pro) for Caltech 2019 scene
            # NOTE: Only tested on ang20190624t230039 - may not generalize
            a = 0.615388
            b = 0.384612
        else:
            # Original Green et al. (1989) - more general
            a = 0.571
            b = 0.429

    elif method == '0.94um':
        # 0.94 µm method (alternative)
        # Absorption channel: ~940 nm
        idx_abs = find_band_indices(wavelengths, 923, 962)

        # Reference channels
        idx_ref1 = find_band_indices(wavelengths, 855, 876)   # 0.865 µm
        idx_ref2 = find_band_indices(wavelengths, 1019, 1039) # 1.025 µm

        # Weighting coefficients
        a = 0.539
        b = 0.461
    else:
        raise ValueError(f"Unknown method: {method}. Use '1.14um' or '0.94um'")

    # Average radiance in each spectral window
    L_abs = radiance[:, :, idx_abs].mean(axis=2)
    L_ref1 = radiance[:, :, idx_ref1].mean(axis=2)
    L_ref2 = radiance[:, :, idx_ref2].mean(axis=2)

    # Continuum interpolation
    L_continuum = a * L_ref1 + b * L_ref2

    # Band ratio (CIBR)
    # Avoid division by zero
    h2o_ratio = np.divide(L_abs, L_continuum,
                          out=np.ones_like(L_abs),
                          where=L_continuum > 0)

    return h2o_ratio


def compute_ndvi(
    reflectance: np.ndarray,
    wavelengths: np.ndarray
) -> np.ndarray:
    """
    Compute Normalized Difference Vegetation Index (NDVI)

    Parameters
    ----------
    reflectance : ndarray
        Surface reflectance (ny, nx, nbands)
    wavelengths : ndarray
        Wavelength vector (nbands,) in nm

    Returns
    -------
    ndvi : ndarray
        NDVI values (ny, nx), range [-1, 1]
    """
    # Red: ~660 nm, NIR: ~850 nm
    idx_red = find_band_indices(wavelengths, 655, 665)
    idx_nir = find_band_indices(wavelengths, 845, 855)

    if len(idx_red) == 0 or len(idx_nir) == 0:
        raise ValueError("Could not find red or NIR bands for NDVI calculation")

    red = reflectance[:, :, idx_red[len(idx_red)//2]]
    nir = reflectance[:, :, idx_nir[len(idx_nir)//2]]

    # NDVI formula
    ndvi = np.divide(nir - red, nir + red,
                     out=np.zeros_like(nir),
                     where=(nir + red) > 0)

    return ndvi


def aerosol_ddv(
    reflectance: np.ndarray,
    wavelengths: np.ndarray,
    ndvi_threshold: float = 0.8,
    ratio_range: Tuple[float, float] = (1.99, 2.01),
    dilation_iterations: int = 6
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Aerosol/visibility estimation using Dark Dense Vegetation (DDV) method

    Based on Kaufman et al. (1997) and used in FLAASH/M4AC.

    Parameters
    ----------
    reflectance : ndarray
        Apparent surface reflectance (ny, nx, nbands)
    wavelengths : ndarray
        Wavelength vector (nbands,) in nm
    ndvi_threshold : float, optional
        NDVI threshold for vegetation (default: 0.8)
    ratio_range : tuple of float, optional
        Acceptable range for ρ(2100nm)/ρ(660nm) ratio (default: 1.99-2.01)
    dilation_iterations : int, optional
        Morphological dilation iterations (default: 6)

    Returns
    -------
    visibility_map : ndarray
        Visibility estimates in km (ny, nx)
        Zero where DDV pixels not found
    ddv_mask : ndarray
        Boolean mask of DDV pixels (ny, nx)

    References
    ----------
    Kaufman, Y.J., Wald, A.E., Remer, L.A., Gao, B-C., Li, R-R., Flynn, L.,
    1997. The MODIS 2.1-μm channel-correlation with visible reflectance
    for use in remote sensing of aerosol. IEEE Trans. Geosci. Remote
    Sensing 35(5), 1286-1298.

    Notes
    -----
    The DDV method assumes:
    - Dense vegetation has ρ(2100nm)/ρ(660nm) ≈ 2.0
    - Dark pixels at 550 nm correlate with aerosol optical depth
    - Visibility ≈ ln(50) / (ρ(550nm) + 0.01159)
    """
    ny, nx = reflectance.shape[:2]

    # Find band indices
    idx_660 = find_band_indices(wavelengths, 655, 665)
    idx_2100 = find_band_indices(wavelengths, 2090, 2110)
    idx_550 = find_band_indices(wavelengths, 545, 555)

    if len(idx_660) == 0 or len(idx_2100) == 0 or len(idx_550) == 0:
        raise ValueError("Could not find required bands for DDV method")

    # Use center wavelength of each band
    rho_660 = reflectance[:, :, idx_660[len(idx_660)//2]]
    rho_2100 = reflectance[:, :, idx_2100[len(idx_2100)//2]]
    rho_550 = reflectance[:, :, idx_550[len(idx_550)//2]]

    # Compute NDVI
    ndvi = compute_ndvi(reflectance, wavelengths)

    # Create vegetation mask (NDVI > threshold)
    veg_mask = ndvi > ndvi_threshold

    # Compute spectral ratio
    ratio_2100_660 = np.divide(rho_2100, rho_660,
                                out=np.zeros_like(rho_2100),
                                where=rho_660 > 0)

    # Create ratio mask (should be ~2.0 for vegetation)
    ratio_mask = (ratio_2100_660 >= ratio_range[0]) & (ratio_2100_660 <= ratio_range[1])

    # Apply morphological dilation to expand regions
    veg_mask_dilated = ndimage.binary_dilation(veg_mask, iterations=dilation_iterations)
    ratio_mask_dilated = ndimage.binary_dilation(ratio_mask, iterations=dilation_iterations)

    # DDV pixels are intersection of both masks
    ddv_mask = veg_mask_dilated & ratio_mask_dilated

    # Calculate visibility from dark pixels
    # Formula from M4AC: vis = ln(50) / (ρ(550nm) + 0.01159)
    visibility_map = np.zeros((ny, nx), dtype=np.float32)
    visibility_map[ddv_mask] = np.log(50.0) / (rho_550[ddv_mask] + 0.01159)

    return visibility_map, ddv_mask


def estimate_scene_visibility(
    reflectance: np.ndarray,
    wavelengths: np.ndarray,
    default_visibility: float = 35.0
) -> float:
    """
    Estimate mean scene visibility using DDV method

    Parameters
    ----------
    reflectance : ndarray
        Apparent surface reflectance (ny, nx, nbands)
    wavelengths : ndarray
        Wavelength vector (nbands,) in nm
    default_visibility : float, optional
        Default visibility if DDV pixels not found (default: 35.0 km)

    Returns
    -------
    visibility : float
        Mean scene visibility in km
    """
    try:
        visibility_map, ddv_mask = aerosol_ddv(reflectance, wavelengths)

        # Get valid DDV pixels
        ddv_pixels = visibility_map[ddv_mask]

        if len(ddv_pixels) == 0:
            print(f"Warning: No DDV pixels found. Using default visibility: {default_visibility} km")
            return default_visibility

        # Return mean visibility from DDV pixels
        mean_visibility = np.mean(ddv_pixels)

        # Sanity check (reasonable visibility range: 5-100 km)
        if mean_visibility < 5.0 or mean_visibility > 100.0:
            print(f"Warning: DDV visibility {mean_visibility:.1f} km outside reasonable range.")
            print(f"Using default visibility: {default_visibility} km")
            return default_visibility

        return float(mean_visibility)

    except Exception as e:
        print(f"Warning: DDV method failed ({e}). Using default visibility: {default_visibility} km")
        return default_visibility
