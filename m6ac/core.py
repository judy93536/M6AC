"""
Core M6AC atmospheric correction class
"""
import numpy as np
import spectral.io.envi as envi
from pathlib import Path
from typing import Tuple, Optional, Union

from .retrieval import water_vapor_cibr, aerosol_ddv
from .lut import LUTInterpolator
from .utils import load_aviris_ng, save_reflectance, mask_bad_bands


class M6AC:
    """
    MODTRAN6-based Atmospheric Correction for AVIRIS-NG

    Parameters
    ----------
    modtran_path : str or Path
        Path to MODTRAN6 installation directory
    sensor_filter : str or Path
        Path to AVIRIS-NG filter file (.flt)
    lut_path : str or Path, optional
        Path to pre-computed LUT file
    """

    def __init__(
        self,
        modtran_path: Union[str, Path],
        sensor_filter: Union[str, Path],
        lut_path: Optional[Union[str, Path]] = None
    ):
        self.modtran_path = Path(modtran_path)
        self.sensor_filter = Path(sensor_filter)
        self.lut = None

        if lut_path is not None:
            self.load_lut(lut_path)

        # Bad bands for AVIRIS-NG (water absorption)
        # 1350-1450nm, 1800-1950nm, 2450-2600nm
        self.bad_band_ranges = [
            (1350, 1450),
            (1800, 1950),
            (2450, 2600)
        ]

    def load_lut(self, lut_path: Union[str, Path]):
        """Load pre-computed MODTRAN6 LUT"""
        self.lut = LUTInterpolator(lut_path)

    def load_aviris_ng(
        self,
        radiance_path: Union[str, Path],
        border: int = 30
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load AVIRIS-NG radiance image

        Parameters
        ----------
        radiance_path : str or Path
            Path to AVIRIS-NG radiance image (without extension)
        border : int, optional
            Border pixels to remove (default: 30)

        Returns
        -------
        radiance : ndarray
            Radiance image (ny, nx, nbands)
        wavelengths : ndarray
            Wavelength vector (nbands,)
        """
        return load_aviris_ng(radiance_path, border)

    def retrieve_parameters(
        self,
        radiance: np.ndarray,
        wavelengths: np.ndarray,
        solar_zenith: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Retrieve atmospheric parameters from radiance

        Parameters
        ----------
        radiance : ndarray
            At-sensor radiance (ny, nx, nbands)
        wavelengths : ndarray
            Wavelength vector (nbands,)
        solar_zenith : float
            Solar zenith angle (degrees)

        Returns
        -------
        h2o_map : ndarray
            Water vapor column density map (ny, nx) in g/cm²
        vis_map : ndarray
            Visibility map (ny, nx) in km
        """
        # Convert radiance to apparent reflectance for retrieval
        mu_s = np.cos(np.deg2rad(solar_zenith))
        apparent_rfl = np.pi * radiance / (self.lut.solar_irradiance[None, None, :] * mu_s)

        # CIBR water vapor retrieval
        h2o_map = water_vapor_cibr(apparent_rfl, wavelengths)

        # DDV aerosol/visibility retrieval
        vis_map = aerosol_ddv(apparent_rfl, wavelengths, h2o_map)

        return h2o_map, vis_map

    def correct(
        self,
        radiance: np.ndarray,
        h2o_map: np.ndarray,
        vis_map: np.ndarray,
        wavelengths: Optional[np.ndarray] = None,
        mask_bad: bool = True
    ) -> np.ndarray:
        """
        Apply atmospheric correction

        Parameters
        ----------
        radiance : ndarray
            At-sensor radiance (ny, nx, nbands)
        h2o_map : ndarray
            Water vapor map (ny, nx)
        vis_map : ndarray
            Visibility map (ny, nx)
        wavelengths : ndarray, optional
            Wavelength vector for bad band masking
        mask_bad : bool, optional
            Whether to mask bad bands (default: True)

        Returns
        -------
        reflectance : ndarray
            Surface reflectance (ny, nx, nbands)
        """
        if self.lut is None:
            raise ValueError("LUT not loaded. Call load_lut() first.")

        # Interpolate LUT for each pixel
        reflectance = self.lut.interpolate(radiance, h2o_map, vis_map)

        # Mask bad bands
        if mask_bad and wavelengths is not None:
            reflectance = mask_bad_bands(
                reflectance,
                wavelengths,
                self.bad_band_ranges
            )

        return reflectance

    def process_scene(
        self,
        radiance_path: Union[str, Path],
        output_path: Union[str, Path],
        solar_zenith: float,
        border: int = 30
    ) -> dict:
        """
        Complete atmospheric correction workflow

        Parameters
        ----------
        radiance_path : str or Path
            Path to input radiance image
        output_path : str or Path
            Path for output reflectance image
        solar_zenith : float
            Solar zenith angle (degrees)
        border : int, optional
            Border pixels to remove (default: 30)

        Returns
        -------
        results : dict
            Dictionary containing reflectance, h2o_map, vis_map
        """
        # Load radiance
        print("Loading radiance...")
        radiance, wavelengths = self.load_aviris_ng(radiance_path, border)

        # Retrieve atmospheric parameters
        print("Retrieving atmospheric parameters...")
        h2o_map, vis_map = self.retrieve_parameters(
            radiance, wavelengths, solar_zenith
        )

        # Apply correction
        print("Applying atmospheric correction...")
        reflectance = self.correct(
            radiance, h2o_map, vis_map, wavelengths
        )

        # Save output
        print(f"Saving reflectance to {output_path}...")
        save_reflectance(
            reflectance,
            output_path,
            radiance_path,
            wavelengths
        )

        return {
            'reflectance': reflectance,
            'h2o_map': h2o_map,
            'vis_map': vis_map,
            'wavelengths': wavelengths
        }
