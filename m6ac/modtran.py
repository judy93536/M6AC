"""
MODTRAN6 interface using pymodtran API

Provides Python wrapper for running MODTRAN6 and extracting results
for atmospheric correction LUT generation.
"""
import sys
import re
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

# Add MODTRAN6 Python API to path
MODTRAN6_PYTHON_PATH = "/home/judy/MODTRAN6.0/developer/python"
if MODTRAN6_PYTHON_PATH not in sys.path:
    sys.path.insert(0, MODTRAN6_PYTHON_PATH)

import pymodtran as pymod


def wgs84_to_modtran_longitude(wgs84_lon: float) -> float:
    """
    Convert WGS-84 longitude to MODTRAN6 convention

    WGS-84 convention: -180 to +180 degrees
        - Negative values = West of Greenwich
        - Positive values = East of Greenwich

    MODTRAN6 convention: 0 to 360 degrees measured WEST from Greenwich
        - 0° = Greenwich Meridian
        - 90° = 90° West
        - 180° = International Date Line (going west)
        - 270° = 90° East (360° - 90°)
        - 360° = Greenwich (wraps around)

    Parameters
    ----------
    wgs84_lon : float
        Longitude in WGS-84 convention (-180 to +180)

    Returns
    -------
    modtran_lon : float
        Longitude in MODTRAN6 convention (0 to 360, degrees west)

    Examples
    --------
    >>> # California (western hemisphere)
    >>> wgs84_to_modtran_longitude(-116.98)
    116.98

    >>> # Egypt (eastern hemisphere)
    >>> wgs84_to_modtran_longitude(30.0)
    330.0

    >>> # Prime Meridian
    >>> wgs84_to_modtran_longitude(0.0)
    0.0
    """
    if wgs84_lon < 0:
        # Western hemisphere: degrees west of Greenwich
        # Simply take absolute value
        return abs(wgs84_lon)
    else:
        # Eastern hemisphere: convert to degrees west
        # 30° East = 330° West (360 - 30)
        return 360.0 - wgs84_lon


class ModtranRunner:
    """
    MODTRAN6 runner using pymodtran API

    Configures and executes MODTRAN6 for atmospheric correction LUT generation

    Parameters
    ----------
    sensor_filter_path : str or Path, optional
        Path to sensor spectral response filter (.flt file)
        For AVIRIS-NG: /raid/AVIRIS_NG/data/tools/aviris_ng.flt
    """

    def __init__(self, sensor_filter_path: Optional[Union[str, Path]] = None):
        self.sensor_filter_path = str(sensor_filter_path) if sensor_filter_path else None
        self.modlib = None
        self.case_name = None
        self.h2o_column_gcm2 = None
        self.solar_zenith_deg = None

    def configure_case(
        self,
        sensor_altitude_km: float,
        ground_altitude_km: float,
        solar_zenith_deg: float,
        solar_azimuth_deg: Optional[float] = None,
        latitude: float = 34.0,
        longitude: float = -118.0,
        day_of_year: int = 175,
        utc_time: float = 23.0,
        h2o_column_gcm2: float = 2.0,
        visibility_km: float = 35.0,
        co2_ppmv: float = 420.0,
        surface_albedo: float = 0.0,
        wl_min_nm: float = 350.0,
        wl_max_nm: float = 2550.0,
        dv_cmm1: float = 5.0,
        fwhm_cmm1: float = 10.0,
        atmosphere_model: int = 6,  # 6 = Mid-Latitude Summer
        case_name: str = "M6AC"
    ) -> int:
        """
        Configure MODTRAN6 case for atmospheric correction

        Parameters
        ----------
        sensor_altitude_km : float
            Sensor altitude (km above sea level)
        ground_altitude_km : float
            Ground elevation (km above sea level)
        solar_zenith_deg : float
            Solar zenith angle (degrees, 0-90)
        solar_azimuth_deg : float, optional
            Solar azimuth angle (degrees, 0-360 CW from N)
        latitude : float
            Scene latitude (degrees, WGS-84)
        longitude : float
            Scene longitude (degrees, WGS-84: -180 to +180, negative=west)
            Automatically converted to MODTRAN6 convention (0-360° west)
        day_of_year : int
            Day of year (1-365)
        utc_time : float
            UTC time (decimal hours, 0-24)
        h2o_column_gcm2 : float
            Water vapor column density (g/cm²)
        visibility_km : float
            Visibility (km)
        co2_ppmv : float
            CO2 concentration (ppmv)
        surface_albedo : float
            Lambertian surface albedo (0-1)
        wl_min_nm : float
            Minimum wavelength (nm)
        wl_max_nm : float
            Maximum wavelength (nm)
        dv_cmm1 : float
            Internal spectral resolution (cm⁻¹)
        fwhm_cmm1 : float
            Internal FWHM (cm⁻¹, typically 2×DV)
        atmosphere_model : int
            MODTRAN atmosphere model:
            1=Tropical, 2=Mid-Lat Summer, 3=Mid-Lat Winter,
            4=Sub-Arctic Summer, 5=Sub-Arctic Winter, 6=US Standard
        case_name : str
            Case identifier

        Returns
        -------
        case_idx : int
            Case index (0-based)
        """
        # Initialize MODTRAN library if needed
        if self.modlib is None:
            self.modlib = pymod.ModLib()

        # Create case
        self.modlib.caseCreate(1)
        case_idx = 0
        modin = self.modlib.caseInput(case_idx)

        # Save case parameters for .chn file parsing
        self.case_name = case_name
        self.h2o_column_gcm2 = h2o_column_gcm2
        self.solar_zenith_deg = solar_zenith_deg

        # Basic configuration
        modin.name = case_name
        modin.options.iemsct = pymod.RT_SOLAR_AND_THERMAL

        # Multiple scattering configuration
        # IMULT: Multiple scattering algorithm
        #   0 = RT_NO_MULTIPLE_SCATTER (single scatter only)
        #   1 = RT_DISORT (discrete ordinate, high fidelity)
        #   3 = RT_ISAACS_2STREAM (fast two-stream)
        #   5 = RT_ISAACS_SCALED (DISORT at key points, Isaac's elsewhere)
        modin.options.imult = 1  # Use DISORT for high fidelity

        # Spectral configuration
        # Convert wavelength (nm) to wavenumber (cm⁻¹): ν = 10^7 / λ
        v1 = 1e7 / wl_max_nm  # Min wavenumber (max wavelength)
        v2 = 1e7 / wl_min_nm  # Max wavenumber (min wavelength)

        modin.spectral.v1 = v1
        modin.spectral.v2 = v2
        modin.spectral.dv = dv_cmm1
        modin.spectral.fwhm = fwhm_cmm1

        # Sensor spectral response filter (FILTNM)
        if self.sensor_filter_path:
            modin.spectral.filtnm = self.sensor_filter_path

        # Geometry configuration
        # H1ALT: Sensor altitude (km)
        # H2ALT: Ground altitude (km)
        # IDAY: Day of year
        # GMTIME: GMT time (decimal hours)
        # PARM1: Latitude (degrees N)
        # PARM2: Longitude (degrees west of Greenwich, 0-360)
        # PARM3/PARM4: Second location (same as H2 location for our case)

        # Convert WGS-84 longitude to MODTRAN6 convention
        modtran_longitude = wgs84_to_modtran_longitude(longitude)

        modin.geometry.h1alt = sensor_altitude_km
        modin.geometry.h2alt = ground_altitude_km
        modin.geometry.iday = day_of_year
        modin.geometry.gmtime = utc_time
        modin.geometry.parm1 = latitude
        modin.geometry.parm2 = modtran_longitude  # Degrees west of Greenwich (0-360)
        modin.geometry.parm3 = latitude           # Same location at ground
        modin.geometry.parm4 = modtran_longitude

        # Observer zenith angle (OBSZEN)
        # For downward-looking (H1 > H2): OBSZEN = 180°
        # For upward-looking (H1 < H2): OBSZEN < 90°
        modin.geometry.obszen = 180.0

        # Atmosphere model
        modin.atmos.model = atmosphere_model

        # Water vapor column
        # Use MODTRAN's water vapor scaling
        # h2o_column in g/cm² → scale factor relative to model
        modin.atmos.h2ostr = h2o_column_gcm2  # Direct H2O column
        modin.atmos.h2ounit = 'g'  # g/cm²

        # Visibility (aerosol)
        modin.aerosols.vis = visibility_km

        # CO2 concentration
        modin.atmos.co2mx = co2_ppmv

        # Surface properties
        # SURFTYPE: Surface reflectance type
        #   REFL_CONSTANT (0): Spectrally constant Lambertian reflectance
        #   REFL_LAMBER_MODEL (1): Lambertian spectral reflectance from file
        #   REFL_BRDF (3): BRDF
        # SURREF: Lambertian surface reflectance (albedo)
        modin.surface.surftype = 0  # REFL_CONSTANT
        modin.surface.surref = surface_albedo

        # File output options
        # Enable .chn and tape6 output for inspection
        modin.fileopt.nofile = 0  # Write all output files including .chn

        return case_idx

    def run(self) -> bool:
        """
        Execute MODTRAN6 run

        Returns
        -------
        success : bool
            True if execution succeeded
        """
        if self.modlib is None:
            raise RuntimeError("No MODTRAN case configured. Call configure_case() first.")

        self.modlib.execute()
        status_msg = self.modlib.statusMessage()

        # Check for errors
        if "error" in status_msg.lower() or "fatal" in status_msg.lower():
            print(f"MODTRAN6 execution error: {status_msg}")
            return False

        return True

    def read_chn_file(self, chn_file_path: Optional[Union[str, Path]] = None) -> Dict[str, np.ndarray]:
        """
        Read and parse MODTRAN6 .chn file for sensor-convolved data

        Based on load_chn() from mod6Utils.py. The .chn file contains
        spectral data convolved with the sensor spectral response function
        (FILTNM filter file), resulting in sensor-band-matched outputs.

        Parameters
        ----------
        chn_file_path : str or Path, optional
            Path to .chn file. If not specified, looks for {case_name}.chn
            in current working directory.

        Returns
        -------
        results : dict
            Dictionary containing sensor-convolved data:
            - wavelength: Channel center wavelengths (nm)
            - channel_width: Channel widths (nm)
            - solar_irradiance: Solar irradiance at sensor (µW/cm²/nm)
            - path_radiance: Single scatter path radiance (µW/cm²/nm/sr)
            - mult_scatter: Multiple scatter radiance (µW/cm²/nm/sr)
            - grnd_reflect: Ground reflected radiance (µW/cm²/nm/sr)
            - drct_reflect: Direct reflected radiance (µW/cm²/nm/sr)
            - total_transmission: Total transmission (ground-to-sensor)
            - total_radiance: Total radiance (µW/cm²/nm/sr)
            - solar_scatter: Combined solar scatter (µW/cm²/nm/sr)
            - h2o_column: Water vapor column (g/cm²)
            - solar_zenith: Solar zenith angle (degrees)

        Notes
        -----
        The .chn file format (from MODTRAN6 User Manual):
        - Column 0: Channel center wavelength (nm)
        - Column 4: Total radiance (W/cm²/sr)
        - Column 6: Channel radiance (W/cm²/sr)
        - Column 8: Channel width (nm)
        - Column 14: Multiple scatter solar (W/cm²/sr)
        - Column 15: Single scatter solar (W/cm²/sr)
        - Column 16: Ground reflected (W/cm²/sr)
        - Column 17: Direct reflected (W/cm²/sr)
        - Column 20: Solar at observer (W/cm²)
        - Column 21: Total transmission

        Unit conversions:
        - Divide radiances by channel width to get per-nm values
        - Multiply by 1E6 to convert W → µW
        """
        # Determine .chn file path
        if chn_file_path is None:
            if self.case_name is None:
                raise RuntimeError("No case name available. Run configure_case() first or provide chn_file_path.")
            chn_file_path = Path(f"{self.case_name}.chn")
        else:
            chn_file_path = Path(chn_file_path)

        if not chn_file_path.exists():
            raise FileNotFoundError(f"Channel file not found: {chn_file_path}")

        # Parse .chn file
        wavelength = []
        channel_width = []
        total_transmission = []
        total_radiance = []
        channel_radiance = []
        mult_scatter = []
        sing_scatter = []
        grnd_reflect = []
        drct_reflect = []
        solar_at_obs = []

        with open(chn_file_path, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                # Skip 5-line header
                if i < 5:
                    continue

                # Parse whitespace-separated tokens
                toks = re.findall(r"[\S]+", line.strip())

                if len(toks) < 22:
                    # Skip incomplete lines
                    continue

                # Extract values (with unit conversions)
                wl = float(toks[0])           # Wavelength (nm)
                wid = float(toks[8])          # Channel width (nm)
                tt = float(toks[21])          # Total transmission
                trd = float(toks[4]) * 1E6    # Total radiance: W/cm²/sr → µW/cm²/sr
                crd = float(toks[6]) / wid * 1E6  # Channel radiance: W/cm²/sr → µW/cm²/nm/sr
                ms = float(toks[14]) / wid * 1E6  # Multiple scatter: W/cm²/sr → µW/cm²/nm/sr
                ss = float(toks[15]) / wid * 1E6  # Single scatter: W/cm²/sr → µW/cm²/nm/sr
                gf = float(toks[16]) / wid * 1E6  # Ground reflected: W/cm²/sr → µW/cm²/nm/sr
                df = float(toks[17]) / wid * 1E6  # Direct reflected: W/cm²/sr → µW/cm²/nm/sr
                sobs = float(toks[20]) / wid * 1E6  # Solar at observer: W/cm² → µW/cm²/nm

                wavelength.append(wl)
                channel_width.append(wid)
                total_transmission.append(tt)
                total_radiance.append(trd)
                channel_radiance.append(crd)
                mult_scatter.append(ms)
                sing_scatter.append(ss)
                grnd_reflect.append(gf)
                drct_reflect.append(df)
                solar_at_obs.append(sobs)

        # Convert to numpy arrays
        wavelength = np.array(wavelength)
        channel_width = np.array(channel_width)
        total_transmission = np.array(total_transmission)
        total_radiance = np.array(total_radiance)
        channel_radiance = np.array(channel_radiance)
        mult_scatter = np.array(mult_scatter)
        sing_scatter = np.array(sing_scatter)
        grnd_reflect = np.array(grnd_reflect)
        drct_reflect = np.array(drct_reflect)
        solar_at_obs = np.array(solar_at_obs)

        # Combined solar scatter (multiple + single)
        solar_scatter = mult_scatter + sing_scatter

        # Placeholder spherical albedo (not in .chn file)
        spherical_albedo = np.zeros_like(wavelength)

        results = {
            'wavelength': wavelength,
            'channel_width': channel_width,
            'solar_irradiance': solar_at_obs,
            'path_radiance': sing_scatter,
            'mult_scatter': mult_scatter,
            'grnd_reflect': grnd_reflect,
            'drct_reflect': drct_reflect,
            'total_transmission': total_transmission,
            'spherical_albedo': spherical_albedo,
            'total_radiance': total_radiance,
            'channel_radiance': channel_radiance,
            'solar_scatter': solar_scatter,
            'h2o_column': self.h2o_column_gcm2 if self.h2o_column_gcm2 is not None else 0.0,
            'solar_zenith': self.solar_zenith_deg if self.solar_zenith_deg is not None else 0.0,
        }

        return results

    def get_results(
        self,
        case_idx: int = 0,
        use_chn_file: Optional[bool] = None,
        chn_file_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, np.ndarray]:
        """
        Extract results from MODTRAN6 output

        Parameters
        ----------
        case_idx : int
            Case index to extract (default: 0)
        use_chn_file : bool, optional
            If True, parse .chn file for sensor-convolved data.
            If False, use API for high-resolution unconvolved data.
            If None (default), use .chn file if sensor_filter_path was configured,
            otherwise use API.
        chn_file_path : str or Path, optional
            Path to .chn file (only used if use_chn_file=True)

        Returns
        -------
        results : dict
            Dictionary containing spectral data. Content varies by method:

            .chn file method (sensor-convolved):
            - wavelength: Channel center wavelengths (nm)
            - channel_width: Channel widths (nm)
            - solar_irradiance: Solar irradiance at sensor (µW/cm²/nm)
            - path_radiance: Single scatter path radiance (µW/cm²/nm/sr)
            - mult_scatter: Multiple scatter radiance (µW/cm²/nm/sr)
            - grnd_reflect: Ground reflected radiance (µW/cm²/nm/sr)
            - drct_reflect: Direct reflected radiance (µW/cm²/nm/sr)
            - total_transmission: Total transmission (ground-to-sensor)
            - total_radiance: Total radiance (µW/cm²/nm/sr)
            - channel_radiance: Channel-integrated radiance (µW/cm²/nm/sr)
            - solar_scatter: Combined solar scatter (µW/cm²/nm/sr)
            - spherical_albedo: Spherical albedo (placeholder)
            - h2o_column: Water vapor column (g/cm²)
            - solar_zenith: Solar zenith angle (degrees)

            API method (high-resolution unconvolved):
            - wavelength: Wavelength vector (nm)
            - wavenumber: Wavenumber vector (cm⁻¹)
            - solar_irradiance: Solar irradiance at sensor altitude (µW/cm²/nm)
            - path_radiance: Single scatter path radiance (µW/cm²/nm/sr)
            - mult_scatter: Multiple scatter radiance (µW/cm²/nm/sr)
            - grnd_reflect: Ground reflected radiance (µW/cm²/nm/sr)
            - drct_reflect: Direct reflected radiance (µW/cm²/nm/sr)
            - total_transmission: Total transmission (ground-to-sensor)
            - spherical_albedo: Spherical albedo (placeholder)
            - total_radiance: Total radiance (µW/cm²/nm/sr)

        Notes
        -----
        For atmospheric correction LUT generation with specific sensor bands
        (e.g., AVIRIS-NG), use .chn file method to get sensor-convolved data
        that matches the instrument's spectral response.

        The API method provides high-resolution unconvolved data useful for
        general radiative transfer calculations but not matched to sensor bands.
        """
        # Determine which method to use
        if use_chn_file is None:
            # Auto-detect: use .chn file if sensor filter was configured
            use_chn_file = (self.sensor_filter_path is not None)

        if use_chn_file:
            # Use .chn file for sensor-convolved data
            return self.read_chn_file(chn_file_path)

        # Use API for high-resolution unconvolved data
        if self.modlib is None:
            raise RuntimeError("No MODTRAN case configured.")

        # Get output structure
        modout = self.modlib.caseOutput(case_idx)

        # Get spectral data size
        num_freq = modout.spectra.radiance.num_freq

        # Extract wavenumber (cm⁻¹)
        freq_ptr = pymod.DoubleArray_frompointer(modout.spectra.radiance.freq)
        wavenumber = np.array([freq_ptr[i] for i in range(num_freq)])

        # Convert wavenumber to wavelength (nm)
        wavelength = 1e7 / wavenumber  # λ (nm) = 10^7 / ν (cm⁻¹)

        # Extract radiance data
        # Total radiance (W/cm²/sr/cm⁻¹) → convert to µW/cm²/nm/sr
        tot_rad_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.total_rad)
        total_radiance_wcm2 = np.array([tot_rad_ptr[i] for i in range(num_freq)])

        # Convert units: W/cm²/sr/cm⁻¹ → µW/cm²/nm/sr
        # dν/dλ = -(10^7)/λ² for λ in nm, ν in cm⁻¹
        # L_λ = L_ν × |dν/dλ| = L_ν × (10^7/λ²)
        # Also convert W → µW (×10⁶)
        total_radiance = total_radiance_wcm2 * (1e7 / wavelength**2) * 1e6

        # Extract transmission
        trans_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.tot_trans)
        total_transmission = np.array([trans_ptr[i] for i in range(num_freq)])

        # Path radiance (single scatter solar component)
        # SING_SCAT in MODTRAN API output
        path_rad_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.sing_scat)
        path_radiance_wcm2 = np.array([path_rad_ptr[i] for i in range(num_freq)])
        path_radiance = path_radiance_wcm2 * (1e7 / wavelength**2) * 1e6

        # Multiple scatter radiance
        mult_scat_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.mult_scat)
        mult_scatter_wcm2 = np.array([mult_scat_ptr[i] for i in range(num_freq)])
        mult_scatter = mult_scatter_wcm2 * (1e7 / wavelength**2) * 1e6

        # Ground reflected radiance
        grnd_rflt_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.grnd_rflt)
        grnd_reflect_wcm2 = np.array([grnd_rflt_ptr[i] for i in range(num_freq)])
        grnd_reflect = grnd_reflect_wcm2 * (1e7 / wavelength**2) * 1e6

        # Direct reflected radiance
        drct_rflt_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.drct_rflt)
        drct_reflect_wcm2 = np.array([drct_rflt_ptr[i] for i in range(num_freq)])
        drct_reflect = drct_reflect_wcm2 * (1e7 / wavelength**2) * 1e6

        # Solar irradiance at observer (H1ALT)
        # SOL_AT_OBS in MODTRAN API output
        sol_irrad_ptr = pymod.FloatArray_frompointer(modout.spectra.radiance.sol_at_obs)
        solar_irrad_wcm2 = np.array([sol_irrad_ptr[i] for i in range(num_freq)])
        solar_irradiance = solar_irrad_wcm2 * (1e7 / wavelength**2) * 1e6

        # Spherical albedo
        # For now, approximate from transmission
        # TODO: Check if there's a specific field for spherical albedo
        # in modout.spectra.transmittance
        spherical_albedo = np.zeros_like(wavelength)

        results = {
            'wavelength': wavelength,
            'wavenumber': wavenumber,
            'solar_irradiance': solar_irradiance,
            'path_radiance': path_radiance,
            'mult_scatter': mult_scatter,
            'grnd_reflect': grnd_reflect,
            'drct_reflect': drct_reflect,
            'total_transmission': total_transmission,
            'spherical_albedo': spherical_albedo,
            'total_radiance': total_radiance,
        }

        return results

    def cleanup(self):
        """Clean up MODTRAN library resources"""
        if self.modlib is not None:
            del self.modlib
            self.modlib = None

    def __del__(self):
        """Destructor - ensure cleanup"""
        self.cleanup()


def run_modtran_case(
    sensor_altitude_km: float,
    ground_altitude_km: float,
    solar_zenith_deg: float,
    h2o_column_gcm2: float = 2.0,
    visibility_km: float = 35.0,
    surface_albedo: float = 0.0,
    **kwargs
) -> Dict[str, np.ndarray]:
    """
    Convenience function to run a single MODTRAN6 case

    Parameters
    ----------
    sensor_altitude_km : float
        Sensor altitude (km)
    ground_altitude_km : float
        Ground elevation (km)
    solar_zenith_deg : float
        Solar zenith angle (degrees)
    h2o_column_gcm2 : float
        Water vapor column (g/cm²)
    visibility_km : float
        Visibility (km)
    surface_albedo : float
        Surface albedo (0-1)
    **kwargs
        Additional parameters for configure_case()

    Returns
    -------
    results : dict
        MODTRAN6 output results
    """
    runner = ModtranRunner()

    try:
        # Configure case
        runner.configure_case(
            sensor_altitude_km=sensor_altitude_km,
            ground_altitude_km=ground_altitude_km,
            solar_zenith_deg=solar_zenith_deg,
            h2o_column_gcm2=h2o_column_gcm2,
            visibility_km=visibility_km,
            surface_albedo=surface_albedo,
            **kwargs
        )

        # Run MODTRAN
        success = runner.run()
        if not success:
            raise RuntimeError("MODTRAN6 execution failed")

        # Extract results
        results = runner.get_results()

        return results

    finally:
        runner.cleanup()
