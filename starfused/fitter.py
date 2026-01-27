"""
SED fitting module for single and binary stellar systems.
"""

import numpy as np
import pandas as pd
from .spectra import StellarModel


class SEDFitter:
    """
    Base class for SED fitting with common functionality.

    Parameters
    ----------
    photometry_df : pandas.DataFrame
        Photometric data with columns:
        - 'sed_wave': Wavelength in Angstroms
        - 'sed_flux': Flux in erg/s/cm^2/A
        - 'sed_eflux': Flux uncertainty in erg/s/cm^2/A
        - 'sed_filter': Filter name (used for normalization band matching)
    distance_pc : float
        Distance to the system in parsecs.
    verbose : bool, optional
        Print progress and results (default: True).
    """

    def __init__(self, photometry_df, distance_pc, verbose=True):
        self.photometry_df = photometry_df
        self.distance_pc = distance_pc
        self.verbose = verbose

        # Extract observation data
        self.obs_waves = photometry_df['sed_wave'].values
        self.obs_flux = photometry_df['sed_flux'].values
        self.obs_eflux = photometry_df['sed_eflux'].values

        # Will store fit results
        self.result = None

    def _validate_source_params(self, params, source_name="source"):
        """Validate that required keys are present in source parameters."""
        required_keys = ['modelname', 'teff_min', 'teff_max', 'logg_min', 'logg_max',
                         'metallicity', 'norm_band']
        for key in required_keys:
            if key not in params:
                raise ValueError(f"{source_name} missing required key: {key}")

    def _build_param_grid(self, params):
        """Build temperature and logg grids from parameters."""
        teff_step = params.get('teff_step', 250)
        logg_step = params.get('logg_step', 0.5)

        teff_grid = list(range(
            int(params['teff_min']),
            int(params['teff_max']) + 1,
            int(teff_step)
        ))
        logg_grid = list(np.arange(
            params['logg_min'],
            params['logg_max'] + logg_step / 2,
            logg_step
        ))
        return teff_grid, logg_grid

    def _get_norm_wavelength(self, norm_band, source_name="source"):
        """Get normalization wavelength and flux for a given band."""
        norm_mask = self.photometry_df['sed_filter'].str.contains(
            norm_band, case=False, na=False
        )
        if norm_mask.sum() == 0:
            raise ValueError(
                f"No filter matching '{norm_band}' found for {source_name} normalization"
            )
        norm_wave = self.photometry_df.loc[norm_mask, 'sed_wave'].values[0]
        norm_flux = self.photometry_df.loc[norm_mask, 'sed_flux'].values[0]
        return norm_wave, norm_flux

    def _get_model_flux(self, model, cache, teff, logg, metallicity, wave_grid, norm_wave):
        """Load model and interpolate to wavelength grid, with caching."""
        cache_key = (teff, logg, metallicity)
        if cache_key not in cache:
            try:
                spec = model.model.load_model(teff=teff, logg=logg, metallicity=metallicity)
                spec_wave = spec['wavelength'].values
                spec_flux = spec['flux'].values

                # Interpolate to grid and norm wavelength
                flux_at_grid = np.interp(wave_grid, spec_wave, spec_flux)
                flux_at_norm = np.interp(norm_wave, spec_wave, spec_flux)

                # Handle extrapolation for IR if model doesn't extend far enough
                max_wave = spec_wave.max()
                if max_wave < wave_grid.max():
                    # Fit power law to red end for extrapolation
                    red_mask = spec_wave > 0.8 * max_wave
                    if red_mask.sum() > 2:
                        log_w = np.log10(spec_wave[red_mask])
                        log_f = np.log10(np.maximum(spec_flux[red_mask], 1e-30))
                        slope, intercept = np.polyfit(log_w, log_f, 1)

                        extrap_mask = wave_grid > max_wave
                        flux_at_grid[extrap_mask] = 10**(
                            slope * np.log10(wave_grid[extrap_mask]) + intercept
                        )

                cache[cache_key] = {
                    'flux_at_grid': flux_at_grid,
                    'flux_at_norm': flux_at_norm,
                    'spectrum': spec
                }
            except Exception:
                cache[cache_key] = None
        return cache[cache_key]

    def _calculate_radius(self, norm):
        """Calculate radius from normalization factor and distance."""
        from astropy import units as u

        distance_cm = (self.distance_pc * u.pc).to(u.cm).value
        radius_cm = np.sqrt(norm) * distance_cm

        return {
            'radius_rsun': (radius_cm * u.cm).to(u.R_sun).value,
            'radius_rearth': (radius_cm * u.cm).to(u.R_earth).value,
            'radius_rjup': (radius_cm * u.cm).to(u.R_jup).value,
        }

    def fit(self):
        """Perform the SED fit. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement fit()")


class SingleSEDFitter(SEDFitter):
    """
    Fit a single star SED.

    Parameters
    ----------
    photometry_df : pandas.DataFrame
        Photometric data with columns:
        - 'sed_wave': Wavelength in Angstroms
        - 'sed_flux': Flux in erg/s/cm^2/A
        - 'sed_eflux': Flux uncertainty in erg/s/cm^2/A
        - 'sed_filter': Filter name (used for normalization band matching)
    distance_pc : float
        Distance to the star in parsecs.
    source_params : dict
        Parameters for the stellar model with keys:
        - 'modelname': str - Model grid name ('ck04', 'phoenix', 'koester', 'bt-settl')
        - 'teff_min': float - Minimum effective temperature (K)
        - 'teff_max': float - Maximum effective temperature (K)
        - 'teff_step': float, optional - Temperature step (default: 250 K)
        - 'logg_min': float - Minimum log g
        - 'logg_max': float - Maximum log g
        - 'logg_step': float, optional - log g step (default: 0.5)
        - 'metallicity': float or None - Metallicity [M/H] (None for Koester)
        - 'norm_band': str - Filter name substring for normalization (e.g., 'NUV', '2MASS:H')
    verbose : bool, optional
        Print progress and results (default: True).

    Examples
    --------
    >>> from starfused import Photometry
    >>> from starfused.fitter import SingleSEDFitter
    >>>
    >>> # Get photometry
    >>> phot = Photometry.query(name="HD 12345", filters=['GALEX', 'SDSS', '2MASS', 'WISE'])
    >>> corrected = Photometry.dust_correction(phot)
    >>>
    >>> params = {
    ...     'modelname': 'ck04',
    ...     'teff_min': 5000, 'teff_max': 7000, 'teff_step': 250,
    ...     'logg_min': 3.5, 'logg_max': 5.0, 'logg_step': 0.5,
    ...     'metallicity': 0.0,
    ...     'norm_band': 'SDSS:r'
    ... }
    >>> fitter = SingleSEDFitter(corrected, distance_pc=50, source_params=params)
    >>> result = fitter.fit()
    """

    def __init__(self, photometry_df, distance_pc, source_params, verbose=True):
        super().__init__(photometry_df, distance_pc, verbose)
        self._validate_source_params(source_params, "source_params")
        self.source_params = source_params

        # Initialize model
        self.model = StellarModel(grid=source_params['modelname'])
        self._cache = {}

    def fit(self):
        """
        Perform the single star SED fit.

        Returns
        -------
        dict
            Best-fit results with keys:
            - 'teff': Best-fit effective temperature (K)
            - 'logg': Best-fit log g
            - 'metallicity': Metallicity used
            - 'norm': Normalization factor (R/d)^2
            - 'radius_rsun': Inferred radius in solar radii
            - 'radius_rearth': Inferred radius in Earth radii
            - 'radius_rjup': Inferred radius in Jupiter radii
            - 'spectrum': DataFrame with model spectrum
            - 'model_flux': Model flux at observed wavelengths
            - 'chi2': Best-fit chi-squared
            - 'reduced_chi2': Reduced chi-squared
            - 'n_data': Number of data points
            - 'n_params': Number of free parameters
            - 'distance_pc': Distance used in fit
        """
        # Build parameter grid
        teff_grid, logg_grid = self._build_param_grid(self.source_params)

        # Get normalization wavelength
        norm_wave, norm_flux = self._get_norm_wavelength(
            self.source_params['norm_band'], "source"
        )

        n_models = len(teff_grid) * len(logg_grid)

        if self.verbose:
            print(f"Model ({self.source_params['modelname']}): {len(teff_grid)} Teff x {len(logg_grid)} logg")
            print(f"  Teff: {teff_grid[0]}-{teff_grid[-1]} K, logg: {logg_grid[0]:.1f}-{logg_grid[-1]:.1f}")
            print(f"  Normalization: {self.source_params['norm_band']} at {norm_wave:.0f} A")
            print(f"Total models to search: {n_models}")

        # Grid search
        best_chi2 = np.inf
        best_params = None
        tested = 0

        for teff in teff_grid:
            for logg in logg_grid:
                tested += 1

                # Get model flux
                model_result = self._get_model_flux(
                    self.model, self._cache, teff, logg,
                    self.source_params['metallicity'], self.obs_waves, norm_wave
                )
                if model_result is None:
                    continue

                flux_at_grid = model_result['flux_at_grid']
                flux_at_norm = model_result['flux_at_norm']

                # Anchor normalization to norm_band
                # norm = observed_flux / model_flux at norm wavelength
                if flux_at_norm <= 0:
                    continue
                norm = norm_flux / flux_at_norm

                if norm <= 0:
                    continue

                # Compute model flux and chi-squared
                model_flux = norm * flux_at_grid
                chi2 = np.sum(((self.obs_flux - model_flux) / self.obs_eflux)**2)

                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_params = {
                        'teff': teff,
                        'logg': logg,
                        'norm': norm,
                        'spectrum': model_result['spectrum'],
                        'model_flux': model_flux,
                    }

        if best_params is None:
            raise RuntimeError("No valid fits found. Check parameter ranges and normalization band.")

        # Calculate radius
        radii = self._calculate_radius(best_params['norm'])

        n_data = len(self.obs_flux)
        n_params = 1  # normalization (Teff/logg selected from grid)
        reduced_chi2 = best_chi2 / (n_data - n_params)

        self.result = {
            'teff': best_params['teff'],
            'logg': best_params['logg'],
            'metallicity': self.source_params['metallicity'],
            'norm': best_params['norm'],
            'radius_rsun': radii['radius_rsun'],
            'radius_rearth': radii['radius_rearth'],
            'radius_rjup': radii['radius_rjup'],
            'spectrum': best_params['spectrum'],
            'model_flux': best_params['model_flux'],
            'chi2': best_chi2,
            'reduced_chi2': reduced_chi2,
            'n_data': n_data,
            'n_params': n_params,
            'distance_pc': self.distance_pc,
        }

        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Best-fit Single Star Model:")
            print(f"  Tested {tested} models")
            print(f"  Model: {self.source_params['modelname']}")
            print(f"  Teff    = {self.result['teff']} K")
            print(f"  log g   = {self.result['logg']}")
            print(f"  Radius  = {radii['radius_rsun']:.4f} R_sun = {radii['radius_rearth']:.2f} R_earth")
            print(f"  Reduced chi2 = {reduced_chi2:.2f}")
            print(f"{'='*60}")

        return self.result


class BinarySEDFitter(SEDFitter):
    """
    Fit a binary system SED with two stellar components simultaneously.

    Both sources share the same distance, but each has independent stellar
    parameters and normalization bandpass.

    Parameters
    ----------
    photometry_df : pandas.DataFrame
        Photometric data with columns:
        - 'sed_wave': Wavelength in Angstroms
        - 'sed_flux': Flux in erg/s/cm^2/A
        - 'sed_eflux': Flux uncertainty in erg/s/cm^2/A
        - 'sed_filter': Filter name (used for normalization band matching)
    distance_pc : float
        Distance to the system in parsecs (shared by both sources).
    source1_params : dict
        Parameters for the first source with keys:
        - 'modelname': str - Model grid name ('ck04', 'phoenix', 'koester', 'bt-settl')
        - 'teff_min': float - Minimum effective temperature (K)
        - 'teff_max': float - Maximum effective temperature (K)
        - 'teff_step': float, optional - Temperature step (default: 250 K)
        - 'logg_min': float - Minimum log g
        - 'logg_max': float - Maximum log g
        - 'logg_step': float, optional - log g step (default: 0.5)
        - 'metallicity': float or None - Metallicity [M/H] (None for Koester)
        - 'norm_band': str - Filter name substring for normalization (e.g., 'NUV', '2MASS:H')
    source2_params : dict
        Parameters for the second source (same structure as source1_params).
    verbose : bool, optional
        Print progress and results (default: True).

    Examples
    --------
    >>> from starfused import Photometry
    >>> from starfused.fitter import BinarySEDFitter
    >>>
    >>> # Get photometry
    >>> phot = Photometry.query(name="TIC 12345678", filters=['GALEX', 'SDSS', '2MASS', 'WISE'])
    >>> corrected = Photometry.dust_correction(phot)
    >>>
    >>> source1 = {
    ...     'modelname': 'koester',
    ...     'teff_min': 10000, 'teff_max': 20000, 'teff_step': 250,
    ...     'logg_min': 7.5, 'logg_max': 8.5, 'logg_step': 0.25,
    ...     'metallicity': None,
    ...     'norm_band': 'NUV'
    ... }
    >>> source2 = {
    ...     'modelname': 'bt-settl',
    ...     'teff_min': 1000, 'teff_max': 3000, 'teff_step': 100,
    ...     'logg_min': 4.5, 'logg_max': 5.5, 'logg_step': 0.5,
    ...     'metallicity': 0.0,
    ...     'norm_band': '2MASS:H'
    ... }
    >>> fitter = BinarySEDFitter(corrected, distance_pc=100,
    ...                          source1_params=source1, source2_params=source2)
    >>> result = fitter.fit()
    """

    def __init__(self, photometry_df, distance_pc, source1_params, source2_params, verbose=True):
        super().__init__(photometry_df, distance_pc, verbose)
        self._validate_source_params(source1_params, "source1_params")
        self._validate_source_params(source2_params, "source2_params")
        self.source1_params = source1_params
        self.source2_params = source2_params

        # Initialize models
        self.model1 = StellarModel(grid=source1_params['modelname'])
        self.model2 = StellarModel(grid=source2_params['modelname'])
        self._cache1 = {}
        self._cache2 = {}

    def fit(self):
        """
        Perform the binary SED fit.

        Returns
        -------
        dict
            Best-fit results with keys:
            - 'source1': dict with fitted parameters for source 1
                - 'teff': Best-fit effective temperature (K)
                - 'logg': Best-fit log g
                - 'metallicity': Metallicity used
                - 'norm': Normalization factor (R/d)^2
                - 'radius_rsun': Inferred radius in solar radii
                - 'radius_rearth': Inferred radius in Earth radii
                - 'radius_rjup': Inferred radius in Jupiter radii
                - 'spectrum': DataFrame with model spectrum
            - 'source2': dict with fitted parameters for source 2 (same structure)
            - 'chi2': Best-fit chi-squared
            - 'reduced_chi2': Reduced chi-squared
            - 'n_data': Number of data points
            - 'n_params': Number of free parameters
            - 'combined_flux': Combined model flux at observed wavelengths
            - 'distance_pc': Distance used in fit
        """
        # Build parameter grids
        s1_teff_grid, s1_logg_grid = self._build_param_grid(self.source1_params)
        s2_teff_grid, s2_logg_grid = self._build_param_grid(self.source2_params)

        # Get normalization wavelengths
        norm1_wave, norm1_flux = self._get_norm_wavelength(
            self.source1_params['norm_band'], "source 1"
        )
        norm2_wave, norm2_flux = self._get_norm_wavelength(
            self.source2_params['norm_band'], "source 2"
        )

        n_models = len(s1_teff_grid) * len(s1_logg_grid) * len(s2_teff_grid) * len(s2_logg_grid)

        if self.verbose:
            print(f"Source 1 ({self.source1_params['modelname']}): {len(s1_teff_grid)} Teff x {len(s1_logg_grid)} logg")
            print(f"  Teff: {s1_teff_grid[0]}-{s1_teff_grid[-1]} K, logg: {s1_logg_grid[0]:.1f}-{s1_logg_grid[-1]:.1f}")
            print(f"  Normalization: {self.source1_params['norm_band']} at {norm1_wave:.0f} A")
            print(f"Source 2 ({self.source2_params['modelname']}): {len(s2_teff_grid)} Teff x {len(s2_logg_grid)} logg")
            print(f"  Teff: {s2_teff_grid[0]}-{s2_teff_grid[-1]} K, logg: {s2_logg_grid[0]:.1f}-{s2_logg_grid[-1]:.1f}")
            print(f"  Normalization: {self.source2_params['norm_band']} at {norm2_wave:.0f} A")
            print(f"Total models to search: {n_models}")

        # Grid search
        best_chi2 = np.inf
        best_params = None
        tested = 0

        for s1_teff in s1_teff_grid:
            for s1_logg in s1_logg_grid:
                # Get source 1 model
                s1_result = self._get_model_flux(
                    self.model1, self._cache1, s1_teff, s1_logg,
                    self.source1_params['metallicity'], self.obs_waves, norm1_wave
                )
                if s1_result is None:
                    continue

                # Anchor source 1 normalization to its norm_band
                # norm1 = observed_flux / model_flux at norm wavelength
                s1_flux_at_norm = s1_result['flux_at_norm']
                if s1_flux_at_norm <= 0:
                    continue
                norm1 = norm1_flux / s1_flux_at_norm

                if norm1 <= 0:
                    continue

                for s2_teff in s2_teff_grid:
                    for s2_logg in s2_logg_grid:
                        tested += 1

                        # Get source 2 model
                        s2_result = self._get_model_flux(
                            self.model2, self._cache2, s2_teff, s2_logg,
                            self.source2_params['metallicity'], self.obs_waves, norm2_wave
                        )
                        if s2_result is None:
                            continue

                        # Anchor source 2 normalization to its norm_band
                        # But account for source 1 contribution at that wavelength
                        s1_flux_grid = s1_result['flux_at_grid']
                        s2_flux_grid = s2_result['flux_at_grid']
                        s2_flux_at_norm = s2_result['flux_at_norm']

                        if s2_flux_at_norm <= 0:
                            continue

                        # Source 1 contribution at source 2's norm wavelength
                        s1_at_norm2_wave = np.interp(norm2_wave, self.obs_waves, s1_flux_grid * norm1)

                        # Residual flux at norm2 band after subtracting source 1
                        residual_at_norm2 = norm2_flux - s1_at_norm2_wave

                        # norm2 = residual / model_flux at norm wavelength
                        norm2 = residual_at_norm2 / s2_flux_at_norm

                        # Skip negative normalizations (unphysical - source 1 already exceeds observed)
                        if norm2 <= 0:
                            continue

                        # Compute combined model and chi-squared
                        combined_flux = norm1 * s1_flux_grid + norm2 * s2_flux_grid
                        chi2 = np.sum(((self.obs_flux - combined_flux) / self.obs_eflux)**2)

                        if chi2 < best_chi2:
                            best_chi2 = chi2
                            best_params = {
                                's1_teff': s1_teff,
                                's1_logg': s1_logg,
                                's1_norm': norm1,
                                's1_spectrum': s1_result['spectrum'],
                                's2_teff': s2_teff,
                                's2_logg': s2_logg,
                                's2_norm': norm2,
                                's2_spectrum': s2_result['spectrum'],
                                'combined_flux': combined_flux,
                            }

        if best_params is None:
            raise RuntimeError("No valid fits found. Check parameter ranges and normalization bands.")

        # Calculate radii from normalizations
        s1_radii = self._calculate_radius(best_params['s1_norm'])
        s2_radii = self._calculate_radius(best_params['s2_norm'])

        n_data = len(self.obs_flux)
        n_params = 4  # 2 normalizations + implicit Teff/logg selection
        reduced_chi2 = best_chi2 / (n_data - n_params)

        self.result = {
            'source1': {
                'teff': best_params['s1_teff'],
                'logg': best_params['s1_logg'],
                'metallicity': self.source1_params['metallicity'],
                'norm': best_params['s1_norm'],
                'radius_rsun': s1_radii['radius_rsun'],
                'radius_rearth': s1_radii['radius_rearth'],
                'radius_rjup': s1_radii['radius_rjup'],
                'spectrum': best_params['s1_spectrum'],
            },
            'source2': {
                'teff': best_params['s2_teff'],
                'logg': best_params['s2_logg'],
                'metallicity': self.source2_params['metallicity'],
                'norm': best_params['s2_norm'],
                'radius_rsun': s2_radii['radius_rsun'],
                'radius_rearth': s2_radii['radius_rearth'],
                'radius_rjup': s2_radii['radius_rjup'],
                'spectrum': best_params['s2_spectrum'],
            },
            'chi2': best_chi2,
            'reduced_chi2': reduced_chi2,
            'n_data': n_data,
            'n_params': n_params,
            'combined_flux': best_params['combined_flux'],
            'distance_pc': self.distance_pc,
        }

        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Best-fit Binary Model:")
            print(f"  Tested {tested} model combinations")
            print(f"  Source 1 ({self.source1_params['modelname']}):")
            print(f"    Teff    = {self.result['source1']['teff']} K")
            print(f"    log g   = {self.result['source1']['logg']}")
            print(f"    Radius  = {s1_radii['radius_rsun']:.4f} R_sun = {s1_radii['radius_rearth']:.2f} R_earth")
            print(f"  Source 2 ({self.source2_params['modelname']}):")
            print(f"    Teff    = {self.result['source2']['teff']} K")
            print(f"    log g   = {self.result['source2']['logg']}")
            print(f"    Radius  = {s2_radii['radius_rsun']:.4f} R_sun = {s2_radii['radius_rjup']:.2f} R_jup")
            print(f"  Combined Reduced chi2 = {reduced_chi2:.2f}")
            print(f"{'='*60}")

        return self.result


def fit_single_sed(photometry_df, distance_pc, source_params, verbose=True):
    """
    Convenience function to fit a single star SED.

    See SingleSEDFitter for full documentation.
    """
    fitter = SingleSEDFitter(photometry_df, distance_pc, source_params, verbose)
    return fitter.fit()


def fit_binary_sed(photometry_df, distance_pc, source1_params, source2_params, verbose=True):
    """
    Convenience function to fit a binary system SED.

    See BinarySEDFitter for full documentation.
    """
    fitter = BinarySEDFitter(photometry_df, distance_pc, source1_params, source2_params, verbose)
    return fitter.fit()
