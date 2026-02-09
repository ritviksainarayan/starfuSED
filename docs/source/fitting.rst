==========
SED Fitting
==========

*starfuSED* provides two SED fitting classes: ``SingleSEDFitter`` for single stars and ``BinarySEDFitter`` for binary systems.

How SED Fitting Works
---------------------

The fitting process uses chi-squared minimization:

1. **Grid Search**: For each combination of Teff and log g in the parameter grid:

   a. Load the model spectrum
   b. Interpolate to observation wavelengths
   c. Compute normalization factor from the specified band
   d. Calculate chi-squared: χ² = Σ[(F_obs - F_model)² / σ²]

2. **Best Fit**: Select parameters with minimum chi-squared

3. **Radius Calculation**: From the normalization factor (R/d)²:

   R = sqrt(norm) × d

Single Star Fitting
-------------------

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from starfused import Photometry, SingleSEDFitter, plot_single_sed

   # Get photometry
   phot = Photometry.query(name="HD 12345", filters=['GALEX', 'SDSS', '2MASS', 'WISE'])
   phot = Photometry.dust_correction(phot)

   # Define parameters
   params = {
       'modelname': 'ck04',
       'teff_min': 5000,
       'teff_max': 7000,
       'teff_step': 250,
       'logg_min': 3.5,
       'logg_max': 5.0,
       'logg_step': 0.5,
       'metallicity': 0.0,
       'norm_band': 'SDSS:r'
   }

   # Create fitter and fit
   fitter = SingleSEDFitter(phot, distance_pc=50, source_params=params)
   result = fitter.fit()

Parameter Dictionary
~~~~~~~~~~~~~~~~~~~~

The ``source_params`` dictionary requires:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Required
     - Description
   * - ``modelname``
     - Yes
     - Model grid: ``'ck04'``, ``'phoenix'``, ``'koester'``, ``'bt-settl'``
   * - ``teff_min``
     - Yes
     - Minimum temperature to search (K)
   * - ``teff_max``
     - Yes
     - Maximum temperature to search (K)
   * - ``teff_step``
     - No
     - Temperature step size (default: 250 K)
   * - ``logg_min``
     - Yes
     - Minimum log g to search
   * - ``logg_max``
     - Yes
     - Maximum log g to search
   * - ``logg_step``
     - No
     - log g step size (default: 0.5)
   * - ``metallicity``
     - Yes
     - Fixed metallicity [M/H] (use ``None`` for Koester)
   * - ``norm_band``
     - Yes
     - Filter name for normalization (e.g., ``'SDSS:r'``, ``'2MASS:H'``)

Adaptive Fitting
~~~~~~~~~~~~~~~~

For better precision, use adaptive grid refinement:

.. code-block:: python

   result = fitter.fit_adaptive(
       n_refine=2,        # Number of refinement iterations
       refine_factor=4,   # Step size reduction per iteration
       coarse_factor=4    # Initial coarsening factor
   )

The adaptive method:

1. Starts with a coarse grid (steps × coarse_factor)
2. Finds approximate best-fit
3. Refines grid around best-fit (steps ÷ refine_factor)
4. Repeats n_refine times

Result Dictionary
~~~~~~~~~~~~~~~~~

.. code-block:: python

   result = {
       'teff': 6000,              # Best-fit temperature (K)
       'logg': 4.5,               # Best-fit log g
       'metallicity': 0.0,        # Metallicity used
       'norm': 1.23e-20,          # Normalization factor (R/d)²
       'radius_rsun': 1.05,       # Radius in solar radii
       'radius_rearth': 115.2,    # Radius in Earth radii
       'radius_rjup': 10.5,       # Radius in Jupiter radii
       'spectrum': DataFrame,     # Best-fit model spectrum
       'model_flux': array,       # Model flux at observed wavelengths
       'chi2': 15.3,              # Chi-squared
       'reduced_chi2': 1.2,       # Reduced chi-squared
       'n_data': 15,              # Number of data points
       'n_params': 1,             # Number of free parameters
       'distance_pc': 50,         # Distance used

       # Added by fit_mc() (if called):
       'mc_errors': {             # Asymmetric error bounds
           'teff_err': (250, 500),
           'logg_err': (0.5, 0.5),
           'radius_rsun_err': (0.02, 0.03),
           # ... plus radius_rearth_err, radius_rjup_err
           'n_iter': 98,          # Successful iterations
           'n_failed': 2          # Failed iterations
       },
       'mc_distributions': {      # Raw parameter arrays
           'teff': array,
           'logg': array,
           'radius_rsun': array,
           # ... plus radius_rearth, radius_rjup
       }
   }

Binary System Fitting
---------------------

For systems with two stellar components (e.g., WD + M dwarf):

.. code-block:: python

   from starfused import Photometry, BinarySEDFitter, plot_binary_sed

   phot = Photometry.query(name="TIC 12345678", filters=['GALEX', 'SDSS', '2MASS', 'WISE'])
   phot = Photometry.dust_correction(phot)

   # Hot component (dominates UV)
   source1 = {
       'modelname': 'koester',
       'teff_min': 10000,
       'teff_max': 30000,
       'teff_step': 500,
       'logg_min': 7.0,
       'logg_max': 9.0,
       'logg_step': 0.25,
       'metallicity': None,
       'norm_band': 'GALEX:NUV'
   }

   # Cool component (dominates IR)
   source2 = {
       'modelname': 'bt-settl',
       'teff_min': 2500,
       'teff_max': 4000,
       'teff_step': 100,
       'logg_min': 4.5,
       'logg_max': 5.5,
       'logg_step': 0.5,
       'metallicity': 0.0,
       'norm_band': '2MASS:H'
   }

   fitter = BinarySEDFitter(
       phot,
       distance_pc=100,
       source1_params=source1,
       source2_params=source2
   )

   result = fitter.fit_adaptive()

Binary Normalization Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The binary fitter anchors each component to a different band:

- **Source 1**: Normalized to its ``norm_band`` (typically UV for hot component)
- **Source 2**: Either normalized to its ``norm_band`` or fit via least-squares to the residual after subtracting source 1

This approach works best when the two components dominate at different wavelengths.

Binary Result Dictionary
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   result = {
       'source1': {
           'teff': 15000,
           'logg': 8.0,
           'metallicity': None,
           'norm': 2.3e-22,
           'radius_rsun': 0.012,
           'radius_rearth': 1.3,
           'radius_rjup': 0.12,
           'spectrum': DataFrame
       },
       'source2': {
           'teff': 3200,
           'logg': 5.0,
           'metallicity': 0.0,
           'norm': 1.5e-21,
           'radius_rsun': 0.25,
           'radius_rearth': 27.5,
           'radius_rjup': 2.5,
           'spectrum': DataFrame
       },
       'chi2': 12.5,
       'reduced_chi2': 1.1,
       'combined_flux': array,
       'distance_pc': 100
   }

Monte Carlo Uncertainties
-------------------------

After fitting, you can estimate uncertainties by perturbing the observed fluxes
within their measurement errors and re-fitting. The ``fit_mc()`` method draws
perturbations from a truncated normal distribution (clipped at ±σ_clip) and
reports asymmetric error bounds from the resulting parameter distributions.

.. code-block:: python

   # First, perform the fit
   fitter = SingleSEDFitter(phot, distance_pc=50, source_params=params)
   fitter.fit_adaptive(n_refine=2)

   # Then run Monte Carlo to estimate 3σ uncertainties
   result = fitter.fit_mc(
       n_iter=100,      # Number of perturbation iterations
       sigma_clip=3,    # Perturb within ±3σ
       seed=42          # For reproducibility
   )

   # Asymmetric error bounds
   err = result['mc_errors']
   print(f"Teff = {result['teff']} (+{err['teff_err'][1]:.0f}/-{err['teff_err'][0]:.0f}) K")
   print(f"log g = {result['logg']:.2f} (+{err['logg_err'][1]:.2f}/-{err['logg_err'][0]:.2f})")
   print(f"Radius = {result['radius_rsun']:.4f} (+{err['radius_rsun_err'][1]:.4f}/-{err['radius_rsun_err'][0]:.4f}) R_sun")

   # Raw distributions are available for custom analysis
   teff_distribution = result['mc_distributions']['teff']

The same method works for binary fitting:

.. code-block:: python

   fitter = BinarySEDFitter(phot, distance_pc=100,
                            source1_params=source1, source2_params=source2)
   fitter.fit_adaptive()
   result = fitter.fit_mc(n_iter=100, seed=42)

   # Binary errors are prefixed with s1_ and s2_
   err = result['mc_errors']
   print(f"Source 1 Teff = {result['source1']['teff']} "
         f"(+{err['s1_teff_err'][1]:.0f}/-{err['s1_teff_err'][0]:.0f}) K")
   print(f"Source 2 Teff = {result['source2']['teff']} "
         f"(+{err['s2_teff_err'][1]:.0f}/-{err['s2_teff_err'][0]:.0f}) K")

How It Works
~~~~~~~~~~~~

1. Each iteration draws random perturbations from N(0, 1), clips them to ±σ_clip, and scales by the measurement uncertainties.
2. The perturbed fluxes are fit using the full original parameter grid via ``fit()``.
3. The min/max of each parameter across all successful iterations gives the asymmetric error bounds.
4. Since the perturbations span ±σ_clip σ, the reported bounds represent σ_clip-σ uncertainties.

.. note::

   Model spectra are cached from the initial fit, so only the chi-squared
   evaluation is repeated each iteration. This makes ``fit_mc()`` much faster
   than running ``fit()`` from scratch each time.

Convenience Functions
---------------------

For quick one-liner fitting:

.. code-block:: python

   from starfused import fit_single_sed, fit_binary_sed

   # Single star
   result = fit_single_sed(phot, distance_pc=50, source_params=params)

   # Binary
   result = fit_binary_sed(phot, distance_pc=100,
                           source1_params=source1, source2_params=source2)

Tips for Good Fits
------------------

Choosing the Normalization Band
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Choose a band where the star dominates (not contaminated by binary companion)
- Avoid bands with large uncertainties
- For binaries, pick bands where each component clearly dominates

.. code-block:: python

   # Good: WD dominates in UV, companion dominates in IR
   source1 = {..., 'norm_band': 'GALEX:NUV'}  # WD
   source2 = {..., 'norm_band': '2MASS:H'}    # M dwarf

   # Bad: Both contribute at optical wavelengths
   source1 = {..., 'norm_band': 'SDSS:r'}  # Uncertain contribution
   source2 = {..., 'norm_band': 'SDSS:i'}  # Uncertain contribution

Setting Parameter Ranges
~~~~~~~~~~~~~~~~~~~~~~~~

- Start with wide ranges, then narrow based on initial fits
- Use adaptive fitting to efficiently search large parameter spaces
- Check that best-fit is not at edge of grid (suggests wrong range)

API Reference
-------------

.. autoclass:: starfused.SEDFitter
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.SingleSEDFitter
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.BinarySEDFitter
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: starfused.fit_single_sed

.. autofunction:: starfused.fit_binary_sed
