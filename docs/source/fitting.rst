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

4. **Uncertainty Estimation**: Using Δχ² = 1 method

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
       coarse_factor=4,   # Initial coarsening factor
       compute_errors=True
   )

The adaptive method:

1. Starts with a coarse grid (steps × coarse_factor)
2. Finds approximate best-fit
3. Refines grid around best-fit (steps ÷ refine_factor)
4. Repeats n_refine times
5. Computes uncertainties from Δχ² = 1 contour

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
       'errors': {                # Only if compute_errors=True
           'teff_err': (250, 250),
           'logg_err': (0.5, 0.5),
           'radius_rsun_err': (0.05, 0.05)
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
       'distance_pc': 100,
       'errors': {...}  # If compute_errors=True
   }

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

Interpreting Chi-Squared
~~~~~~~~~~~~~~~~~~~~~~~~

- **Reduced χ² ≈ 1**: Good fit
- **Reduced χ² >> 1**: Poor fit or underestimated uncertainties
- **Reduced χ² << 1**: Overestimated uncertainties or overfitting

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
