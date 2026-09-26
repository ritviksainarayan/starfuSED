===========
SED Fitting
===========

*starfuSED* provides two SED fitting classes: ``SingleSEDFitter`` for single stars and ``BinarySEDFitter`` for binary systems.

How SED Fitting Works
---------------------

The fitting process is a chi-squared grid search over effective temperature and surface gravity at a fixed metallicity:

1. **Grid Search**: For each combination of Teff and log g in the parameter grid:

   a. Load the nearest model spectrum from the native atmosphere grid (models are never interpolated; see :ref:`grid-precision`).
   b. Sample the model at the effective wavelength of each observed band.
   c. Scale the model so that it matches the observed flux in the normalization band.
   d. Calculate chi-squared: χ² = Σ[(F_obs - F_model)² / σ²].

2. **Best Fit**: Select the parameters with the minimum chi-squared.

3. **Radius Calculation**: From the normalization factor (R/d)²:

   R = sqrt(norm) × d

.. note::

   Model fluxes are evaluated at each band's effective wavelength rather than integrated over the filter transmission curve.
   This approximation is less accurate for broad filters and for cool stars with strong molecular features.

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

The ``source_params`` dictionary accepts the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Key
     - Required
     - Description
   * - ``modelname``
     - Yes
     - Model grid: ``'ck04'``, ``'phoenix'``, ``'koester'``, ``'bt-settl'``.
   * - ``teff_min``
     - Yes
     - Minimum temperature to search (K).
   * - ``teff_max``
     - Yes
     - Maximum temperature to search (K).
   * - ``teff_step``
     - No
     - Temperature step size (default: 250 K).
   * - ``logg_min``
     - Yes
     - Minimum log g to search.
   * - ``logg_max``
     - Yes
     - Maximum log g to search.
   * - ``logg_step``
     - No
     - log g step size (default: 0.5).
   * - ``metallicity``
     - Yes
     - Fixed metallicity [M/H] (use ``None`` for Koester).
   * - ``norm_band``
     - Yes
     - Filter name for normalization (e.g., ``'SDSS:r'``, ``'2MASS:H'``). Matched as a case-insensitive substring of ``sed_filter``; the first match is used.

.. tip::

   Choose ``teff_min``, ``logg_min`` and the step sizes so that the search lands on the native nodes of the model grid (see :ref:`grid-precision`).
   Off-node values are replaced by the nearest available model, but the result reports the value that was requested.

Adaptive Fitting
~~~~~~~~~~~~~~~~

To search a wide parameter range more quickly, use adaptive grid refinement:

.. code-block:: python

   result = fitter.fit_adaptive(
       n_refine=2,        # Number of refinement iterations
       refine_factor=4,   # Step size reduction per iteration
       coarse_factor=4    # Initial coarsening factor
   )

The adaptive method:

1. Starts with a coarse grid (steps × ``coarse_factor``).
2. Finds the approximate best fit.
3. Narrows the grid around the best fit and reduces the steps by ``refine_factor``, but never below the ``teff_step`` and ``logg_step`` you specified.
4. Repeats ``n_refine`` times.

Adaptive fitting reduces the number of models evaluated; it does not resolve parameters more finely than your step sizes or the native model grid.

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

Both components share the supplied distance. The binary fitter scales them as follows:

- **Source 1**: Anchored to its ``norm_band`` (typically a UV band for the hot component).
- **Source 2**: In ``fit()``, anchored so that the combined model matches the observed flux in its ``norm_band``, or, if ``'norm_band': None``, scaled by weighted least squares to the residual flux after subtracting source 1. ``fit_adaptive()`` always uses the least-squares scaling for source 2.

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
       'n_data': 15,
       'n_params': 4,
       'combined_flux': array,
       'distance_pc': 100
   }

.. _grid-precision:

Precision and Grid Spacing
--------------------------

As described in Narayan & Soares-Furtado (2026), *starfuSED* compares the photometry only with model spectra at the native nodes of each atmosphere grid and never interpolates between neighboring models.
The fitted Teff and log g are therefore discrete: they cannot be determined more finely than the local spacing of the grid near the solution, which is typically 100–250 K in Teff and 0.25–0.5 dex in log g.
Finer ``teff_step`` or ``logg_step`` values and adaptive refinement do not change this.

We recommend reporting uncertainties on fitted parameters that are no smaller than one grid step of the model used.

Native grid spacing:

.. list-table::
   :header-rows: 1
   :widths: 15 55 30

   * - Grid
     - Teff spacing
     - log g spacing
   * - ``ck04``
     - 250 K (3500–13000 K); 1000 K (13000–50000 K).
     - 0.5 dex.
   * - ``phoenix``
     - 100 K (2000–7000 K); 200 K (7000–12000 K); 500 K (12000–20000 K); 1000 K (20000–70000 K).
     - 0.5 dex.
   * - ``koester``
     - 250 K (5000–40000 K); 1000 K (40000–80000 K).
     - 0.25 dex.
   * - ``bt-settl``
     - 100 K (400–7000 K).
     - 0.5 dex.

Convenience Functions
---------------------

For quick one-liner fitting with the standard grid search:

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

- Choose a band where the star dominates (not contaminated by a binary companion).
- Avoid bands with large photometric errors.
- For binaries, pick bands where each component clearly dominates.

.. code-block:: python

   # Good: WD dominates in UV, companion dominates in IR
   source1 = {..., 'norm_band': 'GALEX:NUV'}  # WD
   source2 = {..., 'norm_band': '2MASS:H'}    # M dwarf

   # Bad: Both contribute at optical wavelengths
   source1 = {..., 'norm_band': 'SDSS:r'}  # Both components contribute
   source2 = {..., 'norm_band': 'SDSS:i'}  # Both components contribute

Setting Parameter Ranges
~~~~~~~~~~~~~~~~~~~~~~~~

- Start with wide ranges, then narrow them based on initial fits.
- Use adaptive fitting to efficiently search large parameter spaces.
- Check that the best fit is not at the edge of the grid (this suggests the range is wrong).

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
