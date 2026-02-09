==========
Quickstart
==========

This guide will walk you through the basic workflow of *starfuSED*: querying photometry, fitting an SED model, and visualizing the results.

Basic Workflow
--------------

The typical workflow involves four steps:

1. **Query photometry** - Retrieve photometric measurements from VizieR
2. **Apply dust correction** - Correct for interstellar extinction
3. **Fit the SED** - Find the best-fitting stellar model
4. **Visualize results** - Create publication-ready plots

Single Star SED Fitting
-----------------------

Here's a complete example of fitting a single star:

.. code-block:: python

   from starfused import Photometry, SingleSEDFitter, plot_single_sed
   import matplotlib.pyplot as plt

   # Step 1: Query photometry by source name
   phot = Photometry.query(
       name="Gaia DR3 573956069112683392",
       filters=['GALEX', 'UVOT', 'SDSS', 'Gaia', '2MASS', 'WISE'],
       radius=1  # search radius in arcseconds
   )

   print(f"Retrieved {len(phot)} photometric measurements")
   print(phot[['sed_filter', 'sed_wave', 'sed_flux', 'sed_eflux']])

   # Step 2: Apply dust extinction correction
   phot_corrected = Photometry.dust_correction(
       phot,
       extinction='mw',      # Milky Way extinction curve
       dustmap='SandF'       # Schlafly & Finkbeiner dust map
   )

   # Step 3: Define model parameters and fit
   params = {
       'modelname': 'ck04',          # Castelli-Kurucz 2004 models
       'teff_min': 5000,             # Minimum temperature to search
       'teff_max': 8000,             # Maximum temperature to search
       'teff_step': 250,             # Temperature step size
       'logg_min': 3.5,              # Minimum log g
       'logg_max': 5.0,              # Maximum log g
       'logg_step': 0.5,             # log g step size
       'metallicity': 0.0,           # Solar metallicity
       'norm_band': 'SDSS:r'         # Normalize model to this band
   }

   fitter = SingleSEDFitter(
       phot_corrected,
       distance_pc=50,               # Distance to the star
       source_params=params
   )

   # Use adaptive fitting for better precision
   result = fitter.fit_adaptive(n_refine=2)

   # Estimate uncertainties via Monte Carlo perturbation
   result = fitter.fit_mc(n_iter=100, sigma_clip=3, seed=42)

   # Step 4: Plot the result
   fig, axes = plot_single_sed(
       result,
       phot_corrected,
       show_residuals=True,
       show_filter_labels=True
   )
   plt.savefig('sed_fit.png', dpi=150, bbox_inches='tight')
   plt.show()

Binary System SED Fitting
-------------------------

For binary systems (e.g., white dwarf + M dwarf), use ``BinarySEDFitter``:

.. code-block:: python

   from starfused import Photometry, BinarySEDFitter, plot_binary_sed

   # Query photometry
   phot = Photometry.query(name="TIC 12345678", filters=['GALEX', 'SDSS', '2MASS', 'WISE'])
   phot_corrected = Photometry.dust_correction(phot)

   # Define parameters for both sources
   source1_params = {
       'modelname': 'koester',       # White dwarf models
       'teff_min': 10000,
       'teff_max': 25000,
       'teff_step': 250,
       'logg_min': 7.5,
       'logg_max': 8.5,
       'logg_step': 0.25,
       'metallicity': None,          # WD models don't use metallicity
       'norm_band': 'GALEX:NUV'      # WD dominates in UV
   }

   source2_params = {
       'modelname': 'bt-settl',      # Cool star/brown dwarf models
       'teff_min': 2500,
       'teff_max': 4000,
       'teff_step': 100,
       'logg_min': 4.5,
       'logg_max': 5.5,
       'logg_step': 0.5,
       'metallicity': 0.0,
       'norm_band': '2MASS:H'        # Cool companion dominates in IR
   }

   # Fit the binary
   fitter = BinarySEDFitter(
       phot_corrected,
       distance_pc=100,
       source1_params=source1_params,
       source2_params=source2_params
   )
   result = fitter.fit_adaptive()

   # Plot with component separation
   fig, axes = plot_binary_sed(
       result,
       phot_corrected,
       show_components=True,
       fill_under=True
   )

Querying by Coordinates
-----------------------

You can also query photometry using coordinates instead of a source name:

.. code-block:: python

   phot = Photometry.query(
       ra=150.0,                     # Right ascension in degrees
       dec=30.0,                     # Declination in degrees
       filters=['GALEX', 'SDSS', '2MASS'],
       radius=2                      # Larger radius for crowded fields
   )

Working with Different Model Grids
----------------------------------

*starfuSED* supports four stellar atmosphere model grids:

.. code-block:: python

   from starfused import StellarModel

   # Castelli-Kurucz 2004 (main sequence stars)
   ck_model = StellarModel(grid='ck04')
   spectrum = ck_model.load_model(teff=5750, logg=4.5, metallicity=0.0)

   # PHOENIX (main sequence, better for cool stars)
   phoenix_model = StellarModel(grid='phoenix')
   spectrum = phoenix_model.load_model(teff=4000, logg=4.5, metallicity=0.0)

   # Koester (white dwarfs)
   koester_model = StellarModel(grid='koester')
   spectrum = koester_model.load_model(teff=15000, logg=8.0, metallicity=None)

   # BT-Settl (cool stars, brown dwarfs)
   btsettl_model = StellarModel(grid='bt-settl')
   spectrum = btsettl_model.load_model(teff=2500, logg=5.0, metallicity=0.0)

Understanding the Fit Results
-----------------------------

The fit result dictionary contains:

.. code-block:: python

   result = fitter.fit_adaptive()

   # Best-fit parameters
   print(f"Temperature: {result['teff']} K")
   print(f"Surface gravity: {result['logg']}")
   print(f"Metallicity: {result['metallicity']}")

   # Inferred physical properties
   print(f"Radius: {result['radius_rsun']:.4f} R_sun")
   print(f"Radius: {result['radius_rearth']:.2f} R_earth")
   print(f"Radius: {result['radius_rjup']:.4f} R_jup")

   # Fit quality metrics
   print(f"Chi-squared: {result['chi2']:.2f}")
   print(f"Reduced chi-squared: {result['reduced_chi2']:.2f}")

   # Monte Carlo uncertainties (after calling fit_mc())
   if 'mc_errors' in result:
       err = result['mc_errors']
       print(f"Teff = {result['teff']} (+{err['teff_err'][1]:.0f}/-{err['teff_err'][0]:.0f}) K")
       print(f"Radius = {result['radius_rsun']:.4f} "
             f"(+{err['radius_rsun_err'][1]:.4f}/-{err['radius_rsun_err'][0]:.4f}) R_sun")

   # Best-fit model spectrum
   spectrum = result['spectrum']
   print(spectrum.head())

Customizing Plots
-----------------

The plotting functions offer extensive customization:

.. code-block:: python

   fig, axes = plot_single_sed(
       result,
       phot_corrected,
       title="My Custom Title",
       figsize=(12, 8),
       xlim=(1000, 50000),           # Wavelength range in Angstroms
       ylim=(1e-18, 1e-14),          # Flux range
       show_filter_labels=True,
       show_residuals=True,
       residual_ylim=(-3, 3),
       fill_under=True,
       fill_alpha=0.2,
       data_kwargs={'color': 'navy', 'markersize': 10},
       model_kwargs={'color': 'crimson', 'linewidth': 2}
   )

Next Steps
----------

- See :doc:`photometry` for detailed photometry handling options
- See :doc:`models` for information about available stellar model grids
- See :doc:`fitting` for advanced fitting techniques
- See :doc:`plotting` for all visualization options
- See :doc:`api` for the complete API reference
