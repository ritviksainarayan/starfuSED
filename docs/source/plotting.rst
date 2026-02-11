========
Plotting
========

*starfuSED* provides publication-ready plotting functions for visualizing SED fits.

Basic Plotting
--------------

Single Star SED
~~~~~~~~~~~~~~~

.. code-block:: python

   from starfused import plot_single_sed

   fig, axes = plot_single_sed(result, photometry_df)

Binary System SED
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from starfused import plot_binary_sed

   fig, axes = plot_binary_sed(result, photometry_df)

Auto-Detection
~~~~~~~~~~~~~~

Use ``plot_sed()`` to automatically detect the fit type:

.. code-block:: python

   from starfused import plot_sed

   # Automatically calls plot_single_sed or plot_binary_sed
   fig, axes = plot_sed(result, photometry_df)

Customization Options
---------------------

Figure and Axes
~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_single_sed(
       result, phot,
       figsize=(14, 9),           # Figure size (width, height)
       xlim=(1000, 100000),       # Wavelength range in Angstroms
       ylim=(1e-18, 1e-14),       # Flux range
       title="Custom Title"       # Override auto-generated title
   )

Residual Panel
~~~~~~~~~~~~~~

.. code-block:: python

   # Show residuals
   fig, (ax_main, ax_resid) = plot_single_sed(
       result, phot,
       show_residuals=True,
       residual_ylim=(-5, 5),      # Y-limits for residuals
       residual_mode='sigma',      # 'sigma' or 'absolute'
       residual_ylabel=r'$\Delta$ ($\sigma$)'
   )

   # Hide residuals
   fig, ax = plot_single_sed(result, phot, show_residuals=False)

Residual modes:

- ``'sigma'``: (observed - model) / uncertainty — normalized residuals
- ``'absolute'``: observed - model — raw flux difference

Data Point Styling
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_single_sed(
       result, phot,
       show_filter_labels=True,    # Annotate points with filter names
       data_kwargs={
           'color': 'navy',
           'markersize': 10,
           'capsize': 4,
           'fmt': 's',             # Square markers
           'label': 'Observations'
       }
   )

Model Styling
~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_single_sed(
       result, phot,
       fill_under=True,            # Fill area under model curve
       fill_alpha=0.2,             # Transparency of fill
       model_kwargs={
           'color': 'crimson',
           'linewidth': 2,
           'alpha': 0.8,
           'label': 'Best-fit Model'
       },
       fill_kwargs={
           'color': 'crimson',
           'alpha': 0.15
       }
   )

Binary-Specific Options
-----------------------

Component Visibility
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_binary_sed(
       result, phot,
       show_components=True,       # Plot individual source spectra
       fill_under=True             # Fill under component curves
   )

Component Colors and Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_binary_sed(
       result, phot,
       source1_label="White Dwarf (15000 K)",
       source2_label="M Dwarf (3200 K)",
       source1_color='blue',
       source2_color='orange',
       combined_color='red',
       source1_kwargs={'linewidth': 2, 'linestyle': '--'},
       source2_kwargs={'linewidth': 2, 'linestyle': ':'},
       combined_kwargs={'linewidth': 3}
   )

Component Fills
~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_binary_sed(
       result, phot,
       fill_under=True,
       source1_fill_kwargs={'alpha': 0.2, 'hatch': '//'},
       source2_fill_kwargs={'alpha': 0.2, 'hatch': '\\\\'}
   )

Using Existing Axes
-------------------

For multi-panel figures or complex layouts:

.. code-block:: python

   import matplotlib.pyplot as plt

   # Create your own figure
   fig = plt.figure(figsize=(16, 10))

   # Without residuals - pass single axes
   ax1 = fig.add_subplot(121)
   plot_single_sed(result1, phot1, ax=ax1, fig=fig, show_residuals=False)

   # With residuals - pass tuple of axes
   ax2 = fig.add_subplot(222)
   ax3 = fig.add_subplot(224)
   plot_single_sed(result2, phot2, ax=(ax2, ax3), fig=fig, show_residuals=True)

   plt.tight_layout()
   plt.savefig('comparison.pdf')

Example Gallery
---------------

Publication-Ready Single Star
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_single_sed(
       result, phot,
       figsize=(10, 7),
       xlim=(1500, 50000),
       show_filter_labels=True,
       show_residuals=True,
       residual_ylim=(-3, 3),
       fill_under=True,
       fill_alpha=0.15,
       data_kwargs={
           'color': 'black',
           'markersize': 8,
           'capsize': 3,
           'elinewidth': 1.5,
           'markeredgewidth': 1.5
       },
       model_kwargs={
           'color': '#E63946',
           'linewidth': 1.5
       }
   )

   # Customize axes labels
   axes[0].set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
   axes[1].set_xlabel(r'Wavelength (Å)', fontsize=12)

   plt.savefig('sed_publication.pdf', dpi=300, bbox_inches='tight')

WD + M Dwarf Binary
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fig, axes = plot_binary_sed(
       result, phot,
       figsize=(12, 8),
       xlim=(1500, 50000),
       show_components=True,
       show_filter_labels=True,
       fill_under=True,
       source1_label=f"WD ({result['source1']['teff']} K)",
       source2_label=f"dM ({result['source2']['teff']} K)",
       source1_color='#457B9D',
       source2_color='#E9C46A',
       combined_color='#2A9D8F'
   )

Comparing Multiple Fits
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(1, 2, figsize=(16, 6))

   # CK04 fit
   plot_single_sed(
       result_ck04, phot,
       ax=axes[0], fig=fig,
       show_residuals=False,
       title="CK04 Model Fit"
   )

   # Phoenix fit
   plot_single_sed(
       result_phoenix, phot,
       ax=axes[1], fig=fig,
       show_residuals=False,
       title="Phoenix Model Fit"
   )

   plt.tight_layout()

Monte Carlo Ridgeline Plot
--------------------------

After running ``fit_mc()``, visualize the parameter distributions with a ridgeline (joy) plot:

.. code-block:: python

   from starfused import plot_mc_ridgeline

   # Single star
   fig, axes = plot_mc_ridgeline(result)

   # Binary — auto-detects and shows both sources
   fig, axes = plot_mc_ridgeline(binary_result)

   # Custom parameters
   fig, axes = plot_mc_ridgeline(result, params=['teff', 'radius_rjup'])

   # Enable optional markers
   fig, axes = plot_mc_ridgeline(result, show_best=True, show_percentiles=True)

Each ridge is a KDE-smoothed density of the MC samples for one parameter.
The function auto-detects single vs binary fits and selects appropriate defaults.

API Reference
-------------

.. autofunction:: starfused.plot_single_sed

.. autofunction:: starfused.plot_binary_sed

.. autofunction:: starfused.plot_sed

.. autofunction:: starfused.plot_mc_ridgeline
