.. image:: starfused.png
   :align: center
   :width: 400px
   :alt: starfuSED

========
starfuSED
========

**Spectral Energy Distribution fitting for stars**

*starfuSED* is a Python package for SED (Spectral Energy Distribution) fitting of single and binary stellar systems. It provides a complete workflow from photometry retrieval to model fitting and visualization, all with in-memory processing.

.. note::

   This project is under active development.

Key Features
------------

- **In-Memory Processing**: No need to download large dust maps or spectral model grids locally.
- **Multiple Stellar Model Grids**: Unified interface to four atmospheric model grids:

  - Castelli-Kurucz 2004 models.
  - PHOENIX stellar atmosphere models.
  - Koester DA white dwarf models.
  - BT-Settl models for cool stars and brown dwarfs.

- **Automatic Photometry Retrieval**: Query VizieR SED service by object name or coordinates.
- **Dust Extinction Correction**: Multiple extinction curves and dust maps supported.
- **Single and Binary SED Fitting**: Chi-squared minimization with adaptive grid refinement.
- **Publication-Ready Plotting**: Customizable SED plots with residual panels.

Quick Example
-------------

.. code-block:: python

   from starfused import Photometry, SingleSEDFitter, plot_single_sed

   # Query photometry for a star
   phot = Photometry.query(
       name="Gaia DR3 573956069112683392",
       filters=['GALEX', 'SDSS', 'Gaia', '2MASS', 'WISE']
   )

   # Apply dust correction
   phot_corrected = Photometry.dust_correction(phot, extinction='mw', dustmap='SandF')

   # Define fitting parameters
   params = {
       'modelname': 'ck04',
       'teff_min': 5000, 'teff_max': 8000,
       'logg_min': 3.5, 'logg_max': 5.0,
       'metallicity': 0.0,
       'norm_band': 'SDSS:r'
   }

   # Fit the SED
   fitter = SingleSEDFitter(phot_corrected, distance_pc=50, source_params=params)
   result = fitter.fit_adaptive()

   # Plot the result
   fig, axes = plot_single_sed(result, phot_corrected)

Documentation
-------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   photometry
   models
   fitting
   plotting

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api

Acknowledgments
---------------

*starfuSED* relies on data and services from:

- `VizieR SED Service <https://vizier.cds.unistra.fr/>`_ for photometric data.
- `STScI CDBS Archive <https://archive.stsci.edu/>`_ for Castelli-Kurucz and PHOENIX models.
- `Spanish Virtual Observatory <https://svo2.cab.inta-csic.es/theory/newov2/index.php>`_ for Koester and BT-Settl models.
- `IRSA Dust Service <https://irsa.ipac.caltech.edu/applications/DUST/>`_ for extinction values.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
