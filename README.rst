.. image:: docs/source/starfused.png
   :align: center
   :width: 400px
   :alt: starfuSED

=========
starfuSED
=========

**Spectral Energy Distribution fitting for stars**

*starfuSED* is a Python package for SED (Spectral Energy Distribution) fitting of single and binary stellar systems. It provides a complete workflow from photometry retrieval to model fitting and visualization, all with in-memory processing.

.. image:: https://img.shields.io/badge/status-in%20development-yellow
   :alt: Development Status

Key Features
------------

- **In-Memory Processing**: No need to download large dust maps or spectral model grids locally.
- **Multiple Stellar Model Grids**: Unified interface to four atmospheric model grids:

  - Castelli-Kurucz 2004 modes.
  - PHOENIX stellar atmosphere models.
  - Koester DA white dwarf models.
  - BT-Settl models for cool stars and brown dwarfs.

- **Automatic Photometry Retrieval**: Query VizieR SED service by object name or coordinates. 
- **Dust Extinction Correction**: Multiple extinction curves and dust maps supported. 
- **Single and Binary SED Fitting**: Chi-squared minimization with adaptive grid refinement. 
- **Publication-Ready Plotting**: Customizable SED plots with residual panels. 

Installation
------------

.. code-block:: bash

   git clone https://github.com/ritviksainarayan/starfuSED.git
   cd starfuSED
   pip install -e .

Or install dependencies only:

.. code-block:: bash

   pip install numpy pandas scipy matplotlib requests astropy astroquery

Quick Start
-----------

.. code-block:: python

   from starfused import Photometry, SingleSEDFitter, plot_single_sed

   # Query photometry
   phot = Photometry.query(
       name="Gaia DR3 573956069112683392",
       filters=['GALEX', 'SDSS', 'Gaia', '2MASS', 'WISE']
   )

   # Apply dust correction
   phot = Photometry.dust_correction(phot, extinction='mw', dustmap='SandF')

   # Define fitting parameters
   params = {
       'modelname': 'ck04',
       'teff_min': 5000, 'teff_max': 8000,
       'logg_min': 3.5, 'logg_max': 5.0,
       'metallicity': 0.0,
       'norm_band': 'SDSS:r'
   }

   # Fit and plot
   fitter = SingleSEDFitter(phot, distance_pc=50, source_params=params)
   result = fitter.fit_adaptive()
   fig, axes = plot_single_sed(result, phot)

Documentation
-------------

Full documentation is available at `Read the Docs <https://starfused.readthedocs.io/>`_.

Requirements
------------

- Python >= 3.8
- numpy >= 1.20
- pandas >= 1.3
- scipy >= 1.7
- matplotlib >= 3.4
- astropy >= 5.0
- astroquery >= 0.4.6
- requests >= 2.25

License
-------

MIT License

Author
------

Ritvik Sai Narayan

Acknowledgments
---------------

This package uses data from:

- `VizieR SED Service <https://vizier.cds.unistra.fr/>`_
- `STScI CDBS Archive <https://archive.stsci.edu/>`_
- `Spanish Virtual Observatory <https://svo2.cab.inta-csic.es/theory/newov2/index.php>`_
- `IRSA Dust Service <https://irsa.ipac.caltech.edu/applications/DUST/>`_
