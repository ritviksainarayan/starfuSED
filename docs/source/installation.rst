============
Installation
============

Requirements
------------

*starfuSED* requires Python 3.8 or later. The following dependencies will be installed automatically:

**Core Dependencies**

- ``numpy`` - Numerical computing
- ``pandas`` - Data manipulation and DataFrames
- ``scipy`` - Scientific computing (interpolation)
- ``matplotlib`` - Plotting and visualization
- ``requests`` - HTTP requests for data retrieval

**Astronomy Dependencies**

- ``astropy`` - Units, coordinates, FITS handling
- ``astroquery`` - Access to astronomical databases (IRSA Dust)

Installing from Source
----------------------

Clone the repository and install in development mode:

.. code-block:: bash

   git clone https://github.com/ritviksainarayan/starfuSED.git
   cd starfuSED
   pip install -e .

Installing Dependencies Only
----------------------------

If you want to install just the dependencies first:

.. code-block:: bash

   pip install numpy pandas scipy matplotlib requests astropy astroquery

Verifying Installation
----------------------

After installation, verify that *starfuSED* is working correctly:

.. code-block:: python

   import starfused
   print(starfused.__version__)

   # Test a simple query
   from starfused import Photometry, StellarModel

   # This should print the version and work without errors
   model = StellarModel(grid='ck04')
   print("Installation successful!")

Network Requirements
--------------------

*starfuSED* requires an internet connection to:

1. Query photometry from VizieR SED service
2. Download stellar atmosphere models from STScI and SVO
3. Retrieve dust extinction values from IRSA

All data is downloaded on-demand and cached in memory during a session. No persistent local storage of large files is required.

Troubleshooting
---------------

**ImportError: No module named 'astroquery'**
    Install astroquery: ``pip install astroquery``

**Connection errors when querying data**
    Ensure you have an active internet connection. Some institutional firewalls may block access to astronomical data services.

**Memory errors with large model grids**
    When fitting with fine parameter grids, memory usage can be high. Consider using the ``fit_adaptive()`` method which starts with a coarse grid and refines progressively.
