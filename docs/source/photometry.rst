==========
Photometry
==========

The ``Photometry`` class provides methods for querying photometric data from VizieR and applying dust extinction corrections.

Querying Photometric Data
-------------------------

The ``Photometry.query()`` method retrieves photometric measurements from the `VizieR SED service <https://vizier.cds.unistra.fr/>`_.

Query by Source Name
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from starfused import Photometry

   # Query using a source name (Gaia, 2MASS, SIMBAD names work)
   phot = Photometry.query(
       name="Gaia DR3 573956069112683392",
       filters=['GALEX', 'UVOT', 'SDSS', 'Gaia', '2MASS', 'WISE'],
       radius=1  # arcseconds
   )

Query by Coordinates
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Query using RA/Dec coordinates
   phot = Photometry.query(
       ra=150.0,    # degrees
       dec=30.0,    # degrees
       filters=['GALEX', 'SDSS', '2MASS'],
       radius=2     # arcseconds
   )

Supported Filters
~~~~~~~~~~~~~~~~~

The ``filters`` parameter accepts a list of instrument/survey names. Available options include:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Filter Name
     - Bands
     - Wavelength Range
   * - ``GALEX``
     - FUV, NUV
     - 1350-2800 Å (UV)
   * - ``UVOT``
     - UVW2, UVM2, UVW1, U, B, V
     - 1900-5500 Å (UV/Optical)
   * - ``SDSS``
     - u, g, r, i, z
     - 3000-11000 Å (Optical)
   * - ``Gaia``
     - G, BP, RP
     - 3300-10500 Å (Optical)
   * - ``2MASS``
     - J, H, Ks
     - 1.2-2.2 μm (Near-IR)
   * - ``WISE``
     - W1, W2, W3, W4
     - 3.4-22 μm (Mid-IR)

Output Format
~~~~~~~~~~~~~

The query returns a pandas DataFrame with the following columns:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``sed_filter``
     - Filter name (e.g., "SDSS:r", "2MASS:J")
   * - ``sed_wave``
     - Effective wavelength in Angstroms
   * - ``sed_freq``
     - Frequency in GHz
   * - ``sed_flux``
     - Flux in erg/s/cm²/Å
   * - ``sed_eflux``
     - Flux uncertainty in erg/s/cm²/Å
   * - ``_RAJ2000``
     - Right ascension (J2000)
   * - ``_DEJ2000``
     - Declination (J2000)

.. note::

   - Flux is automatically converted from Jy to erg/s/cm²/Å
   - Duplicate filter entries are resolved by keeping the measurement with the lowest uncertainty
   - Measurements with zero uncertainty are handled by computing the standard deviation of multiple measurements if available

Dust Extinction Correction
--------------------------

The ``Photometry.dust_correction()`` method corrects photometry for interstellar dust extinction.

.. code-block:: python

   # Apply Milky Way dust correction
   phot_corrected = Photometry.dust_correction(
       phot,
       extinction='mw',
       dustmap='SandF'
   )

Extinction Curves
~~~~~~~~~~~~~~~~~

Available extinction curves:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Option
     - Description
   * - ``mw`` or ``milky_way``
     - Milky Way diffuse ISM (R_V = 3.1)
   * - ``mw_dense``
     - Milky Way dense regions
   * - ``lmc``
     - Large Magellanic Cloud diffuse
   * - ``lmc_30dor``
     - LMC 30 Doradus region
   * - ``smc_bar``
     - Small Magellanic Cloud bar

Dust Maps
~~~~~~~~~

Available dust maps for E(B-V) retrieval:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Option
     - Description
   * - ``SFD``
     - Schlegel, Finkbeiner & Davis (1998)
   * - ``SandF``
     - Schlafly & Finkbeiner (2011) - recommended

How It Works
~~~~~~~~~~~~

1. The source coordinates are extracted from the photometry DataFrame
2. E(B-V) is queried from the `IRSA Dust Service <https://irsa.ipac.caltech.edu/applications/DUST/>`_
3. The extinction curve is downloaded from STScI
4. A_λ (extinction at each wavelength) is computed: A_λ = E(B-V) × (A_V/E(B-V))
5. Fluxes are corrected: F_corrected = F_observed × 10^(0.4 × A_λ)

Example: Complete Photometry Workflow
-------------------------------------

.. code-block:: python

   from starfused import Photometry
   import matplotlib.pyplot as plt

   # Query photometry
   phot = Photometry.query(
       name="HD 209458",
       filters=['GALEX', 'SDSS', 'Gaia', '2MASS', 'WISE']
   )

   # Apply dust correction
   phot_corrected = Photometry.dust_correction(phot, extinction='mw', dustmap='SandF')

   # Compare original vs corrected
   plt.figure(figsize=(10, 6))
   plt.errorbar(phot['sed_wave'], phot['sed_flux'],
                yerr=phot['sed_eflux'], fmt='o', label='Original', alpha=0.7)
   plt.errorbar(phot_corrected['sed_wave'], phot_corrected['sed_flux'],
                yerr=phot_corrected['sed_eflux'], fmt='s', label='Dust-corrected', alpha=0.7)
   plt.xlabel('Wavelength (Å)')
   plt.ylabel('Flux (erg/s/cm²/Å)')
   plt.xscale('log')
   plt.yscale('log')
   plt.legend()
   plt.title('Effect of Dust Correction')
   plt.show()

API Reference
-------------

.. autoclass:: starfused.Photometry
   :members:
   :undoc-members:
   :show-inheritance:
