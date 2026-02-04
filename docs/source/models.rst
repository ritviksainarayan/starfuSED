==============
Stellar Models
==============

*starfuSED* provides access to four different stellar atmosphere model grids through a unified interface. Models are downloaded on-demand and cached in memory during your session.

Available Model Grids
---------------------

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 20 25

   * - Grid
     - Teff Range
     - log g Range
     - [M/H] Range
     - Best For
   * - ``ck04``
     - 3500-50000 K
     - 0.0-5.0
     - -2.5 to +0.5
     - Main sequence stars
   * - ``phoenix``
     - Full range
     - 0.0-5.0
     - -4.0 to +0.5
     - Cool stars, M dwarfs
   * - ``koester``
     - 5000-80000 K
     - 6.5-9.5
     - N/A
     - DA white dwarfs
   * - ``bt-settl``
     - 400-70000 K
     - -0.5 to 6.0
     - -4.0 to +0.5
     - Brown dwarfs, L/T/Y dwarfs

Using the StellarModel Interface
--------------------------------

The ``StellarModel`` class provides a unified interface to all model grids:

.. code-block:: python

   from starfused import StellarModel

   # Initialize with desired grid
   model = StellarModel(grid='ck04')

   # Load a spectrum
   spectrum = model.load_model(
       teff=5750,        # Effective temperature (K)
       logg=4.5,         # Surface gravity
       metallicity=0.0   # [M/H] in dex
   )

   # spectrum is a DataFrame with 'wavelength' and 'flux' columns
   print(spectrum.head())

Finding the Closest Model
~~~~~~~~~~~~~~~~~~~~~~~~~

Model grids have discrete parameter values. Use ``find_model()`` to see what parameters will actually be used:

.. code-block:: python

   model = StellarModel(grid='ck04')

   # Find closest available model
   params = model.find_model(teff=5800, logg=4.3, metallicity=0.1)
   print(params)
   # Output: {'teff': 5750, 'logg': 4.5, 'metallicity': 0.0, ...}

Castelli-Kurucz 2004 (CK04)
---------------------------

The CK04 grid uses ATLAS9 model atmospheres and is the most commonly used for main-sequence stars.

**Source**: `STScI CDBS Archive <https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/>`_

**Parameter Grid**:

- **Teff**: 3500-50000 K in 250 K steps
- **log g**: 0.0-5.0 in 0.5 dex steps
- **[M/H]**: -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5

.. code-block:: python

   model = StellarModel(grid='ck04')

   # Solar-type star
   solar = model.load_model(teff=5750, logg=4.5, metallicity=0.0)

   # Metal-poor giant
   giant = model.load_model(teff=4500, logg=2.0, metallicity=-1.5)

   # Hot star
   hot_star = model.load_model(teff=15000, logg=4.0, metallicity=0.0)

PHOENIX Models
--------------

PHOENIX models are computed with a more detailed treatment of molecular opacities, making them better for cool stars.

**Source**: `STScI CDBS Archive <https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/>`_

**Parameter Grid**:

- **Teff**: Full range with ~100 K resolution
- **log g**: 0.0-5.0 in 0.5 dex steps
- **[M/H]**: -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, 0.0, +0.3, +0.5

.. code-block:: python

   model = StellarModel(grid='phoenix')

   # M dwarf
   m_dwarf = model.load_model(teff=3500, logg=5.0, metallicity=0.0)

   # K giant
   k_giant = model.load_model(teff=4200, logg=2.5, metallicity=-0.5)

Koester White Dwarf Models
--------------------------

Koester models are pure hydrogen (DA) white dwarf atmospheres, ideal for fitting hot compact objects.

**Source**: `Spanish Virtual Observatory <http://svo2.cab.inta-csic.es/theory/newov2/>`_

**Parameter Grid**:

- **Teff**: 5000-40000 K (250 K steps), 40000-80000 K (1000 K steps)
- **log g**: 6.5-9.5 in 0.25 dex steps
- **[M/H]**: Not applicable (pure H atmosphere)

.. code-block:: python

   model = StellarModel(grid='koester')

   # Typical white dwarf
   wd = model.load_model(teff=12000, logg=8.0, metallicity=None)

   # Hot white dwarf
   hot_wd = model.load_model(teff=40000, logg=7.5, metallicity=None)

   # Cool white dwarf
   cool_wd = model.load_model(teff=6000, logg=8.5, metallicity=None)

BT-Settl Models
---------------

BT-Settl models include cloud formation and dust settling physics, making them ideal for very cool objects like brown dwarfs.

**Source**: `Spanish Virtual Observatory <http://svo2.cab.inta-csic.es/theory/newov2/>`_

**Reference**: Allard, Homeier & Freytag (2012)

**Parameter Grid**:

- **Teff**: 400-70000 K (100 K steps for T < 7000 K)
- **log g**: -0.5 to 6.0 in 0.5 dex steps
- **[M/H]**: -4.0 to +0.5

.. code-block:: python

   model = StellarModel(grid='bt-settl')

   # L dwarf
   l_dwarf = model.load_model(teff=1800, logg=5.0, metallicity=0.0)

   # T dwarf
   t_dwarf = model.load_model(teff=1200, logg=5.0, metallicity=0.0)

   # Very low mass star
   vlm_star = model.load_model(teff=2800, logg=5.0, metallicity=0.0)

Directly Accessing Model Classes
--------------------------------

For advanced usage, you can directly instantiate the individual model classes:

.. code-block:: python

   from starfused import CKModel, PhoenixModel, KoesterModel, BTSettlModel

   # Direct access to CK04 models
   ck = CKModel()
   spectrum = ck.load_model(teff=6000, logg=4.0, metallicity=0.0)

   # Validate parameters before loading
   ck.validate_parameters(teff=6000, logg=4.0, metallicity=0.0)  # Returns True

   # Construct URL without downloading
   url = ck.construct_url(teff=6000, logg=4.0, metallicity=0.0)
   print(url)

Plotting Model Spectra
----------------------

.. code-block:: python

   import matplotlib.pyplot as plt
   from starfused import StellarModel

   # Compare different stellar types
   fig, ax = plt.subplots(figsize=(12, 6))

   models_to_plot = [
       ('ck04', 10000, 4.0, 0.0, 'A star (10000 K)'),
       ('ck04', 5750, 4.5, 0.0, 'G star (5750 K)'),
       ('phoenix', 3500, 4.5, 0.0, 'M star (3500 K)'),
       ('koester', 15000, 8.0, None, 'White dwarf (15000 K)'),
   ]

   for grid, teff, logg, met, label in models_to_plot:
       model = StellarModel(grid=grid)
       spec = model.load_model(teff=teff, logg=logg, metallicity=met, verbose=False)
       ax.plot(spec['wavelength'], spec['flux'], label=label, alpha=0.8)

   ax.set_xlabel('Wavelength (Å)')
   ax.set_ylabel('Flux (erg/s/cm²/Å)')
   ax.set_xlim(1000, 25000)
   ax.set_yscale('log')
   ax.legend()
   ax.set_title('Comparison of Stellar Model Spectra')
   plt.show()

Choosing the Right Model Grid
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Use Case
     - Recommended Grid
   * - Main sequence F/G/K stars
     - ``ck04`` (fast, well-tested)
   * - M dwarfs (Teff < 4000 K)
     - ``phoenix`` (better molecular opacities)
   * - DA white dwarfs
     - ``koester`` (purpose-built)
   * - Brown dwarfs / L/T/Y dwarfs
     - ``bt-settl`` (includes cloud physics)
   * - Very metal-poor stars
     - ``phoenix`` (extends to [M/H] = -4)
   * - Binary with WD + cool companion
     - ``koester`` + ``bt-settl``

API Reference
-------------

.. autoclass:: starfused.StellarModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.BaseModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.CKModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.PhoenixModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.KoesterModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: starfused.BTSettlModel
   :members:
   :undoc-members:
   :show-inheritance:
