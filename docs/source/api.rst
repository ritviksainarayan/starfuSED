=============
API Reference
=============

This page contains the complete API documentation for *starfuSED*.

Package Contents
----------------

.. currentmodule:: starfused

The ``starfused`` package exports the following:

**Photometry**

.. autosummary::
   :toctree: generated

   Photometry

**Stellar Models**

.. autosummary::
   :toctree: generated

   StellarModel
   BaseModel
   CKModel
   PhoenixModel
   KoesterModel
   BTSettlModel

**SED Fitting**

.. autosummary::
   :toctree: generated

   SEDFitter
   SingleSEDFitter
   BinarySEDFitter
   fit_single_sed
   fit_binary_sed

**Plotting**

.. autosummary::
   :toctree: generated

   plot_sed
   plot_single_sed
   plot_binary_sed
   plot_mc_ridgeline

Photometry Module
-----------------

.. automodule:: starfused.preprocess
   :members:
   :undoc-members:
   :show-inheritance:

Stellar Models Module
---------------------

.. automodule:: starfused.spectra
   :members:
   :undoc-members:
   :show-inheritance:

Fitting Module
--------------

.. automodule:: starfused.fitter
   :members:
   :undoc-members:
   :show-inheritance:

Plotting Module
---------------

.. automodule:: starfused.plotter
   :members:
   :undoc-members:
   :show-inheritance:
