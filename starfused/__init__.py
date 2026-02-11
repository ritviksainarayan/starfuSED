"""
starfused - SED fitting.

"""

from .preprocess import Photometry
from .spectra import StellarModel, BaseModel, CKModel, PhoenixModel, KoesterModel, BTSettlModel
from .fitter import fit_binary_sed, fit_single_sed, SEDFitter, SingleSEDFitter, BinarySEDFitter
from .plotter import plot_sed, plot_single_sed, plot_binary_sed, plot_mc_ridgeline

__version__ = "0.1.0"

__all__ = [
    'Photometry',
    'StellarModel',
    'BaseModel',
    'CKModel',
    'PhoenixModel',
    'KoesterModel',
    'BTSettlModel',
    'fit_binary_sed',
    'fit_single_sed',
    'SEDFitter',
    'SingleSEDFitter',
    'BinarySEDFitter',
    'plot_sed',
    'plot_single_sed',
    'plot_binary_sed',
    'plot_mc_ridgeline',
]