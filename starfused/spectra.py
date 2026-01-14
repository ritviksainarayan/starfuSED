import pandas as pd
import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
import requests
from astropy.io import fits
from io import BytesIO
from scipy.interpolate import interp1d
from abc import ABC, abstractmethod
import re


class BaseModel(ABC):
    """
    Abstract base class for stellar atmosphere model grids.

    Parameters
    ----------
    base_url : str
        Base URL for the model grid repository.
    """

    def __init__(self, base_url):
        self.base_url = base_url

    @abstractmethod
    def validate_parameters(self, teff, logg, metallicity):
        """
        Validate that stellar parameters are within the grid's range.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        bool
            True if parameters are valid.

        Raises
        ------
        ValueError
            If any parameter is out of range for this grid.
        """
        pass

    @abstractmethod
    def find_model(self, teff, logg, metallicity):
        """
        Find the closest available model in the grid.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        dict
            Dictionary containing the closest grid point parameters and file information.
        """
        pass

    @abstractmethod
    def construct_url(self, teff, logg, metallicity):
        """
        Construct the URL for downloading a model spectrum.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        str
            Complete URL to the model file.
        """
        pass

    def load_model(self, teff, logg, metallicity):
        """
        Load a stellar atmosphere model spectrum.

        This method finds the closest model in the grid, downloads it,
        and returns the spectrum as a DataFrame.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' (Angstroms) and 'flux' (erg/s/cm^2/A).

        Raises
        ------
        RuntimeError
            If the model cannot be downloaded or parsed.
        """

        params = self.find_model(teff, logg, metallicity)
        print(f"Loading model: Teff={params['teff']}K, log_g={params['logg']}, [M/H]={params['metallicity']}")

        url = self.construct_url(params['teff'], params['logg'], params['metallicity'])

        try:
            response = requests.get(url)
            response.raise_for_status()

            df = self.parse_fits(BytesIO(response.content), logg)

            return df
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {url}: {e}")

    @abstractmethod
    def parse_fits(self, fitsfile, logg):
        """
        Parse a FITS file containing stellar atmosphere model data.

        Parameters
        ----------
        fitsfile : file-like object
            FITS file content (typically BytesIO object).
        logg : float
            Surface gravity (log g) to extract from the file.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' (Angstroms) and 'flux' (erg/s/cm^2/A).
        """
        pass

class CKModel(BaseModel):
    """
    Castelli-Kurucz 2004 ATLAS9 stellar atmosphere model grid.

    This class provides access to the Castelli & Kurucz (2004) ATLAS9 model atmospheres
    hosted at STScI. The grid covers a wide range of stellar parameters suitable for
    modeling main-sequence through evolved stars.
    """

    METALLICITIES = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]
    LOGG_RANGE = (0.0, 5.0)
    TEFF_RANGE = (3500, 50000)
    TEFF_STEP = 250
    LOGG_STEP = 0.5



    def __init__(self):
        """Initialize CKModel with STScI repository URL."""
        super().__init__("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/")

    def _fehtodir(self, metallicity):
        """
        Convert metallicity to directory name format.

        Parameters
        ----------
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        str
            Directory name (e.g., 'ckp00' for [M/H]=0.0, 'ckm25' for [M/H]=-2.5).
        """

        dirval = int(np.abs(metallicity) * 10)

        if metallicity >= 0:
            return f"ckp{dirval:02d}"
        else:
            return f"ckm{dirval:02d}"
        
    def validate_parameters(self, teff, logg, metallicity):
        """
        Validate stellar parameters against CK04 grid ranges.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin (3500-50000 K).
        logg : float
            Surface gravity (0.0-5.0).
        metallicity : float
            Metallicity (-2.5 to +0.5 dex).

        Returns
        -------
        bool
            True if all parameters are valid.

        Raises
        ------
        ValueError
            If any parameter is outside the grid range.
        """

        if not self.TEFF_RANGE[0] <= teff <= self.TEFF_RANGE[1]:
            raise ValueError(f"Teff {teff}K is out of range {self.TEFF_RANGE}")
        if not self.LOGG_RANGE[0] <= logg <= self.LOGG_RANGE[1]:
            raise ValueError(f"logg {logg} is out of range {self.LOGG_RANGE}")
        if metallicity < min(self.METALLICITIES) or metallicity > max(self.METALLICITIES):
            raise ValueError(f"Metallicity {metallicity} is out of range [{min(self.METALLICITIES)}, {max(self.METALLICITIES)}]")

        return True
    
        
    def find_model(self, teff, logg, metallicity):
        """
        Find the closest CK04 model by searching the online directory.

        This method uses regex pattern matching on the directory listing to discover
        available models dynamically, rather than relying on hardcoded grid points.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        dict
            Dictionary with keys:
            - 'teff': Closest available Teff
            - 'logg': Requested log g (matched within FITS file)
            - 'metallicity': Closest available metallicity
            - 'dirname': Directory name for this metallicity
            - 'filename': FITS filename

        Raises
        ------
        ValueError
            If no models are found in the directory.
        """

        self.validate_parameters(teff, logg, metallicity)

        closest_metallicity = min(self.METALLICITIES, key=lambda x: abs(x - metallicity))
        dirname = self._fehtodir(closest_metallicity)
        dir_url = f"{self.base_url}/{dirname}"

        response = requests.get(dir_url)
        response.raise_for_status()

        matches = re.findall(rf'{dirname}_(\d+)\.fits', response.text)


        if not matches:
            raise ValueError(f"No models found in {dirname}")

        available_teffs = [int(m) for m in matches]
        closest_teff = min(available_teffs, key=lambda x: abs(x - teff))

        return {
            'teff': closest_teff,
            'logg': logg, # match in model
            'metallicity': closest_metallicity,
            'dirname': dirname,
            'filename': f"{dirname}_{closest_teff}.fits"
        }
    
    def construct_url(self, teff, logg, metallicity):
        """
        Construct full URL to CK04 model FITS file.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        str
            Complete URL to the FITS file.
        """

        params = self.find_model(teff, logg, metallicity)
        return f"{self.base_url}{params['dirname']}/{params['filename']}"
    
    def parse_fits(self, fitsfile, logg):
        """
        Parse CK04 FITS file and extract spectrum for requested log g.

        CK04 FITS files have a unique structure where each file contains spectra
        for multiple log g values stored as separate columns (g00, g05, g10, ..., g50).

        Parameters
        ----------
        fitsfile : file-like object
            FITS file content (BytesIO object).
        logg : float
            Surface gravity to extract. If exact match not found, uses closest available.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns:
            - 'wavelength': Wavelength in Angstroms
            - 'flux': Flux in erg/s/cm^2/A
        """

        model_data = fits.open(fitsfile)

        data = model_data[1].data
        wavelength = data['WAVELENGTH']  # Angstroms

        # column format: g00, g05, g10, g20, g25, g30, g35, g40, g45, g50

        logg_col = f"g{int(round(logg * 10)):02d}"

        if logg_col in data.dtype.names:
            flux = data[logg_col]
        else:
            available_loggs = []
            for colname in data.dtype.names:
                if colname.startswith('g') and colname[1:].isdigit():
                    available_loggs.append((int(colname[1:]) / 10, colname))

            _, closest_col = min(available_loggs, key=lambda x: abs(x[0] - logg))
            flux = data[closest_col]

        # Convert to native byte order to avoid endianness issues
        wavelength = np.asarray(wavelength, dtype=np.float64)
        flux = np.asarray(flux, dtype=np.float32)

        return pd.DataFrame({
            'wavelength': wavelength,
            'flux': flux
        })

class StellarModel:
    """
    Unified interface for loading stellar atmosphere models from different grids.

    Parameters
    ----------
    grid : str, optional
        Model grid name. Default is 'ck04'.
        Available options: 'ck04'

    Raises
    ------
    ValueError
        If the specified grid is not recognized.
    """

    AVALIABLE_MODELS = {
        'ck04': CKModel,
    }

    def __init__(self, grid='ck04'):

        if grid.lower() not in self.AVALIABLE_MODELS:
            raise ValueError(f"Model grid '{grid}' not recognized. Available options are: {list(self.AVALIABLE_MODELS.keys())}")

        self.model = self.AVALIABLE_MODELS[grid.lower()]()

    def load_model(self, teff, logg, metallicity):
        """
        Load a stellar atmosphere model spectrum.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' (Angstroms) and 'flux' (erg/s/cm^2/A).
        """

        return self.model.load_model(teff, logg, metallicity)

    def find_model(self, teff, logg, metallicity):
        """
        Find the closest available model parameters in the grid.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float
            Metallicity [M/H] in dex.

        Returns
        -------
        dict
            Dictionary containing the closest grid point parameters.
        """

        return self.model.find_model(teff, logg, metallicity)
