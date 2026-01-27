import pandas as pd
import numpy as np
import requests
from astropy.io import fits
from io import BytesIO
from abc import ABC, abstractmethod
import re


class BaseModel(ABC):
    """
    Abstract base class for stellar atmosphere model grids.

    This class provides a common interface for loading stellar atmosphere models
    from different grids.

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
    """

    METALLICITIES = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]
    LOGG_RANGE = (0.0, 5.0)
    TEFF_RANGE = (3500, 50000)
    TEFF_STEP = 250
    LOGG_STEP = 0.5



    def __init__(self):
        """Initialize CKModel with STScI repository URL and empty cache."""
        super().__init__("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/")
        self._temperature_cache = {}  # Cache: {dirname: [list of available temperatures]}

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
            Directory name.
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
        Find the closest CK04 model.

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

        # Check cache first
        if dirname not in self._temperature_cache:
            # Fetch directory listing only if not cached
            dir_url = f"{self.base_url}/{dirname}"
            response = requests.get(dir_url)
            response.raise_for_status()

            matches = re.findall(rf'{dirname}_(\d+)\.fits', response.text)

            if not matches:
                raise ValueError(f"No models found in {dirname}")

            # Cache the available temperatures
            self._temperature_cache[dirname] = [int(m) for m in matches]

        # Get temperatures from cache
        available_teffs = self._temperature_cache[dirname]

        # Find closest temperature
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
    
class PhoenixModel(BaseModel):
    """
    Phoenix stellar atmosphere model grid.
    """

    METALLICITIES = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, 0.0, 0.3, 0.5]
    LOGG_RANGE = (0.0, 5.0)


    def __init__(self):
        """Initialize PhoenixModel with STScI repository URL and empty cache."""
        super().__init__("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix")
        self._temperature_cache = {}  

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
            Directory name.
        """

        dirval = int(np.abs(metallicity) * 10)

        if metallicity > 0:
            return f"phoenixp{dirval:02d}"
        else:
            return f"phoenixm{dirval:02d}"

    def validate_parameters(self, teff, logg, metallicity):
        """
        Validate stellar parameters against Phoenix grid ranges.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin (will be validated against available models).
        logg : float
            Surface gravity (0.0-5.0 in 0.5 steps).
        metallicity : float
            Metallicity (must match available values).

        Returns
        -------
        bool
            True if all parameters are valid.

        Raises
        ------
        ValueError
            If any parameter is outside the grid range.
        """

        if not self.LOGG_RANGE[0] <= logg <= self.LOGG_RANGE[1]:
            raise ValueError(f"logg {logg} is out of range {self.LOGG_RANGE}")

        if metallicity < min(self.METALLICITIES) or metallicity > max(self.METALLICITIES):
            raise ValueError(f"Metallicity {metallicity} is out of range [{min(self.METALLICITIES)}, {max(self.METALLICITIES)}]")

        return True

    def find_model(self, teff, logg, metallicity):
        """
        Find the closest Phoenix model.

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

        # Find closest metallicity
        closest_metallicity = min(self.METALLICITIES, key=lambda x: abs(x - metallicity))
        dirname = self._fehtodir(closest_metallicity)

        # Check cache first
        if dirname not in self._temperature_cache:
            # Fetch directory listing only if not cached
            dir_url = f"{self.base_url}/{dirname}"
            response = requests.get(dir_url)
            response.raise_for_status()

            # Extract all temperature values from filenames
            matches = re.findall(rf'{dirname}_(\d+)\.fits', response.text)

            if not matches:
                raise ValueError(f"No models found in {dirname}")

            # Cache the available temperatures
            self._temperature_cache[dirname] = [int(m) for m in matches]

        # Get temperatures from cache
        available_teffs = self._temperature_cache[dirname]

        # Find closest temperature
        closest_teff = min(available_teffs, key=lambda x: abs(x - teff))

        return {
            'teff': closest_teff,
            'logg': logg,  # Will be matched in FITS file
            'metallicity': closest_metallicity,
            'dirname': dirname,
            'filename': f"{dirname}_{closest_teff}.fits"
        }

    def construct_url(self, teff, logg, metallicity):
        """
        Construct full URL to Phoenix model FITS file.

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
        return f"{self.base_url}/{params['dirname']}/{params['filename']}"

    def parse_fits(self, fitsfile, logg):
        """
        Parse Phoenix FITS file and extract spectrum for requested log g.

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

        # Column format: g00, g05, g10, g15, g20, g25, g30, g35, g40, g45, g50

        logg_col = f"g{int(round(logg * 10)):02d}"

        if logg_col in data.dtype.names:
            flux = data[logg_col]
        else:
            # Find closest available log g
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


class KoesterModel(BaseModel):
    """
    Koester DA White Dwarf model grid from the Spanish Virtual Observatory. 
    """

    TEFF_VALUES = list(range(5000, 40000, 250)) + list(range(40000, 80001, 1000))
    LOGG_VALUES = [6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5]

    def __init__(self):
        """Initialize KoesterModel with SVO repository URL and empty cache."""
        super().__init__("http://svo2.cab.inta-csic.es/theory/newov2")
        self._fid_cache = {}  # Cache: {(teff, logg): fid}
        self._session = None

    def _get_session(self):
        """Get or create a requests session for persistent connections."""
        if self._session is None:
            self._session = requests.Session()
        return self._session

    def validate_parameters(self, teff, logg, metallicity=None):
        """
        Validate stellar parameters against Koester grid ranges.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin (5000-80000 K).
        logg : float
            Surface gravity (6.5-9.5).
        metallicity : float, optional
            Ignored for Koester models.

        Returns
        -------
        bool
            True if all parameters are valid.

        Raises
        ------
        ValueError
            If any parameter is outside the grid range.
        """
        if not 5000 <= teff <= 80000:
            raise ValueError(f"Teff {teff}K is out of range (5000, 80000)")
        if not 6.5 <= logg <= 9.5:
            raise ValueError(f"logg {logg} is out of range (6.5, 9.5)")
        return True

    def _query_fid(self, teff, logg):
        """
        Query SVO to get the file ID for a specific Teff/logg combination.

        Parameters
        ----------
        teff : int
            Effective temperature in Kelvin (must be a grid point).
        logg : float
            Surface gravity (must be a grid point).

        Returns
        -------
        str
            File ID for downloading the spectrum.
        """
        cache_key = (teff, logg)
        if cache_key in self._fid_cache:
            return self._fid_cache[cache_key]

        session = self._get_session()
        url = f"{self.base_url}/index.php"

        # Format logg 
        logg_str = str(logg) if logg != int(logg) else str(int(logg))

        data = {
            "models": ",koester2",
            "params[koester2][teff][min]": str(teff),
            "params[koester2][teff][max]": str(teff),
            "params[koester2][logg][min]": logg_str,
            "params[koester2][logg][max]": logg_str,
            "nres": "10",
            "boton": "Search"
        }

        response = session.post(url, data=data)
        response.raise_for_status()

        # Extract fid from response
        matches = re.findall(r'fid=(\d+)', response.text)
        if not matches:
            raise ValueError(f"No model found for Teff={teff}, logg={logg}")

        fid = matches[0]
        self._fid_cache[cache_key] = fid
        return fid

    def find_model(self, teff, logg, metallicity=None):
        """
        Find the closest Koester model to the requested parameters.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float, optional
            Ignored for Koester models (pure H atmospheres).

        Returns
        -------
        dict
            Dictionary with keys:
            - 'teff': Closest available Teff
            - 'logg': Closest available log g
            - 'fid': File ID for downloading
        """
        self.validate_parameters(teff, logg, metallicity)

        # Find closest Teff
        closest_teff = min(self.TEFF_VALUES, key=lambda x: abs(x - teff))

        # Find closest logg
        closest_logg = min(self.LOGG_VALUES, key=lambda x: abs(x - logg))

        # Get the file ID
        fid = self._query_fid(closest_teff, closest_logg)

        return {
            'teff': closest_teff,
            'logg': closest_logg,
            'metallicity': None,  # Pure H atmosphere
            'fid': fid
        }

    def construct_url(self, teff, logg, metallicity=None):
        """
        Construct full URL to Koester model ASCII file.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float, optional
            Ignored for Koester models.

        Returns
        -------
        str
            Complete URL to the ASCII spectrum file.
        """
        params = self.find_model(teff, logg, metallicity)
        return f"{self.base_url}/ssap.php?model=koester2&fid={params['fid']}&format=ascii"

    def load_model(self, teff, logg, metallicity=None):
        """
        Load a Koester DA white dwarf model spectrum.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' (Angstroms) and 'flux' (erg/s/cm^2/A).
        """
        params = self.find_model(teff, logg, metallicity)
        print(f"Loading model: Teff={params['teff']}K, log_g={params['logg']}")

        url = self.construct_url(params['teff'], params['logg'])

        try:
            session = self._get_session()
            response = session.get(url)
            response.raise_for_status()

            return self._parse_ascii(response.text)
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {url}: {e}")

    def parse_fits(self, content, logg):
        """For ABC but unused in this context."""
        return self._parse_ascii(content)

    def _parse_ascii(self, content):
        """
        Parse Koester ASCII spectrum data.

        Parameters
        ----------
        content : str
            ASCII content from SVO.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' and 'flux'.
        """
        lines = content.strip().split('\n')

        wavelengths = []
        fluxes = []

        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))

        return pd.DataFrame({
            'wavelength': wavelengths,
            'flux': fluxes
        })


class BTSettlModel(BaseModel):
    """
    BT-Settl model grid from the Spanish Virtual Observatory.

    BT-Settl models are designed for cool stars and brown dwarfs, including
    cloud formation and settling. They cover very low temperatures suitable
    for L, T, and Y dwarfs.

    Notes
    -----
    - Teff range: 400 - 70000 K (covers brown dwarfs to hot stars)
    - log g range: -0.5 - 6.0
    - Metallicity: -4.0 to +0.5 dex
    - Includes cloud/dust formation physics
    - Note: Not all parameter combinations are available in the grid

    References
    ----------
    Allard, F., Homeier, D., & Freytag, B. 2012, RSPTA, 370, 2765
    """

    # Grid parameters (from SVO interface)
    TEFF_RANGE = (400, 70000)
    LOGG_RANGE = (-0.5, 6.0)
    METALLICITIES = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5]

    def __init__(self):
        """Initialize BTSettlModel with SVO repository URL and empty cache."""
        super().__init__("http://svo2.cab.inta-csic.es/theory/newov2")
        self._fid_cache = {}  # Cache: {(teff, logg, metallicity): fid}
        self._session = None

    def _get_session(self):
        """Get or create a requests session for persistent connections."""
        if self._session is None:
            self._session = requests.Session()
        return self._session

    def validate_parameters(self, teff, logg, metallicity=0.0):
        """
        Validate stellar parameters against BT-Settl grid ranges.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin (400-7000 K).
        logg : float
            Surface gravity (2.5-5.5).
        metallicity : float, optional
            Metallicity [M/H] in dex (-2.5 to +0.5).

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

    def _query_fid(self, teff, logg, metallicity):
        """
        Query SVO to get the file ID for a specific parameter combination.

        Parameters
        ----------
        teff : int
            Effective temperature in Kelvin.
        logg : float
            Surface gravity.
        metallicity : float
            Metallicity [M/H].

        Returns
        -------
        str
            File ID for downloading the spectrum.
        """
        cache_key = (teff, logg, metallicity)
        if cache_key in self._fid_cache:
            return self._fid_cache[cache_key]

        session = self._get_session()
        url = f"{self.base_url}/index.php"

        # Format parameters
        logg_str = str(logg) if logg != int(logg) else str(int(logg))
        met_str = str(metallicity) if metallicity != int(metallicity) else str(int(metallicity))

        data = {
            "models": ",bt-settl",
            "params[bt-settl][teff][min]": str(teff),
            "params[bt-settl][teff][max]": str(teff),
            "params[bt-settl][logg][min]": logg_str,
            "params[bt-settl][logg][max]": logg_str,
            "params[bt-settl][meta][min]": met_str,
            "params[bt-settl][meta][max]": met_str,
            "nres": "10",
            "boton": "Search"
        }

        response = session.post(url, data=data)
        response.raise_for_status()

        # Extract fid from response
        matches = re.findall(r'fid=(\d+)', response.text)
        if not matches:
            raise ValueError(f"No model found for Teff={teff}, logg={logg}, [M/H]={metallicity}")

        fid = matches[0]
        self._fid_cache[cache_key] = fid
        return fid

    def find_model(self, teff, logg, metallicity=0.0):
        """
        Find the closest BT-Settl model to the requested parameters.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float, optional
            Metallicity [M/H] in dex (default: 0.0).

        Returns
        -------
        dict
            Dictionary with keys:
            - 'teff': Closest available Teff
            - 'logg': Closest available log g
            - 'metallicity': Closest available metallicity
            - 'fid': File ID for downloading
        """
        self.validate_parameters(teff, logg, metallicity)

        # BT-Settl Teff grid: 400-2000 (100K steps), 2000-7000 (100K steps)
        if teff < 2000:
            teff_grid = list(range(400, 2100, 100))
        else:
            teff_grid = list(range(2000, 7100, 100))

        # Find closest Teff
        closest_teff = min(teff_grid, key=lambda x: abs(x - teff))

        # logg grid: -0.5 to 6.0 in 0.5 steps (per SVO website)
        logg_grid = [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
        closest_logg = min(logg_grid, key=lambda x: abs(x - logg))

        # Find closest metallicity
        closest_met = min(self.METALLICITIES, key=lambda x: abs(x - metallicity))

        # Get the file ID
        fid = self._query_fid(closest_teff, closest_logg, closest_met)

        return {
            'teff': closest_teff,
            'logg': closest_logg,
            'metallicity': closest_met,
            'fid': fid
        }

    def construct_url(self, teff, logg, metallicity=0.0):
        """
        Construct full URL to BT-Settl model ASCII file.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin.
        logg : float
            Surface gravity (log g).
        metallicity : float, optional
            Metallicity [M/H] in dex.

        Returns
        -------
        str
            Complete URL to the ASCII spectrum file.
        """
        params = self.find_model(teff, logg, metallicity)
        return f"{self.base_url}/ssap.php?model=bt-settl&fid={params['fid']}&format=ascii"

    def load_model(self, teff, logg, metallicity=0.0):
        """
        Load a BT-Settl brown dwarf/cool star model spectrum.

        Parameters
        ----------
        teff : float
            Effective temperature in Kelvin (400-7000 K).
        logg : float
            Surface gravity (log g) in cgs units.
        metallicity : float, optional
            Metallicity [M/H] in dex (default: 0.0).

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' (Angstroms) and 'flux' (erg/s/cm^2/A).
        """
        params = self.find_model(teff, logg, metallicity)
        print(f"Loading model: Teff={params['teff']}K, log_g={params['logg']}, [M/H]={params['metallicity']}")

        url = self.construct_url(params['teff'], params['logg'], params['metallicity'])

        try:
            session = self._get_session()
            response = session.get(url)
            response.raise_for_status()

            return self._parse_ascii(response.text)
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {url}: {e}")

    def parse_fits(self, content, logg):
        """Required by ABC but not used - BT-Settl uses ASCII format from SVO."""
        return self._parse_ascii(content)

    def _parse_ascii(self, content):
        """
        Parse BT-Settl ASCII spectrum data.

        Parameters
        ----------
        content : str
            ASCII content from SVO (whitespace-delimited wavelength, flux).

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns 'wavelength' and 'flux'.
        """
        lines = content.strip().split('\n')

        wavelengths = []
        fluxes = []

        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                wavelengths.append(float(parts[0]))
                fluxes.append(float(parts[1]))

        return pd.DataFrame({
            'wavelength': wavelengths,
            'flux': fluxes
        })


class StellarModel:
    """
    Unified interface for loading stellar atmosphere models from different grids.

    Parameters
    ----------
    grid : str, optional

    Raises
    ------
    ValueError
        If the specified grid is not recognized.
    """

    AVALIABLE_MODELS = {
        'ck04': CKModel,
        'phoenix': PhoenixModel,
        'koester': KoesterModel,
        'bt-settl': BTSettlModel,
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



