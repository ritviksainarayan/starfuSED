"""
    Preprocess photometric data from the Vizier SED service.
"""

import pandas as pd
import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa.irsa_dust import IrsaDust
import requests
from astropy.io import fits
from io import BytesIO
from scipy.interpolate import interp1d

EXTINCTION_CURVES = {
    'mw'        : 'milkyway_diffuse_001.fits',
    'milky_way' : 'milkyway_diffuse_001.fits',
    'mw_dense'  : 'milkyway_dense_001.fits',
    'lmc'       : 'lmc_diffuse_001.fits',
    'lmc_30dor' : 'lmc_30dorshell_001.fits',
    'smc_bar'   : 'smc_bar_001.fits',
}

class Photometry:



    def query(name=None, ra=None, dec=None, filters=['GALEX','UVOT','SDSS','Gaia','2MASS','WISE'], radius=1):
        """
        Query VizieR SED service by object name or coordinates. 

        Parameters
        ----------
        name : str, optional
            Source name.

        ra : float, optional
            Right ascension in degrees.

        dec : float, optional
            Declination in degrees.

        filters : list of str
            List of instrument names to filter by. Filters are applied to the 'sed_filter' column. 
            Example: ['2MASS', 'WISE', 'GALEX']

        radius : float, optional
            Radius in arcsec.

        Returns
        -------
        pandas.DataFrame 
            Filtered DataFrame containing only rows with specified instruments.
                - for duplicate filters, only the row with the smallest flux uncertainty is kept.
                - flux is converted from Jy to erg/s/cm^2/A.
        """

        if name is not None:
            df = Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={name.replace(' ', '')}&-c.rs={radius}").to_pandas().replace(np.nan, 0.0).dropna(subset='sed_filter').drop(columns=['_ID','_time','_etime'])
        elif ra is not None and dec is not None:
            df = Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={ra}+{dec}&-c.rs={radius}").to_pandas().replace(np.nan, 0.0).dropna(subset='sed_filter').drop(columns=['_ID','_time','_etime'])
        else:
            raise ValueError("Either name or (ra, dec) must be provided.")
        
        # load only specified filters
        df = df.loc[df['sed_filter'].str.contains('|'.join(filters), case=True, na=False)].reset_index(drop=True)

        # if uncertainty is 0, then take std of multiple measurements if available
        for filtername in df.loc[df.sed_eflux == 0,'sed_filter'].unique():
            if len(df.loc[df.sed_filter == filtername]) > 1:
                std = df.loc[df.sed_filter == filtername,'sed_flux'].std()
                df.loc[df.sed_filter == filtername,'sed_eflux'] = std

        # duplicate filter removal
        df = df.loc[(df.groupby('sed_filter')['sed_eflux'].idxmin())].query('sed_eflux > 0').sort_values('sed_freq',ascending=False).reset_index(drop=True)  
    
        # flux conversion from Jy to erg/s/cm^2/A
        df['sed_flux'] = df['sed_flux'].values * u.Jy.to(u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(df['sed_freq'].values * u.GHz))
        df['sed_eflux'] = df['sed_eflux'].values * u.Jy.to(u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(df['sed_freq'].values * u.GHz))
        df['sed_wave'] = np.round((c / (df['sed_freq'].values * u.GHz)).to(u.AA).value, 0)

        
        return df
    
    def dust_correction(df, extinction='mw', dustmap='SandF'):
        """
        Apply dust extinction correction to photometric data. 
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing photometric data with columns 'sed_wave', 'sed_flux', 'sed_eflux', '_RAJ2000', '_DEJ2000'.
        extinction : str, optional
            Extinction curve to use. Options are 'mw' or 'milky_way', 'mw_dense', 'lmc', 'lmc_30dor', 'smc_bar'.
        dustmap : str, optional
            Dust map to use for E(B-V) retrieval. Options are 'SFD' or 'SandF'.

        Returns
        -------
        pandas.DataFrame 
            DataFrame with dust-corrected 'sed_flux' and 'sed_eflux' columns.
        """

        ra = df['_RAJ2000'].values[0]
        dec = df['_DEJ2000'].values[0]

        # query coordinate
        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs', equinox='J2000.0')

        dusttable = IrsaDust.get_query_table(coord, radius=2*u.deg, section='ebv')
        
        if dustmap.upper() == 'SFD':
            ebv = dusttable['ext SFD mean'][0]
        elif dustmap.upper() == 'SANDF':
            ebv = dusttable['ext SandF mean'][0]
        else:
            raise ValueError("dustmap must be either 'SFD' or 'SandF'")

        print(f"Retrived E(B-V) = {ebv:.3f} from {dustmap} dust map.")
        
        if extinction.lower() not in EXTINCTION_CURVES.keys():
            raise ValueError(f"Extinction curve '{extinction}' not recognized. Available options are: {list(EXTINCTION_CURVES.keys())}")
        
        url = f"https://archive.stsci.edu/hlsps/reference-atlases/cdbs/extinction/{EXTINCTION_CURVES[extinction.lower()]}"
        
        response = requests.get(url)

        fitsdata = fits.open(BytesIO(response.content))

        fitswave = (1 / fitsdata[1].data['WAVELENGTH']) * 1e4 # in angstroms
        fits_av_ebv = fitsdata[1].data['Av/E(B-V)']

        av_ebv_obs = interp1d(fitswave, fits_av_ebv, kind='linear', bounds_error=False, fill_value="extrapolate")(df['sed_wave'].values)
        a_lambda = ebv * av_ebv_obs

        # correction factor = 10**(0.4 * A_lambda)

        df_corrected = df.copy()
        df_corrected['sed_flux']  *= 10**(0.4 * a_lambda)
        df_corrected['sed_eflux'] *= 10**(0.4 * a_lambda)

        return df_corrected


