#%%
"""
    Preprocess photometric data from the Vizier SED service.
"""

import pandas as pd
import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
import requests

class PhotometricQuery:

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
    



# %%
