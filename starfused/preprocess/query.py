"""
Photometric data query module for astronomical sources.
Supports querying from across multiple catalogs via VizieR.
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.exceptions import TimeoutError
import warnings


class PhotometricQuery:
    """
    Query photometric data for astronomical sources across multiple surveys if they are available.

    Supports UV to IR photometry from major surveys including:
    - UV: GALEX, SWIFT, HST
    - Optical: SDSS, Gaia, Pan-STARRS, HST
    - Near-IR: 2MASS, HST, JWST
    - Mid-IR: WISE, JWST

    All effective wavelength data comes from SVO Filter Profile Service (http://svo2.cab.inta-csic.es)
    """

    # Catalog definitions with VizieR catalog IDs and bandpass information (using SVO Filter Profile Service)
    CATALOGS = {
        'SWIFT': {
            'vizier_id': 'II/339/uvotssc1',  # SWIFT UVOT Serendipitous Source Catalogue
            'bands': ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V'],
            'flux_cols': ['UVW2mag', 'UVM2mag', 'UVW1mag', 'Umag', 'Bmag', 'Vmag'],
            'err_cols': ['e_UVW2mag', 'e_UVM2mag', 'e_UVW1mag', 'e_Umag', 'e_Bmag', 'e_Vmag'],
            'wavelengths': [2054.61, 2246.43, 2580.74, 3501.25, 4329.89, 5402.49] * u.AA
        },
        'GALEX': {
            'vizier_id': 'II/335/galex_ais',
            'bands': ['FUV', 'NUV'],
            'flux_cols': ['FUV', 'NUV'],
            'err_cols': ['e_FUV', 'e_NUV'],
            'wavelengths': [1548.85, 2303.37] * u.AA
        },
        'HST_WFC3_UVIS': {
            'vizier_id': 'II/369/hscv2',  # Hubble Source Catalog v3
            'bands': ['F218W', 'F225W', 'F275W', 'F336W', 'F438W', 'F606W', 'F814W'],
            'flux_cols': ['F218Wmag', 'F225Wmag', 'F275Wmag', 'F336Wmag', 'F438Wmag', 'F606Wmag', 'F814Wmag'],
            'err_cols': ['e_F218Wmag', 'e_F225Wmag', 'e_F275Wmag', 'e_F336Wmag', 'e_F438Wmag', 'e_F606Wmag', 'e_F814Wmag'],
            'wavelengths': [2222.48, 2372.81, 2720.03, 3358.95, 4323.35, 5782.20, 7964.25] * u.AA
        },
        'HST_WFC3_IR': {
            'vizier_id': 'II/369/hscv2',
            'bands': ['F098M', 'F105W', 'F110W', 'F125W', 'F140W', 'F160W'],
            'flux_cols': ['F098Mmag', 'F105Wmag', 'F110Wmag', 'F125Wmag', 'F140Wmag', 'F160Wmag'],
            'err_cols': ['e_F098Mmag', 'e_F105Wmag', 'e_F110Wmag', 'e_F125Wmag', 'e_F140Wmag', 'e_F160Wmag'],
            'wavelengths': [9826.81, 10430.83, 11200.52, 12363.55, 13734.66, 15278.47] * u.AA
        },
        'HST_ACS_WFC': {
            'vizier_id': 'II/369/hscv2',
            'bands': ['F435W', 'F475W', 'F606W', 'F814W', 'F850LP'],
            'flux_cols': ['F435Wmag', 'F475Wmag', 'F606Wmag', 'F814Wmag', 'F850LPmag'],
            'err_cols': ['e_F435Wmag', 'e_F475Wmag', 'e_F606Wmag', 'e_F814Wmag', 'e_F850LPmag'],
            'wavelengths': [4341.62, 4708.87, 5809.26, 7973.39, 9004.99] * u.AA
        },
        'JWST_NIRCam': {
            'vizier_id': 'J/A+A/690/A240',  # ASTRODEEP-JWST catalog
            'bands': ['F070W', 'F090W', 'F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W'],
            'flux_cols': ['F070Wmag', 'F090Wmag', 'F115Wmag', 'F150Wmag', 'F200Wmag', 'F277Wmag', 'F356Wmag', 'F444Wmag'],
            'err_cols': ['e_F070Wmag', 'e_F090Wmag', 'e_F115Wmag', 'e_F150Wmag', 'e_F200Wmag', 'e_F277Wmag', 'e_F356Wmag', 'e_F444Wmag'],
            'wavelengths': [6988.43, 8984.98, 11433.62, 14872.56, 19680.41, 27278.58, 35287.04, 43504.26] * u.AA
        },
        'Gaia': {
            'vizier_id': 'I/355/gaiadr3',  # Gaia DR3
            'bands': ['G', 'BP', 'RP'],
            'flux_cols': ['Gmag', 'BPmag', 'RPmag'],
            'err_cols': ['e_Gmag', 'e_BPmag', 'e_RPmag'],
            'wavelengths': [5822.39, 5035.75, 7619.96] * u.AA
        },
        'SDSS': {
            'vizier_id': 'V/154/sdss16',
            'bands': ['u', 'g', 'r', 'i', 'z'],
            'flux_cols': ['umag', 'gmag', 'rmag', 'imag', 'zmag'],
            'err_cols': ['e_umag', 'e_gmag', 'e_rmag', 'e_imag', 'e_zmag'],
            'wavelengths': [3543, 4770, 6231, 7625, 9134] * u.AA
        },
        'PanSTARRS': {
            'vizier_id': 'II/349/ps1',
            'bands': ['g', 'r', 'i', 'z', 'y'],
            'flux_cols': ['gmag', 'rmag', 'imag', 'zmag', 'ymag'],
            'err_cols': ['e_gmag', 'e_rmag', 'e_imag', 'e_zmag', 'e_ymag'],
            'wavelengths': [4866, 6215, 7545, 8679, 9633] * u.AA
        },
        '2MASS': {
            'vizier_id': 'II/246/out',
            'bands': ['J', 'H', 'K'],
            'flux_cols': ['Jmag', 'Hmag', 'Kmag'],
            'err_cols': ['e_Jmag', 'e_Hmag', 'e_Kmag'],
            'wavelengths': [12350, 16620, 21590] * u.AA
        },
        'WISE': {
            'vizier_id': 'II/328/allwise',
            'bands': ['W1', 'W2', 'W3', 'W4'],
            'flux_cols': ['W1mag', 'W2mag', 'W3mag', 'W4mag'],
            'err_cols': ['e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag'],
            'wavelengths': [34000, 46000, 120000, 220000] * u.AA
        }
    }

    def __init__(self, row_limit=10):
        """
        Initialize the photometric query interface.

        Parameters
        ----------
        row_limit : int, optional
            Maximum number of rows to return per catalog query (default: 10)
        """
        self.vizier = Vizier(row_limit=row_limit, columns=['**'])
        self.vizier.ROW_LIMIT = row_limit

    def query_source(self, coord, radius=5*u.arcsec, catalogs='all'):
        """
        Query photometric data for a source.

        Parameters
        ----------
        coord : SkyCoord or tuple
            Source coordinates. Can be SkyCoord object or (ra, dec) tuple in degrees.
        radius : Quantity, optional
            Search radius (default: 5 arcsec)
        catalogs : str or list, optional
            Catalogs to query. Use 'all' for all catalogs, or provide list of
            catalog names (e.g., ['GALEX', '2MASS', 'WISE'])

        Returns
        -------
        dict
            Dictionary with catalog names as keys and photometric data as values.
            Each entry contains 'bands', 'mags', 'errors', 'wavelengths'.
        """
        # Convert coordinate if necessary
        if not isinstance(coord, SkyCoord):
            coord = SkyCoord(ra=coord[0]*u.deg, dec=coord[1]*u.deg, frame='icrs')

        # Determine which catalogs to query
        if catalogs == 'all':
            catalogs_to_query = list(self.CATALOGS.keys())
        else:
            catalogs_to_query = catalogs if isinstance(catalogs, list) else [catalogs]

        results = {}

        for cat_name in catalogs_to_query:
            if cat_name not in self.CATALOGS:
                warnings.warn(f"Unknown catalog: {cat_name}. Skipping.")
                continue

            cat_info = self.CATALOGS[cat_name]

            try:
                # Query VizieR
                result = self.vizier.query_region(
                    coord,
                    radius=radius,
                    catalog=cat_info['vizier_id']
                )

                if len(result) == 0 or len(result[0]) == 0:
                    print(f"No data found in {cat_name}")
                    continue

                # Extract the closest match
                table = result[0]
                closest_idx = 0  # Vizier returns sorted by distance

                # Extract photometry
                phot_data = self._extract_photometry(table[closest_idx], cat_info)

                if phot_data is not None:
                    results[cat_name] = phot_data
                    print(f"Found {len(phot_data['bands'])} bands in {cat_name}")

            except TimeoutError:
                warnings.warn(f"Timeout querying {cat_name}")
            except Exception as e:
                warnings.warn(f"Error querying {cat_name}: {str(e)}")

        return results

    def _extract_photometry(self, row, cat_info):
        """
        Extract photometric measurements from a catalog row.

        Parameters
        ----------
        row : astropy.table.Row
            Single row from query result
        cat_info : dict
            Catalog information dictionary

        Returns
        -------
        dict or None
            Dictionary with photometric data or None if no valid data
        """
        bands = []
        mags = []
        errors = []
        wavelengths = []

        is_flux_catalog = 'flux_unit' in cat_info

        for i, band in enumerate(cat_info['bands']):
            flux_col = cat_info['flux_cols'][i]
            err_col = cat_info['err_cols'][i]

            # Check if data exists and is not masked
            if flux_col in row.colnames and not np.ma.is_masked(row[flux_col]):
                flux_val = row[flux_col]

                # Skip if NaN or zero
                if np.isfinite(flux_val) and flux_val != 0:
                    bands.append(band)

                    if is_flux_catalog:
                        # Convert flux to magnitude (assuming Jy)
                        # Using AB magnitude system: m = -2.5*log10(f_Jy) + 8.90
                        mags.append(-2.5 * np.log10(flux_val) + 8.90)
                    else:
                        mags.append(float(flux_val))

                    # Get error if available
                    if err_col in row.colnames and not np.ma.is_masked(row[err_col]):
                        errors.append(float(row[err_col]))
                    else:
                        errors.append(np.nan)

                    wavelengths.append(cat_info['wavelengths'][i])

        if len(bands) == 0:
            return None

        return {
            'bands': bands,
            'mags': np.array(mags),
            'errors': np.array(errors),
            'wavelengths': np.array(wavelengths)
        }

    def get_sed(self, coord, radius=5*u.arcsec, catalogs='all'):
        """
        Get complete SED (Spectral Energy Distribution) for a source.

        Queries all specified catalogs and combines into a single SED sorted by wavelength.

        Parameters
        ----------
        coord : SkyCoord or tuple
            Source coordinates
        radius : Quantity, optional
            Search radius (default: 5 arcsec)
        catalogs : str or list, optional
            Catalogs to query (default: 'all')

        Returns
        -------
        dict
            Combined SED with keys: 'wavelengths', 'mags', 'errors', 'bands', 'catalogs'
        """
        results = self.query_source(coord, radius=radius, catalogs=catalogs)

        if not results:
            return None

        # Combine all photometry
        all_wavelengths = []
        all_mags = []
        all_errors = []
        all_bands = []
        all_catalogs = []

        for cat_name, phot_data in results.items():
            # Convert wavelengths to common unit (Angstroms) for sorting
            wls = phot_data['wavelengths']
            if hasattr(wls, 'unit'):
                wls_aa = wls.to(u.AA).value
            else:
                wls_aa = wls  # Already in Angstroms

            all_wavelengths.extend(wls_aa)
            all_mags.extend(phot_data['mags'])
            all_errors.extend(phot_data['errors'])
            all_bands.extend(phot_data['bands'])
            all_catalogs.extend([cat_name] * len(phot_data['bands']))

        # Sort by wavelength
        sorted_idx = np.argsort(all_wavelengths)

        return {
            'wavelengths': np.array(all_wavelengths)[sorted_idx] * u.AA,
            'mags': np.array(all_mags)[sorted_idx],
            'errors': np.array(all_errors)[sorted_idx],
            'bands': np.array(all_bands)[sorted_idx],
            'catalogs': np.array(all_catalogs)[sorted_idx]
        }

    def query_by_name(self, name, radius=1*u.arcsec, catalogs='all'):
        """
        Query photometric data by source name (uses Simbad for name resolution).

        Parameters
        ----------
        name : str
            Source name (e.g., 'Vega', 'M42', 'NGC 1234')
        radius : Quantity, optional
            Search radius (default: 5 arcsec)
        catalogs : str or list, optional
            Catalogs to query (default: 'all')

        Returns
        -------
        dict
            Photometric data dictionary
        """
        coord = SkyCoord.from_name(name)
        print(f"Resolved {name} to RA={coord.ra.deg:.6f}, Dec={coord.dec.deg:.6f}")
        return self.query_source(coord, radius=radius, catalogs=catalogs)
