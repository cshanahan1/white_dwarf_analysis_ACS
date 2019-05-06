from astropy.io import ascii, fits
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, RectangularAperture
from astropy.wcs import WCS
from rad_prof import RadialProfile
from background_median import aperture_stats_tbl
import numpy as np
import os
import glob
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u


"""
Produces aperutre photometry catalogs, written for processing ACS parallels of
white dwarf fields. Reads in source catalogs that contain x, y position of
detected (point and extended) sources in these images.

For each input drc image, a cosmic ray mask (.mask.fits) should be in the same
directory. This is used to reject CRs detected as sources in regions where
there is only one image of coverage.

By setting 'gaia_match' = True, this works with images that have been aligned
to the GAIA frame with tweakreg. 'BLAH.py' does the alignment procedure as well
as downloads the GAIA catalog for each field into a file called 'gaia.cat'. The
final photometry catalog will have x, y source positions detected in the image,
sky positions in the GAIA frame as well as in the original frame before
alginment. Sources that have also been detected by GAIA will have those
positions as well. Sources with no GAIA match will have the value -999.
"""

def _add_fwhm_to_phot_table(data, phot_table):

	"""Fits a radial profile to every object in the table and returns back
    the table with the column 'fwhm' added. Uses functions in rad_prof.py.

    Parameters
    ----------
    data : 2d array
        Image array.
    phot_table : astropy table
        Photometry table with source position columns 'xcenter' and 'ycenter'.
    Returns
    -------
    phot_table : astropy table
        Photometry table with the column 'fwhm' added. """

	fwhms = []
	for i, val in enumerate(phot_table):
		x, y = phot_table['xcenter'][i], phot_table['ycenter'][i]
		rad_prof = RadialProfile(x, y, data)
		fwhms.append(rad_prof.fwhm)
	phot_table['fwhm'] = fwhms
	return phot_table

def _calculate_sky(data, phot_table, r_in, r_out):

    """Adds sigma-clipped median background and standard deviation in specified
    aperture to photometry table. Uses functions in background_median.py.

    Parameters
    ----------
    data : 2d array
        Image array.
    phot_table : astropy table
        Photometry table with source position columns 'xcenter' and 'ycenter'.
    Returns
    -------
    phot_table : astropy table
        Photometry table with the columns 'back' and 'back_std' added. """

    source_positions = [(phot_table['xcenter'][i], phot_table['ycenter'][i]) \
                        for i, x in enumerate(phot_table)]
    back_apertures = CircularAnnulus(source_positions, r_in, r_out)
    back_tbl = aperture_stats_tbl(data, back_apertures)
    median_sky = back_tbl['aperture_median']
    back_std = back_tbl['aperture_std']

    phot_table['back'] = median_sky
    phot_table['back_std'] = back_std

    return phot_table

def _compute_phot_err_daophot(flux, back, back_rms, phot_ap_area,
                             sky_ap_area, epadu=1.0):

    """Calculate measurement errors in the same manner as IRAF/DAOphot. The
    error terms in this model represent Poisson noise from the source, Poisson
    noise in the sky and readout noise, and error in the sky measurement
    err = sqrt ((flux - back*phot_ap_area)/epadu + phot_ap_area * backrms**2 +
    phot_ap_area**2 * backrms**2 / sky_ap_area)
    Where flux is the aperture sum, back is the per-pixel background level,
    epadu is the conversion factor between e- and adu (gain), phot_ap_area is
    the area of the photometric aperture, backrms is the uncertainty in the sky,
    and sky_ap_area is the sky aperture area. Note that the flux/background
    in the above equation are in ADU, but the input to this function is in e-.
    The conversion is done internally.
    Parameters
    ----------
    flux : float
        Non-sky subtracted flux of star, in electrons.
    back : float
        Per-pixel background level.
    back_rms : float
        Background RMS.
    phot_ap_area : int or float
        Area of photometric aperture.
    sky_ap_area : int or float
        Area of sky aperture.
    epdau : float
        CCD gain, default 1.0.
     Returns
     --------
     errs : tuple
        (error in insturmental magnitudes, error in flux)
    """
    #convert values input in e- to ADU, as equation expects
    #sky subtract flux to isolate poisson noise from source
    flux = (flux - back*phot_ap_area) / epadu
    back = back / epadu
    back_rms = back_rms / epadu

    err1 = (flux/epadu)
    err2 = (phot_ap_area*back_rms**2)
    err3 = (phot_ap_area**2 * back_rms**2) / sky_ap_area

    flux_err_adu = np.sqrt(np.abs(err1 + err2 + err3)) # in ADU
    flux_err = flux_err_adu * epadu

    return flux_err

def _run_aperture_photometry(data, source_catalog,
                             aperture_radii=[2., 3., 5., 7.,]):

    source_positions = [(source_catalog['x'][i], source_catalog['y'][i]) for i,\
                        x in enumerate(source_catalog)]
    aps = [CircularAperture(source_positions, r=rad) for rad in aperture_radii]
    phot_table = Table(aperture_photometry(data, aps))
    phot_table.remove_column('id')

    table_order = []
    for i, rad in enumerate(aperture_radii):
        rad = str(int(rad))
        phot_table.rename_column('aperture_sum_{}'.format(str(i)), \
                                 'countrate_{}'.format(rad))
        table_order.append('countrate_{}'.format(rad))

    table_order = ['xcenter','ycenter'] + table_order

    phot_table = phot_table[table_order]

    return(phot_table)

def _add_cr_flag_to_tab(cr_map, phot_table, ctx, ndrizim):
    nimages_array = np.zeros(ctx.shape)
    for i in range(ndrizim + 1):
        nimages_array =  nimages_array + (np.bitwise_and(ctx, 2**i) / 2**i)
    cr_flags = []
    for i, row in enumerate(phot_table):
        x, y = row['xcenter'], row['ycenter']
        ap = RectangularAperture((x, y), w=2, h=2)
        ap_phot = aperture_photometry(nimages_array, ap)
        nim_avg = ap_phot['aperture_sum'][0]/ap.area()
        if nim_avg > 1.0:
            cr_flags.append(0)
        else:
            #now use same aperture on CR map
            ap_phot_cr = aperture_photometry(cr_map, ap)
            sum_cr = ap_phot_cr['aperture_sum'][0]
            avg_cr = sum_cr/ap.area()
            if avg_cr > 0:
                cr_flags.append(1)
            else:
                cr_flags.append(0)

    phot_table['cr_flag'] = cr_flags
    return phot_table

def make_photometry_catalogs(input_file, source_catalog_path, gaia_match=True,
                             aperture_radii=[2., 4., 6., 8., 10.]):

    #data, headers
    hdu = fits.open(input_file)
    data = hdu[1].data
    ctx_img = hdu[3].data
    hdr0, hdr1 = hdu[0].header, hdu[1].header

    cr_map_file = input_file.replace('.fits', '.convmask.fits')
    cr_map = fits.getdata(cr_map_file)

    #read in catalog with x, y source positions
    source_catalog = ascii.read(source_catalog_path, format='csv')

    #aperture photometry
    phot_table = _run_aperture_photometry(data, source_catalog, aperture_radii = aperture_radii)

    #calculate sky and add to photometry table
    phot_table = _calculate_sky(data, phot_table, r_in = 10., r_out = 15.)

    #add errors to table

    for ap_radius in aperture_radii:
        flux = phot_table['countrate_{}'.format(str(int(ap_radius)))]
        back = phot_table['back']
        back_rms = phot_table['back_std']
        phot_ap_area = np.pi * ap_radius**2.
        sky_ap_area = (np.pi * 15**2.) - (np.pi * 10**2.)
        err = _compute_phot_err_daophot(flux, back, back_rms, phot_ap_area,
                                        sky_ap_area)
        phot_table['err_{}'.format(str(int(ap_radius)))] = err

    #add FWHM to photometry table
    phot_table = _add_fwhm_to_phot_table(data, phot_table)

    #reorder table columns
    new_order = phot_table.colnames
    new_order.remove('fwhm')
    new_order.insert(2, 'fwhm')
    phot_table = phot_table[new_order]

    #add cosmic ray flag to sources in nim=1 areas
    _add_cr_flag_to_tab(cr_map, phot_table, ctx_img, hdr0['ndrizim'])

    if gaia_match:
        #gaia catalog from astroquery
        gaia_cat = ascii.read(os.path.dirname(input_file) + '/gaia.cat')
        ra_gaia_cat, dec_gaia_cat = gaia_cat['ra'], gaia_cat['dec']

        #tweaked wcs positions after alignment to gaia
        gaia_wcs = WCS(hdr1)
        ra_gaia_wcs, dec_gaia_wcs = gaia_wcs.wcs_pix2world(phot_table['xcenter'], phot_table['ycenter'], 1)

        phot_table['RA_GAIA_WCS'] = ra_gaia_wcs
        phot_table['DEC_GAIA_WCS'] = dec_gaia_wcs

        sources_img = SkyCoord(ra=ra_gaia_wcs*u.degree, dec=dec_gaia_wcs*u.degree)
        sources_gaia_cat = SkyCoord(ra=ra_gaia_cat*u.degree, dec=dec_gaia_cat*u.degree)

        idxx, d2ds, d3ds = sources_img.match_to_catalog_sky(sources_gaia_cat)
        d2ds = d2ds.deg
        close_matches_gaia = idxx[d2ds <= 1e-5]
        close_matches_acs = np.where(d2ds <= 1e-5)[0]

        #add gaia catalog ra/dec columns to table
        phot_table['GAIA_CAT_RA'] = np.full(len(phot_table), -999.)
        phot_table['GAIA_CAT_DEC'] = np.full(len(phot_table), -999.)
        phot_table['GAIA_CAT_RA'][close_matches_acs] =sources_gaia_cat.ra.deg[close_matches_gaia]
        phot_table['GAIA_CAT_DEC'][close_matches_acs] = sources_gaia_cat.dec.deg[close_matches_gaia]

    # phot_table_output_path = os.path.dirname(input_file) + '/photometry_catalog_.dat'
    # ascii.write(phot_table, phot_table_output_path, format='csv', overwrite=True)
    # print('Wrote {}.'.format(phot_table_output_path))

    #write table of sources with no CRs
    phot_table_no_crs = phot_table[phot_table['cr_flag'] == 0]
    phot_table_output_path = os.path.dirname(input_file) + '/photometry_catalog_noCRs_FINAL.dat'
    ascii.write(phot_table_no_crs, phot_table_output_path, format='csv', overwrite=True)
    print('Wrote {}.'.format(phot_table_output_path))


if __name__ == '__main__':
    
    # input_file = '/Users/cshanahan/Desktop/WD_acs/data/cmd_p9/F775W/jcmdp9020_drc.fits'
    # source_catalog_path = '/Users/cshanahan/Desktop/WD_acs/data/cmd_p9/F775W/sources_segmap_deblended.csv'
    # make_photometry_catalogs(input_file, source_catalog_path, gaia_match=True)

    input_files = glob.glob('/Users/cshanahan/Desktop/WD_acs/data/*/*/*drc.fits')
    for input_file in input_files:
        print(input_file)
        source_catalog_path = input_file.replace(os.path.basename(input_file), 'inspected_source_cat.csv')
        make_photometry_catalogs(input_file, source_catalog_path)
