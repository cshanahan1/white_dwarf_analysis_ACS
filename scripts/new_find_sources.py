
from astropy.io import fits
import glob
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import matplotlib.pyplot as plt
from photutils import CircularAnnulus, CircularAperture, RectangularAperture, DAOStarFinder, Background2D
from photutils import aperture_photometry,  deblend_sources, detect_sources, detect_threshold, source_properties
import os
from ginga.util import zscale
from astropy.wcs import WCS
from wfc3_photometry import photometry_tools
import cosmics

def detect_sources_segm(data, threshold, npixels, kernel_fwhm = 1.8):
    """Use segmentation map to detect sources in image."""

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold=threshold, npixels=npixels, filter_kernel=kernel)
    segm_cat = source_properties(data, segm).to_table()
    return(segm, segm_cat, kernel)

def filter_bad_sources(ctx_img, xpos, ypos, frac_cover_min, n_exp_total, val = True):

    #boolean array indicating if that pixel is covered by ANY image.
    coverage_arr = np.nan_to_num(ctx_img).astype(bool).astype(int)
    low_coverage_flags = []

    for i, xx in enumerate(xpos):
        if val == False:
            yy = ypos[i]
        else:
            xx, yy = xx.value, ypos[i].value
        ap = RectangularAperture((xx, yy), w=100, h=100)
        ap_phot = aperture_photometry(coverage_arr, ap)
        summ = ap_phot['aperture_sum'][0]
        if summ/ap.area() > frac_cover_min:
            low_coverage_flags.append(1)
        else:
            low_coverage_flags.append(0)

    return low_coverage_flags

def write_source_catalog(source_catalog, output_file_path, wcss):
    """Writes catalog of output sources as csv for use with ds9."""
    with open(output_file_path, 'w') as f:
        f.write('x, y\n')
        for i, val in enumerate(source_catalog):
            ra, dec = wcss.wcs_pix2world(source_catalog['xcentroid'][i].value, source_catalog['ycentroid'][i].value, 1)
            fwhmm = source_catalog['fwhm'][i]
            f.write(str(ra)+', '+str(dec)+', '+str(fwhmm)+'\n')

def add_fwhm_to_photcat(data, photcat):
	"""Uses wfc3photometry tool to fit a radial profile to every object in
	in the photometry table and adds the FWHM as a column."""

	fwhms = []
	for i, val in enumerate(photcat):
		x, y = photcat['xcentroid'][i].value, photcat['ycentroid'][i].value
		rad_prof = photometry_tools.rad_prof.RadialProfile(x, y, data)
		fwhms.append(rad_prof.fwhm)
	photcat['fwhm'] = fwhms
	return photcat

def find_sources(ifile):

    test_file = ifile
    test_dir = test_file.replace(os.path.basename(test_file),'')
    data = fits.open(test_file)[1].data
    hdr = fits.open(test_file)[1].header
    ctx = fits.open(test_file)[3].data

    ndrizim = fits.getval(test_file, 'ndrizim')

    print(test_file)
    back_img = Background2D(data, box_size = 50)

    threshold_img = detect_threshold(data, snr=2, background=back_img.background)
    segm, segm_cat, kernel = detect_sources_segm(data, threshold_img, 100, kernel_fwhm = 2.5)

    bad_sources = filter_bad_sources(ctx, segm_cat['xcentroid'], segm_cat['ycentroid'], 0.98, ndrizim)
    segm_cat['bad_sources'] = bad_sources
    segm_cat = segm_cat[segm_cat['bad_sources'] == 1]
    wcss = WCS(fits.open(test_file)[1].header)

    #fit fwhm
    segm_cat = add_fwhm_to_photcat(data, segm_cat)

    write_source_catalog(segm_cat, test_dir + 'sources_segm.csv', wcss)

    #deblend sources
    segm_deblend = deblend_sources(data, segm, npixels=75, filter_kernel=kernel, nlevels=3,
                                   contrast=0.5)
    segm_deblend_cat = source_properties(data, segm_deblend).to_table()

    #fit fwhm
    segm_deblend_cat = add_fwhm_to_photcat(data, segm_deblend_cat)

    bad_sources = filter_bad_sources(ctx, segm_deblend_cat['xcentroid'], segm_deblend_cat['ycentroid'], 0.98, ndrizim)
    segm_deblend_cat['bad_sources'] = bad_sources
    segm_deblend_cat = segm_deblend_cat[segm_deblend_cat['bad_sources'] == 1]

    segm_deblend_cat = segm_deblend_cat[segm_deblend_cat['fwhm'] > 1.8]
    write_source_catalog(segm_deblend_cat, test_dir + 'sources_segm_deblended.csv', wcss)

    write_source_catalog(segm_deblend_cat, test_dir + 'sources_segm_deblended_fwhm.csv', wcss)

    print('donezo')

if __name__ == '__main__':
    ifiles = glob.glob('/Users/cshanahan/Desktop/clean_desktop/WD_acs/data/dgq*/*/*.fits')
    for i, f in enumerate(ifiles):
        print(i, 'of {}'.format(len(ifiles)))
        find_sources(f)
