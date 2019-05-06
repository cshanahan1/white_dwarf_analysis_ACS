
from astropy.io import fits
import glob
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import matplotlib.pyplot as plt
from photutils import CircularAnnulus, CircularAperture, RectangularAperture, DAOStarFinder, Background2D
from photutils import DAOStarFinder, aperture_photometry,  deblend_sources, detect_sources, detect_threshold, source_properties
import os
from ginga.util import zscale
from astropy.wcs import WCS

def detect_sources_segm(data, threshold, npixels, kernel_fwhm = 1.8):
    """Use segmentation map to detect sources in image."""

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold=threshold, npixels=npixels, filter_kernel=kernel)
    segm_cat = source_properties(data, segm).to_table()
    return(segm, segm_cat, kernel)

def detect_sources_daofind(data, threshold=10, fwhm=1.8):
    daofind = DAOStarFinder(threshold, fwhm)
    daofind_sources = daofind(data)
    return daofind_sources

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
            f.write(str(source_catalog['xcentroid'][i].value)+', '+str(source_catalog['ycentroid'][i].value)+'\n')

def find_sources(ifile):

    #get stuff from file
    test_file = ifile
    test_dir = test_file.replace(os.path.basename(test_file),'')
    data = fits.open(test_file)[1].data
    hdr = fits.open(test_file)[1].header
    ctx = fits.open(test_file)[3].data
    ndrizim = fits.getval(test_file, 'ndrizim')
    wcss = WCS(fits.open(test_file)[1].header)
    print('Detecting sources in ', test_file)

    #make background image
    print('making background image')
    back_img = Background2D(data, box_size = 100)

    #make threshold image
    print('making threshold image')
    threshold_img = detect_threshold(data, snr=2, background=back_img.background)

    #if only 1 image, treat differently
    if ndrizim == 2:
        print('ndrizim = 2')
        segm, segm_cat, kernel = detect_sources_segm(data, threshold_img, 150, kernel_fwhm = 2.5)
    else:
        segm, segm_cat, kernel = detect_sources_segm(data, threshold_img, 75, kernel_fwhm = 2.5)


    #remove sources that fall near areas of image with 0 coverage, aka 'edges'
    bad_sources = filter_bad_sources(ctx, segm_cat['xcentroid'], segm_cat['ycentroid'], 0.98, ndrizim)
    segm_cat['bad_sources'] = bad_sources
    segm_cat = segm_cat[segm_cat['bad_sources'] == 1]

    #write catalog of sources after filtering out sources near edge only
    write_source_catalog(segm_cat, test_dir + 'sources_segmap.csv', wcss)

    #deblend sources
    segm_deblend = deblend_sources(data, segm, npixels=60, filter_kernel=kernel, nlevels=3,
                                   contrast=0.4)
    segm_deblend_cat = source_properties(data, segm_deblend).to_table()



    #remove sources from deblended catalog that fall near areas of image with 0 coverage, aka 'edges'
    bad_sources = filter_bad_sources(ctx, segm_deblend_cat['xcentroid'], segm_deblend_cat['ycentroid'], 0.98, ndrizim)
    segm_deblend_cat['bad_sources'] = bad_sources
    segm_deblend_cat = segm_deblend_cat[segm_deblend_cat['bad_sources'] == 1]

    #write catalog of deblended sources
    write_source_catalog(segm_deblend_cat, test_dir + 'sources_segmap_deblended.csv', wcss)

    #detect point sources with daofind, filter sources on edge and write catalog
    daofind_sources = detect_sources_daofind(data)
    print(len(daofind_sources))

    bad_sources = filter_bad_sources(ctx, daofind_sources['xcentroid'], daofind_sources['ycentroid'], 0.98, ndrizim, val=False)
    daofind_sources['bad_sources'] = bad_sources
    daofind_sources = daofind_sources[daofind_sources['bad_sources'] == 1]

    with open(test_dir + 'sources_daofind.csv', 'w') as f:
        f.write('x, y\n')
        for i, val in enumerate(daofind_sources):
            f.write(str(daofind_sources['xcentroid'][i])+', '+str(daofind_sources['ycentroid'][i])+'\n')

if __name__ == '__main__':
    test_files = glob.glob('/Users/cshanahan/Desktop/WD_acs/data/cmd_p9/F775W/*drc.fits')
    for i, f in enumerate(test_files):
        print(i, ' of ', len(test_files))
        find_sources(f)
