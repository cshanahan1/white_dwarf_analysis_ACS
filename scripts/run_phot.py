#create photometry catalog from drizzled file.
#sources are detetected with image segmentation in photutils
import glob
import os
from multiprocessing import Pool
import sys
sys.path.append('/Users/cshanahan/Desktop/WD_acs')

from astropy.io import ascii, fits
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from ginga.util import zscale
import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAnnulus, CircularAperture, RectangularAperture
from photutils import aperture_photometry, detect_sources, source_properties
from wfc3_photometry import photometry_tools

def detect_sources_segm(data, threshold):
	"""Use segmentation map to detect sources in image."""

	sigma = 1.8 * gaussian_fwhm_to_sigma
	kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
	kernel.normalize()
	segm = detect_sources(data, threshold=threshold, npixels=75, filter_kernel=kernel)
	segm_cat = source_properties(data, segm).to_table()
	return segm_cat

def filter_bad_sources(ctx, xpos, ypos, n_exp_min, n_exp_tot):
	"""Filters out sources that fall in areas of low coverage. Uses the context
	   image to create an image that shows, pixel wise, the number of images
	   that went into that pixel. For each source, the average number of images
	   per pixel in a 5x5 pixel box around the source is calculated. The source
	   is tossed out if that average is below n_ext_min. N_exp_tot is the number
	   of images. A full frame exposure counts as 2, one for each chip."""

	nimages_array = np.zeros(ctx.shape)
	for i in range(n_exp_tot + 1):
		nimages_array =  nimages_array + (np.bitwise_and(ctx, 2**i) / 2**i)

	bad_detection_flags = []
	for i, xx in enumerate(xpos):
		yy = ypos[i]
		xx, yy = xx.value, yy.value
		ap = RectangularAperture((xx, yy), w=5, h=5.)
		ap_phot = aperture_photometry(nimages_array, ap)
		summ = ap_phot['aperture_sum'][0]
		if summ/ap.area() > n_exp_min:
			bad_detection_flags.append(1)
		else:
			bad_detection_flags.append(0)
	return bad_detection_flags

def output_detected_sources_plot(data, source_catalog, save=True,
								 output_file_path=None, show=False):
	"""Outputs a png of science extension with detected sources overplotted."""

	z1, z2 = zscale.zscale(data)
	plt.scatter(source_catalog['xcentroid'], source_catalog['ycentroid'],
				facecolor='None', edgecolor='r')
	plt.imshow(data, origin='lower', vmin=z1, vmax=z2, cmap='Greys_r')

	if save:
		plt.savefig(output_file_path)
	if show:
		plt.show()

def write_source_catalog(source_catalog, output_file_path):
	"""Writes catalog of output sources as csv for use with ds9."""

	with open(output_file_path, 'w') as f:
		f.write('x, y\n')
		for i, val in enumerate(source_catalog):
			f.write(str(source_catalog['xcentroid'][i])+', '+str(source_catalog['ycentroid'][i])+'\n')

def circular_aperture_photometry(data, source_catalog, r_ap, r_in_sky, r_out_sky):
	"""Aperture photometry with circular aperture of radius r. Aperture sums
	are not background subtracted. """

	apertures = CircularAperture((source_catalog['xcentroid'],
								source_catalog['ycentroid']),
								r=r_ap)
	sky_apertures = CircularAnnulus((source_catalog['xcentroid'],
									source_catalog['ycentroid']),
									r_in_sky,
									r_out_sky)
	medsky = photometry_tools.aperture_stats_tbl(data, sky_apertures)['aperture_median']
	photcat = aperture_photometry(data, apertures)
	photcat['med_sky_annulus'] = medsky
	return(photcat)

def add_fwhm_to_photcat(data, photcat):
	"""Uses wfc3photometry tool to fit a radial profile to every object in
	in the photometry table and adds the FWHM as a column."""

	fwhms = []
	for i, val in enumerate(photcat):
		x, y = photcat['xcenter'][i].value, photcat['ycenter'][i].value
		rad_prof = photometry_tools.rad_prof.RadialProfile(x, y, data)
		fwhms.append(rad_prof.fwhm)
	photcat['fwhm'] = fwhms
	return photcat

def main_run_phot(input_file, show_plots = False, save_plots = True):
	"""Main function to run source finding, photometry, and generation of
	output products."""

	output_dir = os.path.join(os.path.dirname(input_file),'')
	input_file_basename = os.path.basename(input_file).replace('.fits','')
	output_basepath = output_dir+input_file_basename

	#open file
	hdu = fits.open(input_file)
	hd0, hdr, data, ctx = hdu[0].header, hdu[1].header, hdu[1].data, hdu[3].data

	#run source detection
	threshold = np.percentile(data[~np.isnan(data)], 99)
	print(os.path.basename(input_file), threshold)
	source_catalog = detect_sources_segm(data, threshold)

	#filter bad souces
	n_exp_tot = hd0['ndrizim']
	n_exp_min = 0.80

	bad_sources = filter_bad_sources(ctx, source_catalog['xcentroid'],
		source_catalog['ycentroid'], n_exp_min,
		n_exp_tot)
	source_catalog['bad_sources'] = bad_sources

	source_catalog = source_catalog[source_catalog['bad_sources'] == 1]

	write_source_catalog(source_catalog,
			input_file.replace(os.path.basename(input_file),'sources.csv'))
	source_plot_path = output_basepath + '_detected_sources.png'
	output_detected_sources_plot(data, source_catalog, save=save_plots,
								 output_file_path=source_plot_path,
								 show=show_plots)
	#photometry
	photcat = circular_aperture_photometry(data, source_catalog, r_ap=10.,
								r_in_sky=20., r_out_sky=25.)
	#add FWHM to photcat
	photcat = add_fwhm_to_photcat(data, photcat)
	ascii.write(photcat, output_basepath+'_photcat.dat', overwrite = True)

if __name__ == '__main__':
	input_files = glob.glob('/Users/cshanahan/Desktop/clean_desktop/WD_acs/data/*/*/*drc.fits')
	for f in input_files:
		main_run_phot(f)
	# p = Pool(8)
	# p.map(main_run_phot, input_files)
