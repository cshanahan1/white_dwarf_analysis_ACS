#create photometry catalog from drizzled file.
#sources are detetected with image segmentation in photutils
import glob
import os
import sys
sys.path.append('/Users/cshanahan/Desktop/WD_acs')

from astropy.io import ascii, fits
from astropy.stats import sigma_clipped_stats
from ginga.util import zscale
import numpy as np
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, Background2D, DAOStarFinder, detect_threshold
from wfc3_photometry import photometry_tools

def detect_sources_daostarfinder(data, fwhm = 2.0, threshold = 1.0):
    print('*** Detecting sources in image ***')
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold, exclude_border=True)
    source_catalog = daofind(data)
    return(source_catalog)

def detect_sources_segmap(data):
    pass

def mask_sources_border_chipgap_drizzle_ACS(source_catalog, npix = 50):
    """ Excludes sources in a catalog within npix of the border of a
    drizzled ACS image, accounting for the fact that the image is not a perfect
    rectangle anymore after projection. """

    print('*** Filtering source catalog for sources within {} pixels of top, bottom, and chipgap ***'.format(npix))

    #define equations of lines that trace top and bottom edge
    top_line = lambda x : 0.05 * x + 4115
    bottom_line = lambda x : 0.03 * x

    #exclude sources within npix of top
    filtered_cat = source_catalog[np.abs((top_line(source_catalog['xcentroid']))
                    -source_catalog['ycentroid'])>npix]
    #exclude sources within npix of bottom
    filtered_cat = filtered_cat[np.abs((bottom_line(filtered_cat['xcentroid']))
                    -filtered_cat['ycentroid'])>npix]

    #define line that traces chip gap
    chip_line = lambda x: 0.04 * x + 2078
    filtered_cat = filtered_cat[np.abs((chip_line(filtered_cat['xcentroid']))
                    -filtered_cat['ycentroid'])>npix]

    return(filtered_cat)

def output_detected_sources_plot(data, source_catalog,  save=True, output_file_path=None, show=False):

    print('*** Making plot of detected sources ***')
    z1, z2 = zscale.zscale(data)
    plt.scatter(source_catalog['xcentroid'], source_catalog['ycentroid'],
                facecolor='None', edgecolor='r')
    plt.imshow(data, origin='lower', vmin=z1, vmax=z2, cmap='Greys_r')
    plt.title('detected sources')

    if save:
        plt.savefig(output_file_path)
    if show:
        plt.show()

def write_source_catalog(source_catalog, output_file_path):
    """Writes catalog of output sources as csv for use with ds9."""
    print('*** Writing {} ***'.format(output_file_path))
    with open(output_file_path, 'w') as f:
        f.write('x, y\n')
        for i, val in enumerate(source_catalog):
            f.write(str(source_catalog['xcentroid'][i])+', '+str(source_catalog['ycentroid'][i])+'\n')

def circular_aperture_photometry(data, source_catalog, r_ap, r_in_sky, r_out_sky):

    """Aperture photometry with circular aperture of radius r. Aperture sums
    are not background subtracted, """

    print('*** Running aperture photometry ***')
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

    print('*** Measuring FWHM of sources ***')
    fwhms = []
    for i, val in enumerate(photcat):
        x, y = photcat['xcenter'][i].value, photcat['ycenter'][i].value
        rad_prof = photometry_tools.rad_prof.RadialProfile(x, y, data)
        fwhms.append(rad_prof.fwhm)

    photcat['fwhm'] = fwhms
    return photcat

def main_run_phot(input_file, output_dir, output_file_prefix = '',
show_plots = False, save_plots = True, save_photcat = True, save_source_list = True):
    """Main function to run source finding, photometry, and generation of
    output products."""

    input_file_basename = os.path.basename(input_file).replace('.fits','')
    output_basepath = output_dir+output_file_prefix+'_'+input_file_basename

    #open file
    hdu = fits.open(input_file)
    hdr, data = hdu[1].header, hdu[1].data

    #run source detection
    print('calculating threshold for source detection')
    #ean, median, std = sigma_clipped_stats(data, sigma=2.5)
    threshold = 0.1*(max(data.flatten()))
    print('threshold = ', threshold)
    source_catalog = detect_sources_daostarfinder(data, fwhm=2.0, threshold=threshold)
    source_catalog = mask_sources_border_chipgap_drizzle_ACS(source_catalog, npix=60)
    if save_source_list:
        write_source_catalog(source_catalog, input_file.replace(os.path.basename(input_file),'sources.csv'))
    source_plot_path = output_basepath + '_detected_sources.png'
    output_detected_sources_plot(data, source_catalog, save=save_plots,
                                 output_file_path=source_plot_path,
                                 show=show_plots)
    #photometry
    photcat = circular_aperture_photometry(data, source_catalog, r_ap=10.,
                                r_in_sky=20., r_out_sky=25.)
    #add FWHM to photcat
    photcat = add_fwhm_to_photcat(data, photcat)

    #write photcat, if write = true
    if save_photcat:
        print('*** Saving photometry catalog ***')
        ascii.write(photcat, output_basepath+'_photcat.dat')

if __name__ == '__main__':

    input_dirs = glob.glob('/Users/cshanahan/Desktop/WD_acs/data/cmd*/*')

    for i, dirr in enumerate(input_dirs):
        output_dir = '/Users/cshanahan/Desktop/WD_acs/output/mydrizzle_run1/'
        prog_id_visit = dirr.split('/')[-2]
        filt = dirr.split('/')[-1]
        output_file_prefix = prog_id_visit+'_'+filt
        print(output_file_prefix)
        input_file = glob.glob(dirr + '/*drc.fits')[0]
        print('***********************************************************')
        print('******* PROCESSING {} *******'.format(input_file))
        print('******* {} of {} *******'.format(i,len(input_dirs)))
        print('***********************************************************')
        main_run_phot(input_file, output_dir,
                      output_file_prefix = output_file_prefix,
                      show_plots = False,
                      save_plots = True,
                      save_photcat = True,
                      save_source_list = True)
        #now run on pipeline drizzle
        output_dir = '/Users/cshanahan/Desktop/WD_acs/output/pipeline_drizzle/'
        asn_id = fits.open(glob.glob(dirr + '/*flc.fits')[0])[0].header['ASN_TAB'][0:9]
        input_file = '/grp/hst/wfc3s/deustua/ACSparallels/{}_drc.fits'.format(asn_id)
        print('******* PROCESSING {} *******'.format(input_file))
        main_run_phot(input_file, output_dir,
                      output_file_prefix = output_file_prefix,
                      show_plots = False,
                      save_plots = True,
                      save_photcat = True,
                      save_source_list = False)
