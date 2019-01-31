import shutil

from astropy.io import fits
import numpy as np

def make_PAMcorr_image(image, outfile='default', extension=0):
    """Creates the Pixel Area Map (PAM) image.

    Parameters:
        image : string
            Name of FITS file.
        outfile : string
            Name of outfile.
        extension : int
            Extension=0 (UV images) by default.

    Returns:
        outfile : FITS file
            The pixel area map called ``<image rootname>_PAM.fits``

    Outputs:
        nothing
    """

    pamdir = '/Users/cshanahan/Desktop/'

    # -- Parse output filename & save a copy to file (NOTE: if outfile == input image, data is overwritten).
    if (image != outfile):
        if outfile == 'default':
            outfile = image.split('.fits')[0] + '_PAM.fits'
        shutil.copy(image,outfile)
    else:
        print('OVERWRITING DATA FOR IMAGE: '+image+'.')

    # -- Read in fits image and header (assume flt/flc - should handle both full- & sub-arrays)
    prihdr = fits.getheader(outfile)
    fdata = fits.getdata(outfile, ext=extension)   # ext=1 for IR?
    exptime = prihdr['exptime']
    detector = prihdr['detector']

    # -- Cycle through each SCI extension
    hdulist = fits.open(outfile,mode='update')
    for ff in range(len(hdulist)):
        if hdulist[ff].name == 'SCI':
            # -- read in header and data info
            scihdr = hdulist[ff].header
            data = hdulist[ff].data
            if detector == 'IR':
                chip = 1
            elif detector == 'WFC':
                chip = scihdr['CCDCHIP']
            else:
                raise Exception('Detector '+detector+' not covered in our case list.')

            naxis1 = scihdr['NAXIS1']
            naxis2 = scihdr['NAXIS2']
            x0 = int(np.abs(scihdr['LTV1']))
            y0 = int(np.abs(scihdr['LTV2']))
            x1 = int(x0 + naxis1)
            y1 = int(y0 + naxis2)

            # -- apply the PAM
            if detector == 'WFC':
                if chip == 1:
                    pam=fits.getdata(pamdir+'wfc1_pam.fits')

                    hdulist[ff].data = data * pam[y0:y1,x0:x1]

                elif chip == 2:
                    pam=fits.getdata(pamdir+'wfc2_pam.fits')

                    hdulist[ff].data = data * pam[y0:y1,x0:x1]

                else:
                    raise Exception('Chip case not handled.')
                    hdulist.close()
            elif detector == 'IR':
                pam=fits.getdata(pamdir+'ir_wfc3_map.fits')
                hdulist[ff].data = data * pam[y0:y1,x0:x1]
            else:
                raise Exception('Detector '+detector+' not covered in our case list.')
                hdulist.close()
    hdulist.close()

    return outfile
