""" MUST BE RUN IN PYTHON 2.7. """
import glob

from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
import cosmics
from multiprocessing import Pool

def run_lacosmic(data):
    c = cosmics.cosmicsimage(data, gain=1.5, readnoise=1, sigclip = 1.5,
                          sigfrac = 0.2, objlim = 1, verbose=False)
    c.run(maxiter=5, verbose=False)
    mask = c.mask.T

    return mask

def grow_cr(mask):
    kernel = Gaussian2DKernel(0.5)
    conv_mask = convolve(mask, kernel)
    conv_mask = conv_mask.astype(bool).astype(int)

    return conv_mask

def output_mask(ifile, conv_mask):
    hdr = fits.open(ifile)[1].header
    cosmics.tofits(ifile.replace('.fits','.convmask.fits'), conv_mask, hdr)

def main_make_cr_map(ifile):

    data = fits.open(ifile)[1].data
    mask = run_lacosmic(data)
    conv_mask = grow_cr(mask)
    output_mask(ifile, conv_mask)


if __name__ == '__main__':

    ifiles = glob.glob('/Users/cshanahan/Desktop/WD_acs/data/*/*/*drc.fits')
    p = Pool(6)
    p.map(main_make_cr_map, ifiles)
