This repository contains scripts to process ACS parallel fields (drc.fits) files and make photometry catalogs of both point and extended sources in the fields. The following outlines what is done to produce the final photometry catalogs and what they contain. 

1. (sort_data_visit_filter.py) Data (drc files) are sorted in subdirectories based on program/visit and then again by filter. The program/visit ID is taken from the file name. 

2. (make_source_catalog.py) Source detection is run on each file. This process consists of:
	1. Making a 2D background image (I did this in 100x100 pixel segments)
	2. Using this background image to make a pixel-by-pixel threshold image for 	source detection. This was necessary because the background levels and S/N 	vary greatly depending on how many FLCs were drizzled in that particular region of the image, so this works better than using a global threshold value for the 	whole image.
	3: Create a segmentation map of the image. The data is smoothed with a 	gaussian kernel (3x3 pixel kernel size, fwhm of 2.5 pixels). The number of 	connected pixels required to count as a source is also set - I found that 75 pixels does a good job. I used a stricter criteria of 150 pixels for images where there was only one input image, to avoid some of the cosmic rays. 
	4. Sources are deblended (photutils.deblend_sources). The fields aren’t crowded so this just takes care of a few close by sources. 
	5. I use the context image to make an array that represents, pixel by pixel, the number of exposures. Using this image, detected sources that fall close to regions of the image that have 0 exposure coverage (like chip gap due to large dither), are removed from the catalog. 
3. (make_cr_map.py). Run LAcosmic on each image to create a cosmic ray mask (not actually correct the data). Next, convolves this mask with a gaussian kernel to ‘grow’ the cosmic rays out a few pixels. This will create <filename>.mask.fits images. This CR mask will later be used to reject sources that are in regions of n_exp=1 that have a detected CR nearby. 

4. (align_imgs_gaia.py) Download gaia catalog near field from astroquery and align images to gaia with tweakreg. They shouldn’t be misaligned greatly, but this will allow final source positions to be in the gaia frame and for gaia catalog sources to be able to be matched to these catalogs. 

4. (make_phot_table.py). Runs photutils.aperture_photometry on each image (which now has a wcs aligned to gaia). The source catalog generated in step 2 and CR mask generate in step 3 are used. 
	1. Metadata from the file header is obtained. 
	2. Aperture photometry on each source in the input source catalog. Aperture radii 	of 2, 4, 6, 8 and 10 pixels are used (this can be set).
	3. Background level and rms is calculated for each source. I used the sigma-	clipped median in a circular annulus (r_in = 10pix, r_out = 15pix) around the 	source.
	4. Photometric errors are calculated for each source. See make_phot_table._compute_phot_err_daophot for more info. 
	5. A radial profile is fit to each source, and the FWHM of this fit is returned.
	6. The CR mask and context image are used to determine which sources fall in 	areas of only one exposure of coverage AND have a detected CR. These 	sources are given a flag cr_flag=1 in the table to signify they may be a CR.
	7. Sources in the detected source catalog are matched to the downloaded gaia source catalog (astropy.match_catalogs). If there is a matching source, the gaia 	RA and dec will be recorded. If there is no match, these columns 	(GAIA_CAT_RA, 	GAIA_CAT_DEC) are set to -999.
	8. Table is written out.

The photometry table produced in step 4 has (x, y) source positions, countrates and countrate error for each photometric aperture for each source (non background subtracted), median sky and std in the sky (countrate/pixel), as well as RA and Dec of each source in the image (column names RA_GAIA_WCS, DEC_GAIA_WCS since the image wcs has been tweaked to align to gaia), and finally the coordinates of any matching gaia sources (GAIA_CAT_RA and GAIA_CAT_DEC). A subset of this catalog the string ‘_noCRs’ appended to the file name is also written out - this catalog only contains sources with cr_flag = 0 (see step 4.6).

The final processing step - to take these catalogs with countrates and subtract sky level, apply photflam, and do aperture corrections exists in a notebook (finalize_cats.ipynb). I will eventually turn this into a standalone script, but it currently exists in a messy state in this notebook. Once this is run, the final catalogs will contain the following columns:
	x_center, y_center: X, Y position of the source in the drc file.
	flux_<2,4,6,8,10>_apcor : flux of each source in the various aperture sizes. To calculate this, I did:
 flux = ( (countate - (background * area) ) / aperture_correction) * photflam 
	err<2,4,6,8,10>_apcor : photometric errors. I did the same calculation above to the countrate errors. 
	RA_GAIA_WCS, DEC_GAIA_WCS : Ra and Dec for each source in gaia frame.
	GAIA_CAT_RA, GAIA_CAT_DEC: Exact coordinates in gaia catalog of matching source, if there was one. -999 otherwise.

Email me (cshanahan@stsci.edu) if you have any questions!



