import glob

from astropy.coordinates import SkyCoord, match_coordinates_sky


def match_cats_common_sources(cats, n_min_frames):
    """Takes photomeyry catalogs containting sky positions (columns named
       'xcentroid' and 'ycentroid') of objects detected in the same field, or
       fields with some degree of overlap, and returns back each input catalog
       with only rows for those common sources."""


    for i, val1 in enumerate(cats):
        for j, val2 in enumerate(cats):
            if (i != j) & (i<j):
                match_idx, _, _ = match_coordinates_sky(cat, ref_cat)
