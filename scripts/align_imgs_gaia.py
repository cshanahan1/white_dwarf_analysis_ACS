import os
import sys
import glob

for dirrr in glob.glob('/Users/cshanahan/Desktop/WD_acs/data/*/*/'):
    print(dirrr)
    os.chdir(dirrr)
    print(glob.glob('*drc.fits'))
    if os.path.isfile(dirrr + 'tweakreg.log'):
        print('skipping', dirrr)
    else:
        cmd = "python /Users/cshanahan/Desktop/WD_acs/scripts/current_scripts/pabeta/driz_tools/reg_drzs.py -g -i '*drc.fits'"
        os.system(cmd)
