from astropy.io import fits
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
import glob
import os
from stwcs import updatewcs

"""Drizzles images taken in the same filter/visit together."""

def run_updatewcs(input_file):
	"""Check if update WCS has been run already. If not, run it on the file.
	   Leaving this blank for now because images probably don't need to be aligned."""
	pass

def align_images(input_files, ref_image):
	"""Leaving this blank for now. Images are grouped by visit and should be aligned."""
	pass

def construct_driz_input_dict(output, final_bits, driz_sep_bits, 
							  driz_cr_corr, driz_cr_grow, driz_cr_scale, driz_cr_snr):
	
	input_dict = {	'output' : output,
					'final_bits' : final_bits,
					'driz_sep_bits' : driz_sep_bits,
					'driz_cr_corr' : driz_cr_corr,
					'driz_cr_grow' : driz_cr_grow,
					'driz_cr_scale' : driz_cr_scale,
					'driz_cr_snr' : driz_cr_snr,
					'clean' : True,
					'in_memory' : True,
					'configobj' : None,
					'build' : True}
					
	return input_dict
	
def write_param_file_to_directory(input_dict, input_file_directory):
	pass

def run_driz(input_file_paths, input_dict):
	
	print('Drizzling files {}'.format(input_file_paths))
	
	astrodrizzle.AstroDrizzle(input_files, **input_dict)
	
def main_full_run_driz(input_files, final_bits, driz_sep_bits, driz_cr_corr, 
					   driz_cr_grow, driz_cr_scale, driz_cr_snr, align = False, 
					   write_param_file = True):

	base_dir = input_files[0].replace(os.path.basename(input_files[0]),'')
	filter = base_dir.split('/')[-2]
	output = filter+'_combined'
	
	start_dir = os.getcwd()
	os.chdir(base_dir)
	
	input_files = [os.path.basename(f) for f in input_files]
	driz_input_dict = construct_driz_input_dict(output, final_bits, 
												driz_sep_bits, driz_cr_corr, 
												driz_cr_grow, driz_cr_scale, 
												driz_cr_snr)
	if write_param_file:
		with open(base_dir + output + '_params.txt', 'w') as f:
			print('Writing {}.'.format(base_dir + output + '_params.txt'))
			for key in driz_input_dict:
				f.write('{}: {}\n'.format(key, driz_input_dict[key]))
	run_driz(input_files, driz_input_dict)
		
	if align:
		for f in input_file_paths:
				run_updatewcs(f)
		align_images(input_files[1:], input_files[0]) #use 0th as ref
	
	os.chdir(start_dir)
if __name__ == '__main__':
	
	input_dir = '/Users/cshanahan/Desktop/WD_acs/data/cmd_p5/F775W/'
	input_files = glob.glob(input_dir + '*flc.fits') #must be all from same visit/filter

	### set params
	final_bits = 96
	driz_sep_bits = 96
	driz_cr_corr = True
	driz_cr_grow = 1
	driz_cr_scale = (1.2,0.7)
	driz_cr_snr = (3.5,3.0)
	
	main_full_run_driz(input_files, final_bits, driz_sep_bits, driz_cr_corr, 
					   driz_cr_grow, driz_cr_scale, driz_cr_snr, align = False, 
					   write_param_file = True)
			

