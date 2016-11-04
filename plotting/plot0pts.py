""" This is to test the before and after ubercal calibrations.
	Primarily to figure out why we are not getting the factor of 2.

	03/10/2016 Markovic
"""

import matplotlib.pyplot as plt
import numpy as np
import re

def get_survey_config(overlapfile):
	""" Extract from first two lines of ubercal overlap file. """

	with open(overlapfile, 'r') as f:
		ndith, ndet = [int(s) for s in re.split('; |, |=|\s', f.readline()) if s.isdigit()]
		nra, ndec = [int(s) for s in re.split('; |, |=|\s', f.readline()) if s.isdigit()]

	return ndith, ndet, nra, ndec

def vertical_read(indir):
	""" Since J pattern is more vertical, display the calibrations grouped by vertical alighnment."""

	# Get the survey config so we know how the data is arranged
	ndith, ndet, nra, ndec = get_survey_config(indir+'/full-survey-overlaps.txt')

	# Get the data
	array = np.loadtxt(indir+'/calibrations.txt')

	# Check if it is det-to-det or exp-to-exp
	if ndith*nra*ndec == len(array): ndet=1

	# Rearrange one dither layer at a time first
	newarray = np.zeros(array.shape)
	for c in range(ndith):
		
		chunk = array[c*len(array)/ndith:(c+1)*len(array)/ndith]
		full_dither = c*ndet*ndet*ndec*nra

		# Find the i-th column across all te vertically alighned exposures
		for i in range(ndet*nra):

			# Get the j-th element in the i-th column
			for j in range(ndet*ndec):
				
				# which exposure in this dither
				exposure_row = j/ndet
				exposure_col = i/ndet
				exposure_number = nra*exposure_row + exposure_col
				
				# which detector in this exposure
				detector_row = j%ndet
				detector_col = i%ndet
				detector_number = detector_row + ndet*detector_col

				# put together
				total_index = ndet**2*exposure_number + detector_number

				# assign element
				newarray[full_dither + i*ndet*ndec + j] = chunk[total_index]

	# Now ungroup the dithers too so that it is all vertically alightned
	return revectorise(newarray, ndith)

def revectorise(array, blocksize=4):
	""" Rearrange the rows in the array by creating blocks of size blocksize. """

	nblocks = len(array)/blocksize
	assert nblocks*blocksize == len(array), "Size of array not divisible by given blocksize!"

	newarray = np.zeros((array.shape[0],array.shape[1]+1))
	for i in range(nblocks):
		for j in range(blocksize):
			newarray[i*blocksize + j,:-1] = array[j*nblocks + i]

	# Reindex
	newarray[:,-1] = newarray[:,0]
	newarray[:,0] = range(len(newarray))

	return newarray

if __name__=="__main__":

	import argparse, time, os.path

	parser = argparse.ArgumentParser(description="Plot VERB>1 output files of ubercal.")
	parser.add_argument("cals", help="folder containing the bestfit calibrations (calibrations.txt)")
	parser.add_argument("-s", "--save", default=False, action='store_true', help="true for saving into given folder, default: show (i.e. false)")
	args = parser.parse_args()

	# Which column we expect each value to be
	ID=0; INIT=1; FLUX=2; CALS=4; FINAL=5; ORIG_ID=6;

	# Read from file & reorder so that all the exposures of each detector are adjacent 
	if not os.path.isdir(args.cals): raise Exception(args.cals + ' is not a folder!')
	cals = vertical_read(args.cals)

	# Make the plot
	fig, ax = plt.subplots()
	# Awkward gridlines
	for i in range(len(cals)-1):
		plt.axvline(x=i+0.5, color='0.8', linestyle='-', linewidth=1)

	plt.bar(cals[:,0]-0.5, cals[:,INIT], fill=True, width=1, label='initial 0-points', color='0.5', ls='--', lw=0)
	plt.bar(cals[:,0]-0.5, cals[:,FLUX], fill=True, width=1, label='calibrator fluxes', color='g', alpha=0.5, ls='-', lw=0, bottom=cals[:,1])
	plt.axhline(y=np.mean(cals[:,INIT]), ls='--', label='initial mean')
	plt.axhline(y=np.mean(cals[:,INIT])+np.mean(cals[:,2]), ls='-', color='y', label='initial tot mean', lw=1)
	plt.axhline(y=np.mean(cals[:,FINAL]), ls=':', label='final mean')
	plt.plot(cals[:,0], cals[:,FINAL], 'r.', label='calibrated 0-points')
	plt.plot(cals[:,0], cals[:,INIT]+cals[:,CALS]+cals[:,FLUX], color='m', marker='+', label='calibrated fluxes', lw=0)

	# Annotate
	plt.xlim([cals[0,ID]-0.5,cals[-1,ID]+0.5])
	plt.xticks(cals[:,ID])
	ax.tick_params('both', length=0, width=2, which='major')
	ax.set_xticklabels(cals[:,ORIG_ID].astype(int), fontsize=10)
	plt.xlabel('original detector number')
	plt.ylabel('zero points')
	plt.legend(bbox_to_anchor=(1, 1), loc='lower right', ncol=2)
	plt.title(' & '.join(args.cals.split('/')[-3:-1]), x=0.05)

	# Show or save
	if not args.save:
		plt.show();
	else:
		ts = time.strftime('%s')
		plt.savefig(args.cals+"/plot0pts-" + ts +'.pdf',dpi=400,bbox_inches='tight') 
		plt.close()