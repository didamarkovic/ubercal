""" This is to test the before and after ubercal calibrations.
	Primarily to figure out why we are not getting the factor of 2.

	03/10/2016 Markovic
"""

import matplotlib.pyplot as plt
import numpy as np

NDITH = 4

def revectorise(array, blocksize=4):
	""" Rearrange the rows in the array by creating blocks of size blocksize. """

	nblocks = len(array)/blocksize
	assert nblocks*blocksize == len(array), "Size of array not divisible by given blocksize!"

	newarray = np.zeros(array.shape)
	for i in range(nblocks):
		for j in range(blocksize):
			newarray[i*blocksize + j] = array[j*nblocks + i]

	return newarray

if __name__=="__main__":

	import argparse, time, os.path

	parser = argparse.ArgumentParser(description="Plot VERB>1 output files of ubercal.")
	parser.add_argument("cals", help="folder containing the bestfit calibrations (calibrations.txt)")
	parser.add_argument("-s", "--save", default=False, action='store_true', help="true for saving into given folder, default: show (i.e. false)")
	args = parser.parse_args()

	# Read from file & reorder so that all the exposures of each detector are adjacent 
	if not os.path.isdir(args.cals): raise Exception(args.cals + ' is not a folder!')
	cals = revectorise(np.loadtxt(args.cals+'/calibrations.txt'), blocksize=NDITH)
	cals = revectorise(cals, blocksize=int(np.sqrt(len(cals)/NDITH)))
	# re-index too
	cals[:,0] = np.sort(cals[:,0])

	# Make the plot
	plt.bar(cals[:,0]-0.5, cals[:,1], fill=True, width=1, label='initial 0-points', color='0.5', ls='--', lw=0)
	plt.bar(cals[:,0]-0.5, cals[:,2], fill=True, width=1, label='calibrator fluxes', color='g', alpha=0.5, ls='-', lw=0, bottom=cals[:,1])
	plt.axhline(y=np.mean(cals[:,1]), ls='--', label='initial mean')
	plt.axhline(y=np.mean(cals[:,1])+np.mean(cals[:,2]), ls='-', color='y', label='initial tot mean', lw=1)
	plt.axhline(y=np.mean(cals[:,-1]), ls=':', label='final mean')
	plt.plot(cals[:,0], cals[:,-1], 'r.', label='calibrated 0-points')
	plt.plot(cals[:,0], cals[:,-1]+cals[:,2], color='m', marker='+', label='calibrated fluxes', lw=0)

	# Annotate
	plt.xlabel('detector number')
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