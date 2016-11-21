""" Plot the likelihoods for diagnostics.

	03/10/2016 Markovic
""" 

import matplotlib.pyplot as plt
import numpy as np

def prob(chisq):
	""" Return an unnormalised probability given the chisq. """
	return np.exp(-chisq/2)

if __name__=='__main__':

	import argparse, time, os.path

	parser = argparse.ArgumentParser(description="Plot VERB>1 likelihood output file of ubercal.")
	parser.add_argument("likes", help="folder containing the likelihoods (likelihoods.txt)")
	parser.add_argument("-s", "--save", default=False, action='store_true', help="true for saving into given folder, default: show (i.e. false)")
	args = parser.parse_args()

	# Read from file & reorder so that all the exposures of each detector are adjacent 
	if not os.path.isdir(args.likes): raise Exception(args.likes + ' is not a folder!')
	tries = np.loadtxt(args.likes+'/likelihoods.txt')
	chisq = tries[:-1,-2]
	chisq0 = tries[:-1,-1]
	truth = tries[-1,:-2]
	tries = tries[:-1,:-2]

	# Read best fit from other file
	cals = np.loadtxt(args.likes+'/calibrations.txt')

	# Make the plot
	cls = 'kbcgyrm'
	for i, tru in enumerate(truth):
		if i>6: break
		if i==0:
			plt.plot(tries[:,i], prob(chisq), cls[i]+'.', label='attempts')
			plt.axvline(x=tru, ls=':', c=cls[i], label='truth')
			plt.axvline(x=cals[i,2], ls='--', c=cls[i], label='best-fit')
		else:
			plt.plot(tries[:,i], prob(chisq), cls[i]+'.')
			plt.axvline(x=tru, ls=':', c=cls[i])
			plt.axvline(x=cals[i,2], ls='--', c=cls[i])

	# Annotate
	plt.xlabel('calibration')
	plt.ylabel(r'exp($-\chi^2/2$)')
	plt.legend()
	#plt.title(' & '.join(args.likes.split('/')[-3:-1]))

	# Show or save
	if not args.save:
		plt.show();
	else:
		ts = time.strftime('%s')
		plt.savefig(args.likes+"/plotprobs-" + ts +'.eps',dpi=600,bbox_inches='tight') 
		plt.close()