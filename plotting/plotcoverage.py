""" plotcoverage.py
	
	Make a plot for the paper of how coverage varies with dither size.
	A similar older scrip exists god knows where.

	06/01/2016
	KMarkovic
"""

import numpy as np
import time as t
import os
from matplotlib import rc
import matplotlib.pyplot as plt
import plotpolys as pp
import coverage
import os.path as op
from glob import glob
from warnings import warn
from scipy.interpolate import interp1d
import plotlib as l

HERE = op.dirname(op.abspath(__file__)) + '/'
BASELINEX = l.GAPX # arcsec

# Set up default plotting fonts
font = {'family' : 'serif',
	    'weight' : 'normal',
    	'size'   : l.FS}
rc('font', **font)	

def sizes_of_pattern(pattern='J', nstep=1, outpath=op.join(HERE,'outputs/'), verb=False, binary=op.join(HERE,'bin/create-euclid-patch'), sizemax=2.0, sizemin=0.0):
	""" Do the calculation for different sizes. """
	sizes = sizemin + (sizemax - sizemin)/nstep*np.arange(nstep)
	covs = []; tots = []
	for size in sizes:
		dithvec = size*coverage.PATTERNS[pattern]
		outv, lim_rec = coverage.survey_coverage(dithvec, binary, outpath, verb)
		cov, tot = coverage.coverage(outv, lim_rec)
		covs.append(cov); tots.append(tot)
	
	if np.std(np.array(tots)) < 1e-15:
		warn("something is going wrong in the total area calculation!")

	return np.array(sizes), np.array(covs)/tot

def to_file(sizes, coverages, outpath, ts='', pattern=None, filename=None):
	""" Write the coverage percentages for each size to file. """

	if filename is not None: 
		openopt = 'a'
	else:
		openopt = 'w'

	writefile = op.join(outpath,'coverage' + ts +'.txt')
	with open(writefile, openopt) as f:
		if openopt=='w':
			if pattern is not None:
				f.write('# Coverage of the ' + pattern + ' pattern\n')
			else:
				f.write('# Coverage\n')
			f.write('# dither scale: 0-pass 1-pass 2-pass 3-pass 4-pass coverage area percentages\n')
		for size, cov in zip(sizes, coverages):
			f.write(str(size) + ': ' + '%, '.join([str(round(c,4)) for c in cov*100]) + '%\n')

	return writefile

def plot(sizes, coverages, plotpath=None, pattern=None, res=500):
	""" Plot the stacked bar plot of coverages (from old code: sometime pre-version-0.6). """

	# Make sure have Numpy arrays
	sizes = BASELINEX*np.array(sizes)
	coverages = np.array(coverages)
	nsiz, ncov = coverages.shape
	if res < nsiz: res = nsiz

	# First sort
	sort = np.argsort(sizes)
	old_sizes = sizes[sort]
	coverages = coverages[sort]

	# Then prepare to interpolate
	sizes = min(old_sizes) + (max(old_sizes)-min(old_sizes))*np.arange(res+1)/res
	# Then plot
	plt.figure(figsize=(9.3, 7.5))
	h = np.arange(ncov)
	width = (max(sizes) - min(sizes)  + 1e-15)/(len(sizes) - 1 + 1e-15)
	col = ['k','r','c','g','m','b','y']
	stack = 0.0 * sizes
	for hint in range(ncov)[::-1]:

		cov_interp = interp1d(old_sizes, coverages[:,hint], kind='linear')(sizes)

		if hint==0:
			plt.bar(sizes,cov_interp,width,bottom=stack,color=col[hint],label=str(hint)+'-pass',lw=0.,align='center')
		else:
			plt.bar(sizes,cov_interp,width,bottom=stack,color=col[hint],label=str(hint)+'-pass',lw=0.,align='center')
	
		stack += cov_interp

	if pattern is not None:
		xran = max(sizes) - min(sizes)
		plt.text(min(sizes) + xran*0.05, 0.05, coverage.convpat(pattern), family='serif', color='k', fontsize='45')

	plt.axvline(l.DX_TRANSITION, color='k',ls=':')
	plt.legend(bbox_to_anchor=(-0.12, 1.03, 0., 0.), loc="lower left", prop={'size':l.FS}, frameon=False, ncol=5, borderaxespad=0.0, columnspacing=0.5, handletextpad=0.25)
	plt.xlabel(r'step size, $d_x$ ["]'); 
	plt.ylabel('fractions of pixels with n-passes')
	plt.xlim([min(sizes)-0.01, max(sizes)+0.01])
	plt.ylim([0.0, 1.0])
	
	if plotpath is None:
		plotfile = None
		plt.show()
	else:
		plotfile = op.join(plotpath,"plotcoverage" + ts +'.eps')
		plt.savefig(plotfile,dpi=600,bbox_inches='tight')

	return plotfile

def from_file(filename, pattern=None):
	""" Read size-coverage info from file. Check pattern if given. """

	with open(filename, 'r') as f:
		content = f.readlines()
	
	# Check pattern if it is there
	if pattern is not None:
		if check_pattern(content[0], pattern): 
			raise Exception("You are trying to load a file that doesn't correspond to your selected pattern!")

	sizes, coverages = [], []
	for line in content[2:]:
		line = line[:-1].split(':')
		sizes += [float(line[0])]
		coverages.append([float(s[:-1]) for s in line[1].split('%, ')])

	return np.array(sizes), np.array(coverages)/100.0

def check_pattern(patline, pattern):
	""" Check the first line of the output file for the requested pattern. """
	patline = patline[10:-1]
	if len(patline)>1:
		match = not (patline[8:-8] == pattern)
	else:
		match = 2
	return match

if __name__=='__main__':

	import argparse
	p = argparse.ArgumentParser(description="Calculate n-pass coverage for different dither sizes and plot.")
	p.add_argument("--pattern", "-p", default="S", choices=["J", "R", "step", "S", "N", "X", "O", "box"], help="Name of dither pattern.")
	p.add_argument("--maxsize", "-s", default=2.0, type=float, help="Max scale of dither, size=1 means d_x=50'', d_y=100''.")
	p.add_argument("--minsize", "-m", default=0.0, type=float, help="Min scale of dither, size=1 means d_x=50'', d_y=100''.")
	p.add_argument("--nsizes", "-n", default=2, type=int, help="Number of steps between minsize and maxsize to calculate coverage for.")
	p.add_argument("-o", "--outpath", default=op.join(HERE,"../outputs/"), help="where you want your output files and plots.")
	p.add_argument("-c", "--srcpath", default=op.join(HERE,"../ubercal/"), help="directory containing test-calibration source code")
	p.add_argument("-b", "--bin", default=op.join(HERE,'../bin/'), help="Location of binaries if no need for compilation.")
	p.add_argument("-v", "--verb", default=0, type=int, help="verbosity level. 0 for silent.")
	p.add_argument("-d", "--display", action='store_true', help="Set flag to display plot instead of saving.")
	p.add_argument("-r", "--res", default=500, type=int, help="Plotting resolution to interpolate sizes to.")
	p.add_argument("-f", "--force", action='store_true', help="Force a new file.")
	arg = p.parse_args()

	# Since only calculating a 3x3 survey, need to warn if size of dither is too big
	if arg.maxsize > 2.0:
		warn("The coverage calculation is not guaranteed to work for very large dithers!")

	# Don't save plot if it should be displayed
	if arg.display:
		plotpath = None
	else:
		plotpath = arg.outpath

	if arg.pattern is "R":
		arg.pattern = "step"
	if arg.pattern is "O":
		arg.pattern = "box"

	# Make sure we have built the latest version of the c-code
	target = 'create-euclid-patch'
	binary = op.join(arg.bin,target)
	if not op.isfile(binary):
		binary = coverage.build(path=arg.srcpath, target=target, verb=arg.verb)

	# Search through the output path to see if this pattern has already been calculated
	files = glob(op.join(arg.outpath,'coverage*.txt'))
	exists = False
	writefile = None
	if not arg.force:
		for fn in files:
			with open(fn, 'r') as f:
				first_line = f.readline()
			if not check_pattern(first_line, arg.pattern):
				if exists is True:
					print "There is more than one file containing pattern " + arg.pattern + " coverages in the path:"
					print "\t" + arg.outpath
					print "I will make a new file for this calculation as I cannot choose between the existing ones to add to."
					exists = False
					writefile = None
					break
				writefile = fn
				exists = True


	# Pattern sizes range
	rangein = arg.maxsize-arg.minsize
	dsize = rangein/(arg.nsizes-1)
	sizemax = arg.maxsize
	sizemin = arg.minsize

	# If a file with the pattern exists, check whether the requested size range is full/partially/not included
	maxmatch = False; minmatch = False
	sizes = None
	if exists:
		# Open up the file containing all the coverages and read it
		sizes, coverages = from_file(writefile, arg.pattern)
		
		ts = writefile.split('/')[-1][8:-4]

		# Compare the extrema of the file with the requested extrema
		maxmatch = arg.maxsize <= max(sizes)
		minmatch = arg.minsize >= min(sizes)

		# If only one extremum is outside of the previous range, supplement the old writefile:
		if not maxmatch and minmatch:
			exists = False
			dsize = (arg.maxsize - max(sizes))/arg.nsizes
			sizemin = max(sizes) + dsize
			sizemax = sizemax
		elif not minmatch and maxmatch:
			dsize = 0
			exists = False
			sizemax = min(sizes)
		# If requested range is both below and above the calculated range, start afresh:
		elif not (maxmatch and minmatch):
			exists = False
			writefile = None

	# Print out what you end up doing:
	if exists:
		print "Full range already exists for pattern " + arg.pattern + " in file " + writefile + ","
		print "plotting."
	elif writefile is None and not (maxmatch and minmatch) and sizes is not None:
		print "The requested range goes below and above what exists in any of the files " + arg.pattern + "."
		print "Making a new file for pattern " + arg.pattern + " in the given path: \n\t" + arg.outpath + "."
	elif writefile is None:
		print "Making a new file for pattern " + arg.pattern + " in the given path: \n\t" + arg.outpath + "."
	else:
		print "Some of the requested range for pattern " + arg.pattern + " is already in the file: \n\t" + writefile + ","
		print "but we need to supplement that with calculations of sizes between " + str(sizemin-dsize) + " and " + str(sizemax) + "."

	#### End of checks, do the calculation & plotting
	if not exists:
		# If pattern-size doesn't yet exist, do the whole calculation
		if writefile is None: ts = '-' + t.strftime('%s')
		sizes, coverages = sizes_of_pattern(arg.pattern, arg.nsizes, arg.outpath, arg.verb, 
			binary, sizemax + dsize, sizemin)
		writefile = to_file(sizes, coverages, arg.outpath, ts, arg.pattern, filename=writefile)
		print 'Coverage calculated. Results saved in:\n\t' + writefile + ','

	# Now plot just the data from this calculation (not all from file -> do that separately)
	plotfile = plot(sizes, coverages, plotpath, arg.pattern, res=arg.res)

	# Concluding message and plot display/save
	if arg.display:
		print 'plot displayed.'
	else:
		print 'plot saved in: \n\t' + plotfile + '.'
