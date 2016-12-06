""" 
	For plotting the dither tests.
	
	17/03/2015, KMarkovic
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, pandas, glob
import coverage as c

FS = 20
matplotlib.rcParams.update({'font.size': FS})

MINY = 0.0
MAXY = 0.02025

# Get the baseline fcalb and icalb from the test.out file
BASELINE_x = 50.0
BASELINE_d = 311.8033989

DETX = 600 # arcsec
GAPX = 50 # arcsec 
GAPY = 100 # arcsec
NODETX = 4
NODITH = 4
MAXD = 4*np.sqrt((DETX+GAPX)**2+(DETX+GAPY)**2)
DX_TRANSITION = NODETX*(DETX+GAPY)/(NODITH-1)/2
#print "this: ", DX_TRANSITION

ALLPATS = ["baseline", "J", "R", "O", "S", "N", "X"]
C = {
	'J': 'k',
	'S': 'g',
	'R': 'c',
	'N': 'y',
	'O': 'm',
	'X': 'r'
	}
L = {
	'J': '-',
	'S': '-',
	'R': '--',
	'N': '--',
	'O': '-',
	'X': '--'
}
# Alternative pattern names for backward compatibility
C['box'] = C['O']
C['step'] = C['R']
L['box'] = L['O']
L['step'] = L['R']

def find_recent(path):
	""" Get the most recent run. """
	timestamp = 0
	runs = glob.glob(path+'*')
	for run in runs:
		fname = os.path.basename(run)
		irun = int(fname.split("-")[0])
		if irun > timestamp: 
			# Check run has finished
			if os.path.isfile(path + fname + '/run.log'):
				timestamp = irun
	return fname

def find_patterns(rundir, silence=False):
	""" Find all the patterns in rundir (should contain only 1 seed). """
	if not os.path.isdir(rundir): raise Exception("I found no such folder: " + rundir)

	pattern = []
	for fnm in glob.glob(rundir + "/*.dat"):
		pattern.append(os.path.basename(fnm).strip(".dat")) # Get the patterns available

	if not silence: print "I found these patterns: ", pattern, " in ", rundir
	return pattern	

def get_run(outpath, fname=0):

	# Find most recent run
	if not fname: fname = find_recent(outpath)
	print 'reading ' + str(fname)

	# Find all the patterns in the timestamp folder
	rundir = os.path.join(outpath,str(fname))+'/'
	patterns = find_patterns(rundir)

	return rundir, patterns, str(fname)

def get_runs(outpath, rnames=0):
	# Get the means and stdevs of several runs (hpefully many files in a dir)

	print "In " + outpath + ", scanning:"
	if not rnames: rnames = glob.glob(outpath + '*')
	if rnames == []: raise Exception("I found nothing in "+outpath)
	
	patterns = ALLPATS[1:]
	dfs = {}
	tmp = []

	for rname in rnames:	
		run = rname.split('/')[-1]
		runs = {}

		pat = None
		for fnm in glob.glob(rname + "/*.dat"):
			pat = os.path.basename(fnm).strip(".dat")

			if pat in patterns:
				runs[pat] = load_pattern(fnm, recalc=False)
				tmp.append(pat)
		if pat is None: continue
		patterns = list(set(patterns) & set(tmp))
		dfs[run] = runs
	print "Only these patterns found in all folders: " + str(patterns)

	return dfs, patterns

def get_baseline(rundir):
	tmp = load_pattern(rundir + '/baseline.dat')
	return tmp.fcal, tmp.ical

def load_pattern(fname, recalc=True):
	""" Load the pattern from a [].dat file. Recalculate the calibrations from the calibrations.txt
		files if recalc is requested. """

	# Don't recalculate in the exp-to-exp case - it will be wrong!
	if recalc and ('exp' in fname or 'EXP' in fname):
		print "Not recalculating in the exp-to-exp case!"
		recalc = False

	# Get both headers
	with open(fname) as f:
		head = ''.join([next(f)[2:] for i in [0,1]])[:-1]

	# This is the header format that testdith.py saves in
	fhead = ['i','d','area','ical','fcal']
	fhead += ['x_1','y_1','x_2','y_2','x_3','y_3']
	fhead += ['no_pointings_x','no_pointins_y']
	fhead += ['frac_0_dith','frac_1_dith','frac_2_dith','frac_3_dith','frac_4_dith']

	# Construct a new filename for the recalculation result and check if it has already been done
	if recalc:
		recalced_fname = fname[:-4]+'-recalculated.dat'
		if os.path.isfile(recalced_fname): 
			print "Recalculation found in " + recalced_fname + "!"
			fname = recalced_fname
			recalc = False

	# Grab the data from the .dat file
	df = pandas.read_table(fname, sep=' ', names=fhead, skiprows=2)

	# Re-calculate the calibration stats if it is requested and not found yet
	if recalc:

			# Replace the columns with the new calculations
			dx, df['ical'], df['fcal'] = calibrations_stats(fname[:-4])

			# Save into a new file and return
			np.savetxt(recalced_fname, df, fmt='%.6f', header=head)

	return df

def calibrations_stats(path, silence=False):
	""" Call calibration_stats for all the dither configurations in the pattern. """

	# Baseline has no folder:
	if 'baseline' in path:
		baseline = True
		path = os.path.dirname(path)
		paths = glob.glob(path+'/calibrations.txt')
	else:
		baseline = False
		paths = glob.glob(path+'/*/calibrations.txt')

	# Find all the calibrations.txt files in the given path
	if len(paths)==0: 
		raise Exception('No calibrations.txt files found in '+ path + '!')

	sigma_i = np.ones(len(paths))*np.inf
	sigma_f = np.ones(len(paths))*np.inf
	dithersizes = np.zeros(len(paths))
	for i,file in enumerate(paths):

		# Grab the dither size from the filename
		if baseline:
			dithersizes[i] = BASELINE_x
		else:
			dithersizes[i] = float(file.split('/')[-2][11:])
		if not silence: print 'recalculating dx = ' + str(dithersizes[i])
		
		# Get the fov configuration from the header
		with open(file) as f:
			headline = next(f)
		N1, N2 = headline[39:-10].split('x')
		N = int(N1)*int(N2)

		# The columns are:
		# exposure-id initial-zero-point mean-measured-flux used-area calibration-correction calibrated-zero-point
		data = np.loadtxt(file)

		if 'exp' not in headline:
			# Extract the sigma_i, sigma_f from the given calibrations.txt file.
			# Loop over fovs
			nexp = int(len(data)/N)
			newdata = np.zeros([nexp,2])

			# This loop is very slow
			for alpha in range(nexp):
				newdata[alpha,0] = np.mean(data[alpha*N:(alpha+1)*N,1])
				newdata[alpha,1] = np.mean(data[alpha*N:(alpha+1)*N,5])

			sigma_i[i] = np.std(newdata[:,0]) 
			sigma_f[i] = np.std(newdata[:,1])
		else:
			sigma_i[i] = np.std(data[:,1]) 
			sigma_f[i] = np.std(data[:,5])

	# Sort in increasing dithersize
	index = np.argsort(dithersizes)

	return dithersizes[index], sigma_i[index], sigma_f[index]

def plot_vs_x(rundir, patterns=["J", "O", "S", "R"], plotssize = [None], linestyles = ['-', '--'], vsd=False):

	# Whether to plot against dx or total d:
	if vsd:
		baseline = BASELINE_d
	else:
		baseline = BASELINE_x

	# Open all the pattern.dat (except baseline)
	maxx=0; minx=float('inf')
	for p in patterns:
		if p=='baseline': continue

		ts = rundir.strip('/').split('/')[-1]
		fname = os.path.join(rundir, p+'.dat')
		if not os.path.isfile(fname): 
			print('Can find the file for ' + p + ' in the ' + ts + ' run. Skip it.')
			patterns.remove(p)
			continue
		tmp = load_pattern(fname)

		if vsd:
			ordinate = tmp.d
		else:
			ordinate = tmp.x_1

		# Find limits for x-axis
		if max(ordinate)>maxx: maxx = max(ordinate)
		if min(ordinate)<minx: minx = min(ordinate)

		# Plot line for each pattern
		for npoint, style in zip(plotssize,linestyles):
			if npoint is not None: 
				test = tmp.no_pointings_x==npoint
				if np.sum(test)==0: 
					print "No runs with "+str(int(npoint))+" pointings found in "+rundir+"!",
					npoint = min(tmp.no_pointings_x-npoint)+npoint
					test = tmp.no_pointings_x == npoint
					print "Plotting npoint="+str(int(npoint))+" instead."
			else:
				test = range(len(tmp))
			plt.plot(ordinate[test], np.array(tmp.fcal[test]), style, 
					 label = c.convpat(p) + r': '+str(int(npoint))+'x'+str(int(npoint)), c=C[p], lw=2); 

	plt.xlabel(r'step size, $d_x$ ["]'); 
	if not vsd: plt.axvline(DX_TRANSITION, color='k', ls=':')
	if 'baseline' in patterns:
		baseline_y, init_y = get_baseline(rundir)
		plt.axvline(baseline, ymin=0, ymax=(baseline_y-MINY)/(MAXY-MINY), ls='--', c=C['J'], lw=2, zorder=len(patterns))
		plt.scatter(baseline, baseline_y, marker = "o", s = 50, c=C['J'], label="Laureijs et al. (2011)", zorder=len(patterns))
	plt.ylabel(r'final zero-point scatter, $\sigma_f$'); 
	plt.legend(scatterpoints=1,fontsize=FS, loc=1, ncol=1, handletextpad=0, labelspacing=0.25);
	if vsd and max(tmp.d)>MAXD/4.0: maxx = MAXD/4.0
	plt.xlim([minx, maxx])
	plt.ylim([MINY, MAXY])

def plot_vs_d(rundir, patterns=["J", "O", "S", "R"], plotssize = [20, 19], linestyles = ['-', '--']):
	""" Only for backward compatibility. """
	plot_vs_x(rundir, patterns, plotssize, linestyles, vsd=True)

def plot_vs_nx(stats, mask, patterns=ALLPATS, lw=2):
	""" Plot vs survey size """

	# Pull out the stats
	x, mu_f, sig_f, d, mu_i, sig_i = stats

	# Open the figure
	fig = plt.figure()
	
	# Loop through the patterns in the table
	leg = []; lab = ''
	MINY=0.0; MAXY=0.0
	for pat in patterns:

		# If a non-pattern file snuck in - ignore it
		if pat not in mu_f.keys(): continue

		# Get the mean of the sigma... eek. It's just for the plot.
		wmean_i = np.average(mu_i[pat][mask],weights=np.power(sig_i[pat][mask],-2))
		wmean_f = np.average(mu_f[pat][mask],weights=np.power(sig_f[pat][mask],-2))
		print "The mean residual zero-point scatter for " + pat + "-pattern is " + str(wmean_f) + "."

		# Make the plot against the root(number of pointings)
		plt.plot(x[mask], wmean_i+x[mask]*0, c='0.75', ls='--')
		plt.plot(x[mask], wmean_f+x[mask]*0, c=C[pat], ls=':', label=None)
		p, = plt.plot(x[mask], mu_f[pat][mask], c=C[pat], ls=L[pat], label=c.convpat(pat) + ": d = " + str(round(d[pat],1)), lw=lw)

		# For prettification
		leg.append(p)
		plt.fill_between(x[mask], mu_f[pat][mask]-sig_f[pat][mask], mu_f[pat][mask]+sig_f[pat][mask], interpolate=True, facecolor='0.9', edgecolor='0.9')
		MAXY = max(MAXY, np.max(wmean_i))

	# Prettyfy
	plt.xlim([1, max(x)])
	plt.ylim([MINY, MAXY*1.001])
	plt.ylabel(r'final zero-point scatter, $\sigma_f$') 
	plt.xlabel(r'survey size, $\sqrt{N}$-pointings')
	plt.legend(handles=leg, ncol=2, loc=1, fontsize=FS);

	return fig, lab

def save_or_show(seefig, fname=None, fig=None):

	if fig==None: fig = plt.gcf()
	
	# Save or show plot
	if seefig:
		plt.show();
	else:
		plt.savefig(fname,dpi=400,bbox_inches='tight')
		print "Plot saved to: " + fname
	plt.close()

def get_stats(tables, dx=50.0):

	# Assuming patters same in tables (get_tables gives it):
	patterns = tables[tables.keys()[0]].keys()

	# Assuming all are same length, x-axis!
	test = (tables[tables.keys()[0]][patterns[0]].x_1 - dx)**2 < 0.01 # only precise to 1dp!!
	if not sum(test): raise Exception("No dx = " +str(dx)+ " found!")
	nx_vec = tables[tables.keys()[0]][patterns[0]].no_pointings_x[test]  
	stdevs_i = {}
	stdevs_f = {}
	means_i = {}
	means_f = {}
	td_vec={}
	N = len(tables)

	for pat in patterns:

		print pat + ": ",
		mean_i = np.zeros(len(nx_vec))
		mean_f = mean_i*0.0
		sd_V = mean_i*0.0
		tmp = mean_i*0.0
		sig_isq = mean_i*0.0
		sig_fsq = mean_i*0.0
	
		for key in tables.keys():

			print key,

			table = tables[key][pat]
			value_i = np.array(table.ical[test])
			value_f = np.array(table.fcal[test])

			mean_i += value_i/N
			mean_f += value_f/N

		for key in tables.keys():

			table = tables[key][pat]
			value_i = np.array(table.ical[test])
			value_f = np.array(table.fcal[test])

			sig_isq += (value_i - mean_i)**2/N
			sig_fsq += (value_f - mean_f)**2/N

		stdevs_i[pat] = np.sqrt(sig_isq)
		stdevs_f[pat] = np.sqrt(sig_fsq)
		means_i[pat] = mean_i
		means_f[pat] = mean_f
		td_vec[pat] = np.array(table.d[test])[0]

		print ''

	return nx_vec, means_f, stdevs_f, td_vec, means_i, stdevs_i

def find_seedfolders(inpath, tss=[]):
	""" Find the folders that contain a seed each with all the patterns in them."""

	# Find all the right directories
	# http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
	ignoredirs = ['.git', 'mangle', 'bin', 'plotting', 'samples'] + ALLPATS
	foundtss = []
	for root, dirnames, filenames in os.walk(inpath):
		if 'run.log' in filenames: foundtss.append(root)
		for igdir in ignoredirs:
			if igdir in dirnames: dirnames.remove(igdir)
	if len(foundtss)==0: raise Exception('No ubercal run results found in ' + inpath + '!')

	# If a list of timestamps is given to plot, make sure we've found all of them and save paths
	if len(tss)>0:
		tspaths=[]
		for ts in tss:
			for fts in foundtss:
				if ts in fts:
					tspaths.append(fts)
					tss.remove(ts)
		if len(tss)>0: print "Warning: we didn't find the following timestamps in the path you gave: " + ', '.join(tss)
	else:
		tspaths = foundtss

	return tspaths
