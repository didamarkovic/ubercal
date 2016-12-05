""" 
	For plotting the dither tests.
	
	17/03/2015, KMarkovic
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os.path as op
from os.path import basename, isfile, isdir
from glob import glob
import pandas as pd
import coverage as c

FS = 20
matplotlib.rcParams.update({'font.size': FS})

MINY = 0.0
MAXY = 14.5

# Get the baseline fcalb and icalb from the test.out file
BASELINE_x = 50.0
BASELINE_d = 311.8033989

DETX = 612 # arcsec
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

def get_run(outpath, fname=0):

	# Get the most recent run
	if not fname:
		timestamp = 0
		runs = glob(outpath+'*')
		print outpath
		for run in runs:
			#try:
			fname = basename(run)
			irun = int(fname.split("-")[0])
			if irun > timestamp: 
				# Check run has finished
				if isfile(outpath + fname + '/run.log'):
					timestamp = irun
			#except ValueError:
			#	pass
		timestamp = str(timestamp)
	print 'reading ' + str(fname)

	rundir = outpath + str(fname) + '/'
	if not isdir(rundir): raise Exception("I found no such folder: " + rundir)
	pattern = []
	for fnm in glob(rundir + "*.dat"):
		pattern.append(basename(fnm).strip(".dat")) # Get the patterns available
	print "I found these patterns: ", pattern, " in ", rundir

	return rundir, pattern, str(fname)

def get_runs(outpath, rnames=0):
	# Get the means and stdevs of several runs (hpefully many files in a dir)

	print "In " + outpath + ", scanning:"
	if not rnames: rnames = glob(outpath + '*')
	if rnames == []: raise Exception("I found nothing in "+outpath)
	
	patterns = ALLPATS[1:]
	dfs = {}
	tmp = []

	for rname in rnames:	
		run = rname.split('/')[-1]
		runs = {}

		pat = None
		for fnm in glob(rname + "/*.dat"):
			pat = basename(fnm).strip(".dat")

			if pat in patterns:
				runs[pat] = load_pattern(fnm)
				tmp.append(pat)
		if pat is None: continue
		patterns = list(set(patterns) & set(tmp))
		dfs[run] = runs
	print "Only these patterns found in all folders: " + str(patterns)

	return dfs, patterns

def get_baseline(rundir):

	try:
		cf = open(rundir + '/baseline_test.out', "r")
		for line in cf: pass
		# The last line should be like: 
		#	"Initial calibration ical[i], final calibration fcal[i]"
		icalb = float(line.split(", ")[0].split(" ")[2]) 
		fcalb = float(line.split(", ")[1].split(" ")[2])
		cf.close()
	except IOError as e:
		print "No baseline comparison available: IOError:" + str(e)
		exit(1)

	except UnboundLocalError as e:
		print "No baseline comparison available: UnboundLocalError:" + str(e)
		exit(1) 

	return fcalb, icalb

def load_pattern(fname):
	fhead = ['i','d','area','ical','fcal']
	fhead += ['x_1','y_1','x_2','y_2','x_3','y_3']
	fhead += ['no_pointings_x','no_pointins_y']
	fhead += ['frac_0_dith','frac_1_dith','frac_2_dith','frac_3_dith','frac_4_dith']
	return pd.read_table(fname, sep=' ', names=fhead, skiprows=2)

def plot_vs_x(rundir, patterns=["J", "O", "S", "R"], plotssize = [None], linestyles = ['-', '--']):
	# Open all the pattern.dat (except baseline)
	maxx=0; minx=float('inf')
	for p in patterns:
		if p=='baseline': continue

		ts = rundir.strip('/').split('/')[-1]
		fname = op.join(rundir, p+'.dat')
		if not op.isfile(fname): 
			print('Can find the file for ' + p + ' in the ' + ts + ' run. Skip it.')
			patterns.remove(p)
			continue
		tmp = load_pattern(fname)

		# Find limits for x-axis
		if max(tmp.x_1)>maxx: maxx = max(tmp.x_1)
		if min(tmp.x_1)<minx: minx = min(tmp.x_1)

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
			plt.plot(tmp.x_1[test], np.array(tmp.fcal[test]), style, 
					 label = c.convpat(p) + r': '+str(int(npoint))+'x'+str(int(npoint)), c=C[p], lw=2); 
			#plt.plot(tmp.x_1[test], 1 - tmp.frac_1_dith[test] - tmp.frac_2_dith[test] - tmp.frac_3_dith[test] - tmp.frac_4_dith[test], '--', label='0-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_1_dith[test], ':', label='1-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_2_dith[test], ':', label='2-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_3_dith[test], '--', label='3-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_4_dith[test], '--', label='4-pass coverage'); 

	plt.xlabel(r'step size, $d_x$ ["]'); 
	plt.axvline(DX_TRANSITION, color='k', ls=':')
	if 'baseline' in patterns:
		baseline_y, init_y = get_baseline(rundir)
		plt.axvline(BASELINE_x, ymin=0, ymax=(baseline_y-MINY)/(MAXY-MINY), ls='--', c=C['J'], lw=2, zorder=len(patterns))
		plt.scatter(BASELINE_x, baseline_y, marker = "o", s = 50, c=C['J'], label="Laureijs et al. (2011)", zorder=len(patterns))
	plt.ylabel(r'final zero-point scatter, $\sigma_f$'); 
	plt.legend(scatterpoints=1,fontsize=FS,loc=1, ncol=1, handletextpad=0);
	#if max(tmp.d)>(DETX+2*GAPX): 
	#maxx = DETX+2*GAPX
	plt.xlim([minx, maxx])
	plt.ylim([MINY, MAXY])

def plot_vs_d(rundir, patterns=["J", "O", "S", "R"], plotssize = [20, 19], linestyles = ['-', '--']):
	# Open all the pattern.dat (except baseline)
	maxx=0; minx=float('inf')
	for p in patterns:
		if p=='baseline': continue

		ts = rundir.strip('/').split('/')[-1]
		fname = op.join(rundir, p+'.dat')
		if not op.isfile(fname): 
			print('Can find the file for ' + p + ' in the ' + ts + ' run. Skip it.')
			patterns.remove(p)
			continue
		tmp = load_pattern(fname)

		# Find limits for x-axis
		if max(tmp.d)>maxx: maxx = max(tmp.d)
		if min(tmp.d)<minx: minx = min(tmp.d)
		
		# Plot line for each pattern
		for npoint, style in zip(plotssize,linestyles):
			if npoint is not None: 
				test = tmp.no_pointings_x == npoint
				if np.sum(test)==0: 
					print "No runs with "+str(int(npoint))+" pointings found in "+rundir+"!",
					npoint = min(tmp.no_pointings_x-npoint)+npoint
					test = tmp.no_pointings_x == npoint
					print "Plotting npoint="+str(int(npoint))+" instead."
			else:
				test = range(len(tmp))
			plt.plot(tmp.d[test], np.array(tmp.fcal[test]), style, label=c.convpat(p) + r': '+str(int(npoint))+'x'+str(int(npoint)), c=C[p], lw=2); 
		
	plt.xlabel('total distance, $D$ ["]'); 
	#plt.axvline(np.sqrt(5*DX_TRANSITION**2), color='k',ls=':')
	if 'baseline' in patterns:
		baseline_y, init_y = get_baseline(rundir)
		plt.scatter(BASELINE_d, baseline_y, marker = "o", s = 50, c=C['J'], label="Laureijs et al. (2011)", zorder=len(patterns))
		plt.axvline(BASELINE_d, ymin=0, ymax=(baseline_y-MINY)/(MAXY-MINY), ls='--', c=C['J'], lw=2, zorder=len(patterns))
	plt.ylabel(r'final zero-point scatter, $\sigma_f$'); 
	plt.legend(scatterpoints=1, fontsize=FS, loc=1, ncol=1, handletextpad=0);
	if max(tmp.d)>MAXD/4.0: maxx = MAXD/4.0
	plt.xlim([minx, maxx])
	plt.ylim([MINY, MAXY])

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

		# Get the mean of the sigman... eek. It's just for the plot.
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
	