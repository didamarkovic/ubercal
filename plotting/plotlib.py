""" 
	For plotting the dither tests.
	
	17/03/2015, KMarkovic
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from os.path import basename, isfile, isdir
from glob import glob
import pandas as pd
import coverage as c

FS = 14
matplotlib.rcParams.update({'font.size': FS})

MINY = 1
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

C = {
	'J': 'k',
	'S': 'g',
	'X': 'r',
	'R': 'c',
	'step': 'c',
	'O': 'm',
	'box': 'm',
	'N': 'y'
	}

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
	
	patterns = ["J", "S", "step", "box"]
	dfs = {}
	tmp = []
	
	for rname in rnames:	
		run = rname.split('/')[-1]
		runs = {}

		for fnm in glob(rname + "/*.dat"):
			pat = basename(fnm).strip(".dat")

			if pat in patterns:
				runs[pat] = load_pattern(fnm)
				tmp.append(pat)

		patterns = list(set(patterns) & set(tmp))
		dfs[run] = runs
	print "Only these patterns found in all folders: " + str(patterns)

	return dfs, patterns

def get_baseline(rundir):

	try:
		cf = open(rundir + 'baseline_test.out', "r")
		for line in cf: pass
		# The last line should be like: 
		#	"Initial calibration ical[i], final calibration fcal[i]"
		icalb = float(line.split(", ")[0].split(" ")[2]) 
		fcalb = float(line.split(", ")[1].split(" ")[2])
		cf.close()
	except IOError:
		print "No baseline comparison available."
		exit(1)
		#pattern = ["J", "box"]
		#fcalb = 0.0123985 
	except UnboundLocalError:
		print "No baseline comparison available."
		exit(1)
		#pattern = ["J", "box"]
		#fcalb = 0.0123985 

	return fcalb/icalb

def load_pattern(fname):
	fhead = ['i','d','area','ical','fcal']
	fhead += ['x_1','y_1','x_2','y_2','x_3','y_3']
	fhead += ['no_pointings_x','no_pointins_y']
	fhead += ['frac_0_dith','frac_1_dith','frac_2_dith','frac_3_dith','frac_4_dith']
	return pd.read_table(fname, sep=' ', names=fhead, skiprows=2)

def plot_vs_x(rundir, patterns=["J", "box", "S", "step"], plotssize = [None], linestyles = ['-', '--']):
	# Open all the pattern.dat (except baseline)
	maxx=0; minx=float('inf')
	for p in patterns:
		if p=='baseline': continue

		fname = rundir + p + '.dat'
		tmp = load_pattern(fname)

		# Find limits for x-axis
		if max(tmp.x_1)>maxx: maxx = max(tmp.x_1)
		if min(tmp.x_1)<minx: minx = min(tmp.x_1)

		# Plot line for each pattern
		for npoint, style in zip(plotssize,linestyles):
			if npoint is not None: 
				test = tmp.no_pointings_x==npoint
				if np.sum(test)==0: raise Exception("No runs with "+str(npoint)+" pointings found in "+rundir+"!")
			else:
				test = range(len(tmp))
			plt.plot(tmp.x_1[test], np.array(tmp.ical[test])/np.array(tmp.fcal[test]), style, label=c.convpat(p) + r': '+str(npoint)+'x'+str(npoint), c=C[p], lw=2); 
			#plt.plot(tmp.x_1[test], 1 - tmp.frac_1_dith[test] - tmp.frac_2_dith[test] - tmp.frac_3_dith[test] - tmp.frac_4_dith[test], '--', label='0-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_1_dith[test], ':', label='1-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_2_dith[test], ':', label='2-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_3_dith[test], '--', label='3-pass coverage'); 
			#plt.plot(tmp.x_1[test], tmp.frac_4_dith[test], '--', label='4-pass coverage'); 

	plt.xlabel(r'step size, $d_x$ ["]'); 
	plt.axvline(DX_TRANSITION, color='k', ls=':')
	if 'baseline' in patterns:
		baseline_y = 1.0/get_baseline(rundir)
		plt.axvline(BASELINE_x, ymin=0, ymax=(baseline_y-MINY)/(MAXY-MINY), ls='--', c=C['J'], lw=2, zorder=len(patterns))
		plt.scatter(BASELINE_x, baseline_y, marker = "o", s = 50, c=C['J'], label="Laureijs et al. (2011)", zorder=len(patterns))
	plt.ylabel(r'improvement ratio, $q$'); 
	plt.legend(scatterpoints=1,fontsize=FS,loc=1,ncol=3, handletextpad=0);
	#if max(tmp.d)>(DETX+2*GAPX): 
	maxx = DETX+2*GAPX
	plt.xlim([minx, maxx])
	plt.ylim([MINY, MAXY])

def plot_vs_d(rundir, patterns=["J", "box", "S", "step"], plotssize = [20, 19], linestyles = ['-', '--']):
	# Open all the pattern.dat (except baseline)
	maxx=0; minx=float('inf')
	for p in patterns:
		if p=='baseline': continue

		fname = rundir + p + '.dat'
		tmp = load_pattern(fname)

		# Find limits for x-axis
		if max(tmp.d)>maxx: maxx = max(tmp.d)
		if min(tmp.d)<minx: minx = min(tmp.d)
		
		# Plot line for each pattern
		for npoint, style in zip(plotssize,linestyles):
			if npoint is not None: 
				test = tmp.no_pointings_x==npoint
				if np.sum(test)==0: raise Exception("No runs with "+str(npoint)+" pointings found in "+rundir+"!")
			else:
				test = range(len(tmp))
			plt.plot(tmp.d[test], np.array(tmp.ical[test])/np.array(tmp.fcal[test]), style, label=c.convpat(p) + r': '+str(npoint)+'x'+str(npoint), c=C[p], lw=2); 
		
	plt.xlabel('total distance, $D$ ["]'); 
	#plt.axvline(np.sqrt(5*DX_TRANSITION**2), color='k',ls=':')
	if 'baseline' in patterns:
		baseline_y = 1.0/get_baseline(rundir)
		print BASELINE_d
		plt.scatter(BASELINE_d, baseline_y, marker = "o", s = 50, c=C['J'], label="Laureijs et al. (2011)", zorder=len(patterns))
		plt.axvline(BASELINE_d, ymin=0, ymax=(baseline_y-MINY)/(MAXY-MINY), ls='--', c=C['J'], lw=2, zorder=len(patterns))
	plt.ylabel(r'improvement ratio, $q$'); 
	plt.legend(scatterpoints=1, fontsize=FS, loc=2, ncol=3, handletextpad=0);
	if max(tmp.d)>MAXD/4.0: maxx = MAXD/4.0
	plt.xlim([minx, maxx])
	plt.ylim([MINY, MAXY])

def plot_vs_nx(rundir, patterns=["J", "box", "S", "step"], dx=[0.0, 50.0], linestyles = ['-', '--'], ms=10):
	# Open all the pattern.dat (except baseline)
	maxx=0; minn=0
	for p in patterns:
		if p=='baseline': continue

		fname = rundir + p + '.dat'
		tmp = load_pattern(fname)
	
		# Find limits for x-axis
		if max(tmp.no_pointings_x)>maxx: maxx = max(tmp.no_pointings_x)**2
		if min(tmp.no_pointings_x)<minn: minn = min(tmp.no_pointings_x)**2
		
		# Plot line for each pattern
		for dsize, style in zip(dx,linestyles):
			test = tmp.x_1==dsize
			plt.plot(np.array(tmp.no_pointings_x[test]), np.array(tmp.fcal[test])/np.array(tmp.ical[test]), style, ms=ms)#, label=p + '-pattern: '+str(dsize)+'"'); 
		
	plt.xlabel('root number of pointings in survey'); 
	plt.ylabel(r'improvement ratio, $q$'); 
	plt.legend(bbox_to_anchor=(1.0, 0.3), fontsize=FS);
	#plt.xlim([minn, maxx])

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

	# Assuming patters smae in tables (get_tables gives it):
	patterns = tables[tables.keys()[0]].keys()

	# Assuming all are same length, x-axis!
	test = tables[tables.keys()[0]][patterns[0]].x_1==dx
	print tables[tables.keys()[0]][patterns[0]].x_1
	print dx
	nx_vec = tables[tables.keys()[0]][patterns[0]].no_pointings_x[test]  
	stdevs = {}
	means = {}
	td_vec={}

	for pat in patterns:

		print pat
		mean_r = np.zeros(len(nx_vec))
		meansq_r = mean_r*0.0 # Was asigning not copying!!!!!
		sd_V = mean_r*0.0
		tmp = sd_V*0.0
	
		for key in tables.keys():

			table = tables[key][pat]
			ratio = np.array(table.fcal[test])/np.array(table.ical[test])

			mean_r += ratio/len(tables)
			meansq_r += ratio**2/len(tables)

		stdevs[pat] = np.sqrt(meansq_r - mean_r**2)
		means[pat] = mean_r
		td_vec[pat] = table.d[test][0] 

	return nx_vec, means, stdevs, td_vec
