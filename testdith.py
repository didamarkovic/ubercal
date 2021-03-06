"""
	Plotting script for Will's calibration tests.
	Based on his README.bsh bash script.

	16/03/2015 KMarkovic
"""

import os
import numpy as np
import time
from subprocess import call
import traceback
import ubercal

def get_timestamp():
	t = (2016,1,1,0,0,0,0,0,0)
	return str(int(time.time() - time.mktime(t)))

################### GLOBALS: ################

# Path for execution (either input or find)
thispath = os.path.dirname(os.path.abspath(__file__))

# File where execution time and configuration are saved
final_file = 'run.log'
#config_file = 'config.ini'

# Get timestamp for this run
TS = get_timestamp()

# Default patterns to test
PATTERNS = ["J", "S"]	

###############################################

def config_filename(pats):
	fnm = 'config'
	for pat in pats: fnm += "-" + pat
	return fnm + '.ini'

def test_dithers(dx=50.0, pattern=PATTERNS, NX=3, NY=None, nsur = None, totcals=None, mode='test', \
				 seed=None, dy=None, verb=False, rundir=os.path.join(thispath, 'outputs/', TS), \
				 calipath=os.path.join(thispath,'bin'), starfile=os.path.join(thispath,'samples/stars.dat'),
				 cont=False, detbool=ubercal.DETBOOL, ftol=ubercal.FTOL):

	if cont: 
		timestamp = seed
	else:
		timestamp = TS
	
	# If seed has been input, see if there is a run with that seed already
	# 	Should really also check the config is the same!!!: dump pars in run.log
	if not seed or seed=='None': seed = TS
	if int(seed) > 0: seed = "-"+seed
	if verb: print "Running with seed = " + seed + "."
	if not NY: NY = NX

	surveys = [[NX, NY]]
	arrsteps_x = [0]

	# Special settings for quick TESTING MODE:
	if verb: print "Running in mode '" + mode + "'."
	if mode == 'test':
		if not totcals: totcals = 1
		pattern = ["J"]
		verb = True
		seed = '-1'
		if verb: print "Seeting seed to -1, since we re in testing mode."
	elif mode == 'base':
		totcals = 1
		pattern = ["J"]
		verb = True
		if verb: print "Seeting seed to -1, since we re in baseline-only mode."
		seed = '-1'
	elif mode == 'area':
		if not nsur: nsur = 2*NX
		if nsur > NX/2: 
			arrsteps_x = np.arange(nsur+1)[1:]
		else:
			arrsteps_x = np.arange(nsur) + NX - int(nsur/2)
		if nsur > NY/2: 
			arrsteps_y = np.arange(nsur+1)[1:]
		else:
			arrsteps_y = np.arange(nsur)  + NY - int(nsur/2)
		surveys = zip(arrsteps_x,arrsteps_y) 

	# Check if have dy
	if not dy: dy = 2.0*dx
	 
	# At the top of result files will have
	filehead = " with random seed = " + seed + "\n"
	filehead += "i d area ical fcal x_1 y_1 x_2 y_2 x_3 y_3 no_pointings_x no_pointins_y frac_1_dith frac_2_dith frac_3_dith frac_4_dith"
	# Note that the first column is actually the dither factor column
	# It is denoted by the leading # in the header

	# Time the execution
	ti = time.time()

	# First get the baseline and put it into the main run directory (if not testing)
	# 	(This should some day just call a function and below as well!)
	fname = os.path.join(rundir, 'baseline.dat')
	if (mode != 'test' or mode is not 'nobase') and not os.path.isfile(fname):

		if verb: print "\n--------- Getting the baseline... --------"

		bl = ubercal.dith.BASELINE
		r, d = ubercal.dith.totaldither(bl)

		### Run Will's code to create the baseline survey with default dithers
		area, frac = ubercal.dith.create_surveypatch(bl, NX, NY, calipath, rundir, os.path.join(rundir, 'baseline_patch.out'), verb)

		### Run Will's code to test the ubercalibration using the baseline survey
		icalb, fcalb = ubercal.test_calibration(starfile, seed, calipath, rundir, os.path.join(rundir, 'baseline_test.out'), verb, detbool, ftol)

		# Save the above into the baseline output file
		row = [0, d, area] + [icalb,fcalb] + list(bl) + [NX, NY] + frac
		np.savetxt(fname, np.array(row)[None], fmt='%.6f', delimiter=' ', newline='\n', header="J-pattern: baseline" + filehead)

		if verb: print "Baseline saved to " + fname + '.'

	if mode != 'base':
		# Then loop over patterns
		for p in range(len(pattern)):

			if verb: print "\n--------- Doing the " + pattern[p] + "-pattern now. --------"

			# Save the directory for this pattern in the run directory
			patdir = os.path.join(rundir,pattern[p])
			if mode is 'test':
				patdir += '-'
			else:	
				if not os.path.isdir(patdir): 
					call(['mkdir',patdir])
					if verb: print "Created " + patdir

			# Build up an array of dither strategy multiples
			end_dithers = ubercal.dith.getpat(pattern[p], dx, dy)
			
			# Now loop over multiples of dither strategy: 
			# 	resulting calibration gets saved into a file each time => read into array
			dither_arr = []
			for i in range(totcals):

				# Calculate the next dither vector as fraction of final, largest pattern
				new_dithers = np.array(end_dithers)
				if totcals != 1:
					new_dithers = np.array(end_dithers)*float(i)/(totcals-1)

				r, d = ubercal.dith.totaldither(new_dithers.flatten())

				# Create directory for this dither size
				outdir = os.path.join(patdir, 'dithersize-'+str(new_dithers[0,0]))
				if not os.path.isdir(outdir): 
					call(['mkdir', outdir+'/'])
					if verb: print "Created " + outdir
				
				# Create mangle result files
				skippatch = False
				patchfile = os.path.join(outdir, 'patch-'+str(new_dithers[0,0])+'.out')
				#if os.path.isfile(patchfile):
				#	pf = open(patchfile, "r")
				#	line = None
				#	for line in pf: pass
				#	pf.close()
				#	if "total area =" in line: 
				#		skippatch = True
				#		if verb: "The patch was already done, results in " + patchfile + "."
				# The way the code works currently (v0.2), we should never skip making the patch as it is done for the full survey each time!							
				
				for nx, ny in surveys:

					if mode == 'area':
						outdira = os.path.join(outdir, str(nx)+'x'+str(ny))
					else:
						outdira = outdir

					### Construct this survey

					# The following is another awful hack:
					if mode!='area' and os.path.isfile(patchfile): 
						skippatch = True
						line = None
						with open(patchfile, "r") as pf:
							for line in pf: pass
						if line is None: skippatch = False
					else:
						skippatch = False

					# Read or create survey geometry files
					patchfile = os.path.join(outdira, 'patch-'+str(new_dithers[0,0])+'.out')
					if skippatch:
						try:
							area, frac = ubercal.dith.get_area(patchfile, verb)
						except:
							skippatch = False
					if not skippatch:				
						# Run Will's code to create the survey files
						area, frac = ubercal.dith.create_surveypatch(new_dithers.flatten(), 
																	 nx, ny, calipath, outdir, 
																	 patchfile, verb)						

					### Test calibration

					if verb: print "Calibration for a survey with " + str(nx) + "x" + str(ny) + " pointings."

					califile = os.path.join(outdir, 'test-'+str(new_dithers[0,0])+'.out')		
					
					if not os.path.isfile(califile) and not os.path.isfile(os.path.join(outdira, 'test-'+str(new_dithers[0,0])+'.out')):

						if verb: print "No " + califile

						ical, fcal = ubercal.test_calibration(starfile, seed, calipath, outdir, califile, verb, detbool, ftol)
										
					else:
						
						if os.path.isfile(os.path.join(outdira, 'test-'+str(new_dithers[0,0])+'.out')):
							califile = os.path.join(outdira, 'test-'+str(new_dithers[0,0])+'.out')

						print "Reading previous result from " + califile

						ical, fcal = ubercal.get_cals(califile, True, verb)

						# If the file is unfinished, run test-calibration again
						if ical is None or fcal is None:
							ical, fcal = ubercal.test_calibration(starfile, seed, calipath, outdir, califile, verb, detbool, ftol)

					# Append results of this dither size the result array
					row = [int(i), d, area] + [ical,fcal] + list(np.hstack(new_dithers)) + [nx, ny] + frac
					dither_arr.append(row)
					# This will likely be the biggest RAM consumer (but still shouldn't go beyond ~100MB)
					
					if mode == 'area':
						stashprint = False
						if not os.path.isdir(outdira): call(['mkdir', outdira])
						if outdira not in patchfile:
							stashprint = True
							call(['mv', patchfile, outdira])
							call(['mv', os.path.join(outdir,"full-survey-overlaps.txt"), outdira])
							call(['mv', os.path.join(outdir,"full-survey.vrt"), outdira])
							call(['mv', os.path.join(outdir,"full-survey-overlaps.vrt"), outdira])
						if outdira not in califile:	
							stashprint = True
							call(['mv', califile, outdira])
						if verb and stashprint: print "Stashed results to: "+ outdira + '.\n'							
						# If we go back to creating the patch first, the patch files stay the same for different survey sizes, 
						#  so they should not be stashed, but reused for speed!

			# Save all calibrations from each pattern into a file (or append if it exists)
			dither_arr = np.array(dither_arr)
			fname = os.path.join(rundir, pattern[p]+'.dat')
			if verb: print "Saving result array to " + fname + "."
			if os.path.isfile(fname):
				#fa = file(fname, 'a')
				call(['mv', fname, fname+'-bu'])
				if verb: print "Backed up " + fname + "."
			np.savetxt(fname, dither_arr, fmt='%.6f', delimiter=' ', newline='\n', header=pattern[p]+"-pattern" + filehead)#, footer='', comments='# ')

			# Print results to screen
			if verb: print pattern[p] + " pattern"
			if verb: print filehead
			if verb: print row

	# Save a little log file at end (expand info in future)
	duration = str(time.time() - ti)
	nodith = str(totcals*len(pattern))
	outtext = ['Run took ' + duration + ' seconds for ' + nodith + ' different dither configs, for '+ str(len(arrsteps_x))+' survey(s).']
	
	if mode is 'test': outtext.extend(['This was a run in test mode.'])
	if mode is 'base': outtext = ['Run took ' + duration + ' seconds for baseline dither config.']

	outtext.extend(str(ubercal.io.GitEnv()).split('\n'))

	np.savetxt(os.path.join(rundir, final_file), outtext, fmt='# %s')
	
	print outtext[0]

	return timestamp, os.path.join(rundir,final_file)


def run(args, patterns, rundir, seed='-'+TS, carryon=False, verbose=False):

	# Name the directory for this run (i.e. this seed) and create it if it doesn't exist
	if not os.path.isdir(rundir):
		os.makedirs(rundir)
		if verbose: print " and in it " + os.path.basename(rundir) + "."
	else:
		if verbose: print "."

	# Save the configuration
	config_file = os.path.join(rundir, config_filename(patterns))
	head = ""
	mode = "w"
	if carryon and os.path.isfile(config_file): 
		mode = "a"
		print "... CONTINUING INTERRUPTED RUN " + str(seed) + "."
		head += "\n\n# ... CONTINUE THE RUN ...\n"
	elif os.path.isfile(config_file):
		raise Exception("The run already exists in\n\t"+rundir+\
					    "\nIf you'd like to continue it, use the -f flag.")
	elif carryon:
		print "You say I should CONTINUE AN INTERRUPTED RUN, but there is no evidence of any " +\
			  "previous results from seed " + str(seed) + "!"
		print "\t => starting from scratch."

	with open(config_file,mode) as f:
	
		# Write run metadata into the opened configuration file
		f.write(head+"# Run " + seed[1:] + " has configuration:\n")
		for arg, value in sorted(vars(args).items()):
			f.write(str(arg) + '=' + str(value) + '\n')
	
		# Now run:
		exitcode=0
		try:
			[ts, fn] = test_dithers(dx=args.xmax,
									pattern=patterns, 
									NX = args.nx, 
									NY = args.ny, 
									nsur = args.nsurveys, 
									totcals=args.nsizes, 
									mode=args.mode, 
									seed=str(seed), 
									dy=args.ymax, 
									verb=args.verbose, 
									rundir=rundir, 
									calipath=os.path.abspath(args.calipath), 
									starfile=os.path.abspath(args.starpath), 
									cont=args.carryon, 
									detbool=args.detbool, 
									ftol=args.ftol)
		except:
			tb = traceback.format_exc()
			f.write("\n#...INTERRUPTED due to:\n" + tb)
			print "RUN WAS INTERRUPTED... Config and traceback saved to " + \
				  os.path.join(rundir, config_file) + ".\n\t" + tb.splitlines()[-1]
			exitcode=1
		else:
			print "Run " + seed[1:] + " is finished, configuration saved to " + config_file +\
				  ", log in " + os.path.basename(fn) + "."

	return exitcode

######################################################### MAIN ############################################
if __name__=='__main__':
	import argparse

	parser = argparse.ArgumentParser(description="Calculate self-calibration quality for different dither patterns.")
	
	# General options
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-m", "--mode", type=str, default='test', help="mode of execution: leave blank for production, 'test' for a simple test, 'base' for only baseline J-pattern, 'nobase' to skip calculating for the baseline, 'area':study variation with survey size (number of pointings) ")

	# Mutually exclusive mode arguments
	group.add_argument("-t", "--test", action='store_true', default=False, help="mode = 'test'")
	parser.add_argument("-v", "--verbose", default=False, action="store_true")
	parser.add_argument("-f", "--carryon", default=False, action="store_true", help="continue a previous run")
	
	# These should be a non-mutually exclusive group:
	parser.add_argument("-x", "--xmax", type=float, default=50.0, help="x-displatement of first dither")
	parser.add_argument("-y", "--ymax", type=float, default=None, help="y-displatement of first dither")
	parser.add_argument("-p", "--patterns", nargs='+', default=PATTERNS, help="dithering patterns to test: implemented: J, Jsq, R, S, Ssq, N, X, O")
	parser.add_argument("-s", "--seed", type=int, default='-'+TS, help="random seed for stellar density (ignored if -b > 1")
	parser.add_argument("-nx", type=int, default=3, help="number of pointins in survey in x-direction (RA)")
	parser.add_argument("-ny", type=int, default=None, help="number of pointins in survey in y-direction (dec)")

	# These in principle work together, so don't have to be mutually exclusive
	parser.add_argument("-a", "--nsurveys", type=int, default=1, help="how many survey sizes to test")
	parser.add_argument("-b", "--nseeds", type=int, default=1, help="how many seed to run for better stats")
	parser.add_argument("-n", "--nsizes", type=int, default=1, help="number of pattern sizes to run (max set by xmax argument)")

	parser.add_argument("-o", "--outpath", default=thispath, help="where you want or have your 'outputs' folder")
	parser.add_argument("-c", "--calipath", default=thispath+"/bin/", help="directory containing test-calibration binaries")
	parser.add_argument("-st", "--starpath", default=thispath+'/samples/stars.dat', help="file containing a table of stellar populations, densities and their SNR")

	# Some additional options	
	parser.add_argument("-d", "--detbool", action='store_true', default=False, help="vary 0-points detector to detector instead of exposure to exposure")
	parser.add_argument("-e", "--ftol", default=ubercal.FTOL, help="tolerance parameter for the optimisation")
	parser.add_argument("-l", "--compile", default=None, help="compile the C-code before running - input the source code folder")

	# Save the inputs:
	args = parser.parse_args()
	
	# Set the missing redundant settings:
	if args.nsurveys>1:
		args.mode = 'area'
	elif (args.nsizes or args.xmax or args.ymax or args.patterns or args.seed) and args.mode=='test':
		args.mode = 'producion'	

	# Compile first, before looping 
	if args.compile: ubercal.compile(os.path.abspath(args.compile))

	# Generate seed list
	if args.nseeds == 1: seeds = [args.seed]
	else: seeds = -1*np.arange(args.nseeds).astype(int) - int(TS)
	
	# Make directory structure for outputs for this run

	# Create the name for the outputs directory and create it if it does not exist
	if 'outputs' != os.path.basename(os.path.normpath(args.outpath)):
		outdir = os.path.abspath(os.path.join(args.outpath, 'outputs/'))
	else:
		outdir = os.path.abspath(args.outpath)

	# Add a DET or EXP folder depending on what is requested
	if args.detbool:
		outdir = os.path.join(outdir, 'DET/')
	else:
		outdir = os.path.join(outdir, 'EXP/')

	# Make everything if it doesn't already exist
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		if args.verbose: print "Created " + outdir, 

	# Run the test once for each seed
	for s in seeds:
		rundir = os.path.join(outdir, str(abs(s)))
		exitcode = run(args, args.patterns, rundir, str(s), args.carryon, args.verbose)
	
	# Exit
	exit(exitcode)
