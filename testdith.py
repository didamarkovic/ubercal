"""
	Plotting script for Will's calibration tests.
	Based on his README.bsh bash script.

	16/03/2015 KMarkovic
"""

import os
import numpy as np
import time
from subprocess import call

def get_timestamp():
	t = (2015,3,1,0,0,0,0,0,0)
	return str(int(time.time() - time.mktime(t)))

################### GLOBALS: ################

# Path for execution (either input or find)
thispath = os.path.dirname(os.path.abspath(__file__))

# File where execution time and configuration are saved
final_file = 'run.log'
#config_file = 'config.ini'

# Get timestamp for this run
ts = get_timestamp()	

###############################################

def config_file(pats):
	fnm = 'config'
	for pat in pats: fnm += "-" + pat
	return fnm + '.ini'

def test_dithers(dx=50.0, pattern=PATTERNS, NX=3, NY=None, nsur = None, totcals=None, mode='test', seed=None, dy=None, \
	verb=False, rundir=os.path.join(thispath, 'outputs/', ts), calipath=os.path.join(thispath,'bin'), starfile=os.path.join(thispath,'samples/stars.dat'),cont=False):

	if cont: 
		timestamp = seed
	else:
		timestamp = ts
	
	# If seed has been input, see if there is a run with that seed already
	# 	Should really also check the config is the same!!!: dump pars in run.log
	if not seed or seed=='None': seed = ts
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
	elif mode == 'base':
		totcals = 1
		pattern = ["J"]
		verb = True
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

	# Compile first, before looping 
	#call(['make','clean'], cwd = calipath)
	#call(['make','create-euclid-patch'], cwd = calipath)
	#call(['make','test-calibration'], cwd = calipath)

	# Check if have dy
	if not dy: dy = 2.0*dx
	 
	# At the top of result files will have
	filehead = " with random seed = " + seed + "\n"
	filehead += "i d area ical fcal x_1 y_1 x_2 y_2 x_3 y_3 no_pointings_x no_pointins_y frac_1_dith frac_2_dith frac_3_dit frac_4_dith"
	# Note that the first column is actually the dither factor column
	# It is denoted by the leading # in the header

	# Time the execution
	ti = time.time()

	# First get the baseline and put it into the main run directory (if not testing)
	# 	(This should some day just call a function and below as well!)
	fname = os.path.join(rundir, 'baseline.dat')
	if (mode != 'test' or mode is not 'nobase') and not os.path.isfile(fname):

		if verb: print "Getting the baseline..."

		bl = [50.0,100.0,0.0,100.0,0.0,100.0]
		r = np.sqrt((bl[0]+bl[2]+bl[4])**2 + (bl[1]+bl[3]+bl[5])**2)
		d = np.sqrt(bl[0]**2 + bl[1]**2) + np.sqrt(bl[2]**2 + bl[3]**2) + np.sqrt(bl[4]**2 + bl[5]**2)
				
		# Open baseline log files
		pf = open(os.path.join(rundir, 'baseline_patch.out'), "w")
		cf = open(os.path.join(rundir, 'baseline_test.out'), "w")
			
		# Run Will's code with default dithers
		cmd = ['./create-euclid-patch', rundir+'/', str(bl[0]), str(bl[1]), str(bl[2]), str(bl[3]), str(bl[4]), str(bl[5]), str(NX), str(NY)]
		try:
			call(cmd, stdout=pf, cwd = calipath);
		except OSError as e:
			raise Exception('\n\tThe following command failed:\n\t' + ' '.join(cmd) + '\n\tin folder:\n\t' + calipath + '\n\twith the error:\n\t' + str(e))
		cmd = ['./test-calibration', rundir+'/', seed, starfile]
		try:
			call(cmd, stdout=cf, cwd = calipath);
		except OSError as e:
			raise Exception('\n\tThe following command failed:\n\t' + ' '.join(cmd) + '\n\tin folder:\n\t' + calipath + '\n\twith the error:\n\t' + str(e))
		
		# Now get the area in patch info
		line = None
		with open(os.path.join(rundir, 'baseline_patch.out'), "r") as pf:
			frac=[]
			for line in pf: 
				try:
					if "exposures over area" in line:
						frac += [float(line.split(", fraction ")[-1])]
					elif "tot" in line and "area" in line and "=" in line:
						area = float(line.split("= ")[-1])
					else:
						pass
				except ValueError:
					raise Exception("This line makes no sense:\n" + line)
		if line==None: raise Exception("File " + rundir + 'baseline_patch.out' + " is empty!") 

		# Now get last line of the file
		pf.close(); cf.close()
		cf = open(os.path.join(rundir,'baseline_test.out'), "r")
		for line in cf: 
			if verb: pass
		# The last line should be like: 
		#	"Initial calibration ical[i], final calibration fcal[i]"
		if line==None or "calibration" not in line: raise Exception("File " + rundir + 'baseline_test.out' + " is not what I expected!")
		icalb = float(line.split(", ")[0].split(" ")[2]) 
		fcalb = float(line.split(", ")[1].split(" ")[2])
		cf.close()

		row = [0, d, area] + [icalb,fcalb] + list(bl) + [NX, NY] + frac
		np.savetxt(fname, np.array(row)[None], fmt='%.6f', delimiter=' ', newline='\n', header="J-pattern: baseline" + filehead)

		if verb: print "Baseline saved to " + fname

	if mode != 'base':
		# Then loop over patterns
		for p in range(len(pattern)):

			if verb: print "Doing the " + pattern[p] + "-pattern now."

			# Save the directory for this pattern in the run directory
			patdir = os.path.join(rundir,pattern[p])
			if mode is 'test':
				patdir += '-'
			else:	
				if not os.path.isdir(patdir): 
					call(['mkdir',patdir])
					if verb: print "Created " + os.path.dirname(patdir)

			# Build up an array of dither strategy multiples
			# 	it is in arcsec (n.b. one CCD is 2040 pixel = 612", gap about 50", 100")
			if pattern[p]=="J":
				end_dithers = np.array([[dx,dy],[0.0,dy],[0.0,dy]])
			elif pattern[p]=="box":
				end_dithers = np.array([[dx,0.0],[0.0,dy],[-dx,0.0]])
			elif pattern[p]=="step":
				end_dithers = np.array([[dx,dy],[dx,0.0],[dx,0.0]])
			elif pattern[p]=="S":
				end_dithers = np.array([[dx,dy],[0.0,dy],[dx,dy]])
			elif pattern[p]=="X":
				end_dithers = np.array([[dx,dy],[0.0,0.0],[-dx,-dy]])
			elif pattern[p]=="N":
				end_dithers = np.array([[dx,dy],[dx,0.0],[dx,dy]])
			else:
				raise RuntimeError("You did not tell me what the " + pattern[p]+ "-pattern is!\n")
			
			# Now loop over multiples of dither strategy: 
			# 	resulting calibration gets saved into a file each time => read into array
			dither_arr = []
			for i in range(totcals):

				# Calculate the next dither vector as fraction of final, largest pattern
				new_dithers = np.array(end_dithers)
				if totcals != 1:
					new_dithers = np.array(end_dithers)*float(i)/(totcals-1)
				x1 = new_dithers[0,0]
				y1 = new_dithers[0,1]
				x2 = new_dithers[1,0]
				y2 = new_dithers[1,1]
				x3 = new_dithers[2,0]
				y3 = new_dithers[2,1]
				# Calculate the vector magnitude of 4th dither displacement from 1st
				r = np.sqrt((x1+x2+x3)**2 + (y1+y2+y3)**2)
				# Calculate full distance travelled by telescope (in pix)
				d = np.sqrt(x1**2 + y1**2) + np.sqrt(x2**2 + y2**2) + np.sqrt(x3**2 + y3**2)
						
				# Create directory for this dither size
				outdir = os.path.join(patdir, 'dithersize-'+str(x1))
				if not os.path.isdir(outdir): call(['mkdir', outdir+'/'])
				
				# Create mangle result files
				skippatch = False
				patchfile = os.path.join(outdir, 'patch-'+str(x1)+'.out')
				#if os.path.isfile(patchfile):
				#	pf = open(patchfile, "r")
				#	for line in pf: pass
				#	pf.close()
				#	if "total area =" in line: 
				#		skippatch = True
				#		if verb: "The patch was already done, results in " + patchfile + "."
				# The way the code works currently (v0.2), we should never skip making the patch as it is done for the full survey each time!							
				
				for nx, ny in surveys:

					if not skippatch:
						# Run Will's code by calling bash
						pf = open(patchfile, "w")
						call(['./create-euclid-patch', outdir+'/', str(x1), str(y1), str(x2), str(y2), str(x3), str(y3), str(nx), str(ny)], stdout=pf, cwd = calipath)
						pf.close();

					# Now get the area in patch info
					line = None
					with open(patchfile, "r") as pf:
						frac=[]
						line = None
						for line in pf: 
							if "exposures over area" in line:
								frac += [float(line.split(", fraction ")[-1])]
							elif "tot" in line and "area" in line and "=" in line:
								area = float(line.split("= ")[-1])
							else:
								pass
					if line==None: raise Exception("File " + patchfile + " is empty!")

					if verb: print "Calibration for a survey with " + str(nx) + "x" + str(ny) + " pointings."

					outdira = os.path.join(outdir, str(nx)+'x'+str(ny))
					califile = os.path.join(outdir, 'test-'+str(x1)+'.out')		
					if not os.path.isfile(califile):

						if verb: print "No " + califile

						# Create mangle result files
						cf = open(califile, "w")
					
						# Run Will's code by calling bash
						call(['./test-calibration', outdir+'/', seed, starfile], stdout=cf, cwd = calipath)
						cf.close();
						
						# Now get last line of the file (to get initial and final claibrations)
						cf = open(califile, "r")
						for line in cf: pass						
						if line==None or "calibration" not in line: raise Exception("File " + califile + " is not what I expected! I got this line:\n" + line)
					else:
						print "Reading previous result from " + califile
						cf = open(califile, "r")
						for line in cf: pass

						if not "final calibration" in line:

							print califile + " not finished."

							# Create mangle result files
							cf = open(califile, "w")
						
							# Run Will's code by calling bash
							call(['./test-calibration', outdir+'/', seed, starfile], stdout=cf, cwd = calipath)
							cf.close();
							
							cf = open(califile, "r")
							for line in cf: pass
							if line==None or "calibration" not in line: raise Exception("File " + califile + " is not what I expected!")

					# The last line should be like: 
					#	"Initial calibration ical[i], final calibration fcal[i]"
					ical = float(line.split(", ")[0].split(" ")[2]) 
					fcal = float(line.split(", ")[1].split(" ")[2])
					cf.close()

					# Append results of this dither size the result array
					row = [int(i), d, area] + [ical,fcal] + list(np.hstack(new_dithers)) + [nx, ny] + frac
					dither_arr.append(row)
						# This will likely be the biggest RAM consumer (but still shouldn't go beyond ~100MB)
					
					if mode == 'area':
						if not os.path.isdir(outdira): call(['mkdir', outdira])
						if verb: print "Stashing results to: "+outdira	
						call(['mv', califile, outdira])
						call(['mv', patchfile, outdira])
						call(['mv', os.path.join(outdir,"full-survey-overlaps.txt"), outdira])
						call(['mv', os.path.join(outdir,"full-pointing.pol"), outdira])
						call(['mv', os.path.join(outdir,"full-pointing.vrt"), outdira])

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
	np.savetxt(os.path.join(rundir, final_file), outtext, fmt='# %s')
	
	print outtext[0]

	return timestamp, os.path.join(rundir,final_file)


######################################################### MAIN ############################################
if __name__=='__main__':
	import argparse

	parser = argparse.ArgumentParser(description="Calculate self-calibration quality for different dither patterns.")
	
	# General options
	parser.add_argument("-v", "--verbose", default=False, action="store_true")
	parser.add_argument("-f", "--carryon", default=False, action="store_true", help="continue a previous run")
	
	# Mutually exclusive mode arguments
	mode_group = parser.add_mutually_exclusive_group()
	mode_group.add_argument("-m", "--mode", type=str, default='test', help="mode of execution: leave blank for production, 'test' for a simple test, 'base' for only baseline J-pattern, 'nobase' to skip calculating for the baseline, 'area':study variation with survey size (number of pointings) ")
	mode_group.add_argument("-t", "--test", action='store_true', default=False, help="mode = 'test'")
	
	# These should be a non-mutually exclusive group:
	parser.add_argument("-x", "--xmax", type=float, default=50.0, help="x-displatement of first dither")
	parser.add_argument("-y", "--ymax", type=float, default=None, help="y-displatement of first dither")
	parser.add_argument("-p", "--patterns", nargs='+', default=PATTERNS, help="dithering patterns to test: implemented: J, step, S, box")
	parser.add_argument("-s", "--seed", type=int, default='-'+ts, help="random seed for stellar density")
	parser.add_argument("-nx", type=int, default=3, help="number of pointings in survey in x-direction (RA)")
	parser.add_argument("-ny", type=int, default=None, help="number of pointings in survey in y-direction (dec)")
	
	# As below, these two should be mutually exclusive for safety. Setting both of them runs, but hasn't been tested.
	mode_group_2 = parser.add_mutually_exclusive_group()
	mode_group_2.add_argument("-n", "--nsizes", type=int, default=1, help="number of pattern sizes to run (max set by xmax argument)")
	mode_group_2.add_argument("-a", "--surveys", type=int, default=None, help="how many survey sizes to test")

	parser.add_argument("-o", "--outpath", default=thispath, help="where you want or have your 'outputs' folder")
	parser.add_argument("-c", "--calipath", default=thispath+"/bin/", help="directory containing test-calibration binaries")
	parser.add_argument("-st", "--starpath", default=thispath+'/samples/stars.dat', help="file containing a table of stellar populations, densities and their SNR")

	# Save the inputs:
	args = parser.parse_args()

	# Set the missing redundant settings:
	if not args.surveys == None: 
		args.mode = 'area'
	elif (args.nsizes or args.xmax or args.ymax or args.patterns or args.seed) and args.mode=='test':
		args.mode = 'producion'	
	
	# As above, these two should be mutually exclusive for safety. Setting both of them runs, but hasn't been tested.
	if args.nsizes > 1 and args.mode=='area': raise Exception("argument -n/--nsizes not allowed in area mode")

	# Make directory structure for outputs for this run
	rundir = os.path.join(args.outpath, 'outputs/', str(abs(args.seed)))
	if not os.path.isdir(os.path.join(args.outpath,'outputs')):
		call(['mkdir', os.path.dirname(rundir)])
		call(['mkdir', rundir])
		print "Created " + os.path.dirname(rundir)
	# If the given seed has no directory, create one from the timestamp
	elif not os.path.isdir(rundir):
		call(['mkdir', rundir])
		if args.verbose: print "Created " + os.path.dirname(rundir)
			
	# Save the configuration explicitly (no INTERRUPTED at end):
	config_file = config_file(args.patterns)
	f = open(os.path.join(rundir, config_file),'w')
	if args.carryon: 
		print "... CONTINUING INTERRUPTED RUN " + str(args.seed) + "."
		f.write("# This run was interrupted and restarted.\n")
	f.write("# Run " + ts + " has configuration:\n")
	for arg, value in sorted(vars(args).items()):
		f.write(str(arg) + '=' + str(value) + '\n')
	
	# Now run:
	try:
		[ts, fn] = test_dithers(dx=args.xmax,pattern=args.patterns, NX = args.nx, NY = args.ny, 
			nsur = args.surveys, totcals=args.nsizes, \
			mode=args.mode, seed=str(args.seed), dy=args.ymax, 
			verb=args.verbose, rundir=rundir, calipath=args.calipath, starfile=args.starpath, cont=args.carryon)
	except Exception as e:
		f.write("...INTERRUPTED")
		print "RUN WAS INTERRUPTED... (saved to " + config_file + ")."
		raise
	finally:
		f.close()

	print "Run " + ts + " is finished, configuration saved to " + config_file + "."
