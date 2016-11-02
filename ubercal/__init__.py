""" Python functions to run Will's C-code for ubercalibration.

	09/06/2016 KMarkovic
"""
from subprocess import call
import dith, io

FTOL = 1e-6
DETBOOL=1

def test_calibration(starfile, seed=-1, calipath='./', outdir='./', outfile='./test-calibration.out', verb=False):
	""" Runs test-calibration.c to simulate a normally distributed set of calibrations 
		given a Mangle vertices file of the survey. It then finds the best fit and returns
		the zero-point scatter of the calibrated survey.

		Inputs
		------
		...

		Outputs
		-------
		ical - initial calibration scatter
		fcal - final calibration scatter
	"""
	# Run Will's code to simulate the calibrator statistics and do the ubercalibration test
	with open(outfile, "w") as cf:

		cmd = ['./test-calibration', outdir+'/', str(seed), str(starfile), str(FTOL), str(int(DETBOOL)), str(int(verb))]

		if verb: print "In " + calipath + ", run:\n\t" + ' '.join(cmd)
		try:
			call(cmd, stdout=cf, cwd = calipath);
		except OSError as e:
			raise Exception('\n\tThe following command failed:\n\t' + ' '.join(cmd) + '\n\tin folder:\n\t' + calipath + '\n\twith the error:\n\t' + str(e))
		else:
			if verb: print "Ran test-calibration."

	return get_cals(outfile, False, verb)


def get_cals(outfile, may_be_unfinished=False, verb=False):
	""" Read initial and final calibration scatter from the output file of test-calibration.c.

		Input
		-----
		outfile - path+filename to file produced by test-calibration

		Outputs
		-------
		ical - initial calibration scatter
		fcal - final calibration scatter
	"""

	# Get calibration results from last line of the output file
	with open(outfile, "r") as cf:		
		line = None
		for line in cf: 
			if verb: pass
		
	# The last line should be like: 
	#	"Initial calibration ical[i], final calibration fcal[i]"
	if line==None or not "final scatter" in line:
		if may_be_unfinished:
			print outfile + " not finished."
			return None, None
		else:
			raise Exception("File " + outfile + " is not what I expected! I got this line:\n" + str(line))

	# If all is as expected, extract the calibration values
	ical = float(line.split(", ")[0].split(" ")[2]) 
	fcal = float(line.split(", ")[1].split(" ")[2])

	return ical, fcal


def compile(calipath='./'):
	""" Compiles the C-code. Note that it does not compile Mangle, which is 
		assumed to be installed on your system.

		Inputs
		------
		calipath (optional) - path to the source code

		Outputs
		-------
		None
	"""

	print "\n--------- Compiling ubercal... ---------"
	call(['make','clean'], cwd = calipath)
	call(['make','create-euclid-patch'], cwd = calipath)
	call(['make','test-calibration'], cwd = calipath)
	print "--------- Compilation finished. ---------\n"

