""" Functions to do with the dithering pattern and survey strategy.
	Also to run Will's survey construction code.

	09/06/2016 KMarkovic
"""

from subprocess import call
import numpy as np

BASELINE = [50.0,100.0,0.0,100.0,0.0,100.0]


def getpat(pattern_name='S', dx=50, dy=100):
	""" The dither "vector" as a 2d Numpy array given the pattern name

		Input
		-----
		pattern_name - a string name (choose from S, J, N, ...)
		dx - x-displacements in arcsec
		dy - y-displacements in arcsec

		Output
		------
		dithers - a 2d numpy array of dither vectors
	"""

	# Build up an array of dither strategy multiples
	# 	it is in arcsec (n.b. one CCD is 2040 pixel = 612", gap about 50", 100")
	if pattern_name=="J":
		end_dithers = np.array([[dx,dy],[0.0,dy],[0.0,dy]])
	elif pattern_name=="Jsq":
		end_dithers = np.array([[2.*dx,dy],[0.0,dy],[0.0,dy]])
	elif pattern_name=="O" or pattern_name=="box":
		end_dithers = np.array([[dx,0.0],[0.0,dy],[-dx,0.0]])
	elif pattern_name=="R" or pattern_name=="step":
		end_dithers = np.array([[dx,dy],[dx,0.0],[dx,0.0]])
	elif pattern_name=="S":
		end_dithers = np.array([[dx,dy],[0.0,dy],[dx,dy]])
	elif pattern_name=="Ssq":
		end_dithers = np.array([[2.*dx,dy],[0.0,dy],[2.*dx,dy]])
	elif pattern_name=="X":
		end_dithers = np.array([[dx,dy],[0.0,0.0],[-dx,-dy]])
	elif pattern_name=="N":
		end_dithers = np.array([[dx,dy],[dx,0.0],[dx,dy]])
	else:
		raise RuntimeError("I don't know what the " + pattern_name+ "-pattern is!\n")

	return end_dithers


def totaldither(dv):
	""" Calculates the total distance from first to last dither and the
		total telescope path.

		Input
		-----
		dv - dither vector as a list

		Outputs
		-------
		r - total distance
		d - total path
	"""
	# Calculate the vector magnitude of 4th dither displacement from 1st
	r = np.sqrt((dv[0]+dv[2]+dv[4])**2 + (dv[1]+dv[3]+dv[5])**2)
	# Calculate full distance travelled by telescope
	d = np.sqrt(dv[0]**2 + dv[1]**2) + np.sqrt(dv[2]**2 + dv[3]**2) + np.sqrt(dv[4]**2 + dv[5]**2)
	return r, d


def create_surveypatch(dv=BASELINE, nx=3, ny=3, calipath='./', outdir='./', outfile='./create-surveypatch.out', verb=False):
	""" Runs create-euclid-patch.c to create Mangle vertices files for a small
		Euclid-like survey on a flat sky.

		Inputs
		------
		dv - dither vector as a list of x and y dither displacements, e.g. [x1, y1, x2, y2...], default is baseline
		nx - number of pointings in the x-direction
		ny - number of pointintgs in the y-direction
		calipath - path to ubercal source c-code
		outdir - path to save the vertices files to
		outfile - path+filename of log file

		Outputs
		-------
		area - ...
		frac - ...
	"""

	# Create the survey, save the Mangle files
	with open(outfile, "w") as pf:
			
		cmd = ['./create-euclid-patch', outdir+'/', str(dv[0]), str(dv[1]), str(dv[2]), str(dv[3]), str(dv[4]), str(dv[5]), str(nx), str(ny)]
		
		if verb: print "In " + calipath + ", run:\n\t" + ' '.join(cmd)
		try:
			call(cmd, stdout=pf, cwd = calipath);
		except OSError as e:
			raise Exception('\n\tThe following command failed:\n\t' + ' '.join(cmd) + '\n\tin folder:\n\t' + calipath + '\n\twith the error:\n\t' + str(e))
		else:
			if verb: print "Ran create-euclid-patch."

	return get_area(outfile, verb)


def get_area(outfile, verb=False):
	""" Extract the area from the survey patch output file. 

	Input
	-----
	outfile - path+filename to file where create-euclid-patch output was saved

	Output
	------
	area - ...
	frac - ...

	"""

	line = None
	area = None
	with open(outfile, "r") as pf:
		
		frac=[]
		for line in pf: 

			try:
				if "coverage" in line and "-passes" in line:
					frac += [float(line.split(" = ")[-1][:-2])/100.0]
				elif "tot" in line and "area" in line and "=" in line:
					area = float(line.split("= ")[-1])
				else:
					pass
			except ValueError as e:
				raise Exception("This line makes no sense:\n'" + line[:-1] + "'\n" + str(e))
	
	if line is None or area is None: raise Exception("File " + outfile + " is empty!")

	return area, frac
