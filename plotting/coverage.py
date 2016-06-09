"""
	Calculate the n-pass coverage for a given dither pattern of a given size.
	First attempt did this with binary pixel-level masks, but this is very 
	inefficient => write as overlapping squares.

	TODO?: perhaps the plotpolys.Poly class might be useful here?

	KMarkovic, 
	24/11/2014, v0.0
	21/10/2015, v1.0
"""

import numpy as np
import subprocess as sp
import os.path as op
from glob import glob
from warnings import warn

import matplotlib.pyplot as plt
import plotpolys as pp

# Instrument model
NDETX = 4
NDETY = NDETX
XGAP = 50.0
YGAP = 100.0
XDET = 612.0
YDET = XDET
SCALE = 3600.0
PRECISION = 0.1 # in arcsec

PATTERNS = {
	'J': np.array([   XGAP,YGAP,   0.0,YGAP,   0.0, YGAP]),
	'S': np.array([   XGAP,YGAP,   0.0,YGAP,  XGAP, YGAP]),
	'step': np.array([XGAP,YGAP, XGAP,  0.0,  XGAP,  0.0]), # deprecate
	'R': np.array([XGAP,YGAP, XGAP,  0.0,  XGAP,  0.0]),
	'N': np.array([   XGAP,YGAP, XGAP,  0.0,  XGAP, YGAP]),
	'X': np.array([   XGAP,YGAP,  0.0,  0.0, -XGAP,-YGAP]),
	'box': np.array([ XGAP, 0.0,  0.0,100.0, -XGAP,  0.0]), # deprecate
	'O': np.array([ XGAP, 0.0,  0.0,100.0, -XGAP,  0.0])
	}

# Path for execution (either input or find)
HERE = op.dirname(op.abspath(__file__)) + '/'

def convpat(pat):
	if pat=='box': 
		pat='O'
	elif pat=='step':
		pat='R'
	return pat

def build(path, target, verb=False):
	""" Build all the targets in the path, assuming there is a Makefile present there."""

	# Use subprocess to run the make command and keep track of the output and errors
	cmd = ['make', target]
	child = sp.Popen(cmd, cwd = path, stdout = sp.PIPE, stderr = sp.PIPE)
	out, err = child.communicate()
	
	# Raise an error if the build was unsuccessful
	if err:
		raise Exception("\n\t'"+' '.join(cmd)+"'\ncalled from\n\t"+path+"\nexited with the following error:\n"+err)
	elif verb>1:
		print out

	# Return the full path to the binary by reading the output
	try:
		out = [string for string in out.split() if string[0] is not '-' and string[-2] is not '.'][-1]
		out = op.abspath(op.join(path,out))
	except:
		raise Warning("The compiler statement can't be clearly read to find the binary:\n" + out)
	else:
		return out

def create_euclid_patch(dithvec, outpath, nra=1, ndec=1, ndetx=4, binary='create-euclid-patch', verb=False):
	""" Call the ubercal c-code to construct the Euclid mask of pointings and dithers. """

	# Construct the call with the correct input arguments
	bindir, binfile = op.split(binary)
	cmd = ['./'+binfile,outpath] + [str(x) for x in dithvec] + [str(nra), str(ndec), str(ndetx)]
	if verb>1: print ' '.join(cmd)

	# Check that the binary exists
	if not op.isfile(binary):
		raise Exception('Binary ' + binary + " doesn't exist. Perhaps you forgot to run make.")

	# Use subprocess to run the Mangle-calling, mask-creating code
	try:
		child = sp.Popen(cmd, cwd = bindir, stdout = sp.PIPE, stderr = sp.PIPE)
	except OSError:
		raise Exception("I can't find the executable: " + binary)
	out, err = child.communicate()

	# Raise an error if the build was unsuccessful
	if err:
		raise Exception("\n\t'"+' '.join(cmd)+"'\ncalled from\n\t"+bindir+"\nexited with the following error:\n"+err)
	elif verb:
		print out

	out = [outpath+'full-survey.vrt', outpath+'full-survey-overlaps.vrt']
	out.sort(key=len)
	return out

def read_vrt(filename, verb=False, COORD_SCALE=SCALE):
	""" Open the vertex file, read line by line
		Expecting a standard format for these .vrt files => remove the unneeded info
		Read the polygons into vertex arrays """	
	
	minx = float('inf'); maxx = 0.0
	miny = float('inf'); maxy = 0.0
	
	with open(filename, 'r') as f:			
		""" Read one line of the vertex file. """
		
		l = 0
		r = 0
		nonrec = False
		for line in f:
			line = line.strip().split()
			
			if l==0 and "polygons" in line:
				nrec = int(line[0])
				# Create an empty array with nrec rows and 10 columns for:
				#  ID weight dx1 dy1 dx2 dy2 dx3 dy3 dx4 dy4
				rectangles = np.zeros((nrec,10))*np.NAN
			
			elif l==0 and "vertices" not in line:
				if verb>1: print ' '.join(line)
			
			elif 'vertices' in line:
				for word in ['vertices', 'vertices,', '(', 'weight,', 'mid):']: line.remove(word);
				rectangles[r,0] = int(line[0])
				rectangles[r,1] = float(line[2])
				l+=1
			
			elif 'vertices' not in line:
				verts = np.array(line).astype(float)/COORD_SCALE
				try:
					rectangles[r,2:] = verts
				except ValueError as e:
					if not nonrec:
						if verb: print "Non-rectangles found in your .vrt file!\n" + '\t' + filename
						if verb: print "Squashing those of the following that have nverts > 4 and leaving those that have nverts < 4:"
						nonrec = True
					rectangles[r,2:] = union_rectangle(verts)
					if verb: print "Polygon no."+str(int(rectangles[r,0]))+" has " + str(len(verts)/2) + " vertices, squash to: " + '[' + ','.join(format(x, ".3f") for x in rectangles[r,2:]) + ']'
					
				l+=1
				r+=1
	
	# Check that all the array entries have been filled correctly
	check = np.isnan(rectangles)
	assert np.sum(check)==0, "I'm having trouble reading rectangles: \n" + str(np.unique(check.nonzero()[0])) + "\n in " + filename + "!" 
	
	return rectangles

def area_rec(rectangles):
	""" Return an array of areas for the input rectangles (np.array). """
	ndim = len(rectangles.shape)
	if ndim==1:
		width = np.max(rectangles[::2]) - np.min(rectangles[::2])
		height = np.max(rectangles[1::2]) - np.min(rectangles[1::2])	
	elif ndim==2:
		width = np.max(rectangles[:,::2],1) - np.min(rectangles[:,::2],1)
		height = np.max(rectangles[:,1::2],1) - np.min(rectangles[:,1::2],1)
	else:
		raise Exception("Input must be an array of rows of vertices => a 2d numpy array!")
	return np.array(width*height)

def snap(rectangles, P=PRECISION/SCALE, verb=False):
	""" Make coordinates equal if they are within a certain distance of each other. """
	return np.round(rectangles/P)*P

def overlaps(rec1, rec2, verb=False):
	""" Return a boolean (or an array of them if 2-d array is input for rec1).
		True if the given rectangles overlap. rec2 should be a 1-d array of 
		vertices coords. The rectangle arrays should only contain the vertices!"""

	x1 = rec1[:,::2]
	y1 = rec1[:,1::2] 
	x2 = rec2[::2]
	y2 = rec2[1::2]

	# Check which rec1[i] is partially inside rec2
	xmin  = x1 > min(x2)
	xmax  = x1 < max(x2)
	ymin  = y1 > min(y2)
	ymax  = y1 < max(y2)
	xy = xmin*ymin*xmax*ymax
	rec1_in_rec2 = np.sum(xy,1)
	if verb>2: print "1 in 2: "+str(rec1_in_rec2)

	# Check if rec2 is partially inside any of the rec1[i]
	x2 = np.logical_and(np.matrix(np.min(x1,1)).T < x2, np.matrix(np.max(x1,1)).T > x2)
	y2 = np.logical_and(np.matrix(np.min(y1,1)).T < y2, np.matrix(np.max(y1,1)).T > y2)
	xy = np.array(np.logical_and(x2,y2))
	rec2_in_rec1 = np.sum(xy,1)
	if verb>2: print "2 in 1: "+str(rec2_in_rec1)

	# Return the row numbers where all the above are true at least once
	# 	i.e. the element of rec1 where at least one vertex is inside rec2 or v.v.
	olaps = list(np.where(rec1_in_rec2 + rec2_in_rec1)[0])

	# Finally add those that completely identical to rec2
	if verb>2: print _str(rec2) + " <- c.f. to this"
	for i in range(len(rec1)):
		if verb>2:
			print "testing: " + str(i)
			print _str(rec1[i])
		if i not in olaps:
			if inside(rec1[i], rec2):
				olaps = olaps + [i]
				if verb>2: print "\t Fully inside!"
			elif inside(rec2, rec1[i]):
				olaps = olaps + [i]
				if verb>2: print "\t Fully inside!"
			else:
				if verb>2: print "\t Far away."
		else:
			if verb>2: print "\t Partial overlap."

	return sorted(olaps)

def union_rectangle(rectangles):
	""" Returns the all-encompassing ractangle, assuming edges of all input
		rectangles are aligned with the x-y axes. """
	ndim = len(rectangles.shape)
	if ndim==1:
		x1, x2 = [np.min(rectangles[::2]), np.max(rectangles[::2])]
		y1, y2 = [np.min(rectangles[1::2]), np.max(rectangles[1::2])]
	elif ndim==2:
		x1, x2 = [np.min(rectangles[:,::2]), np.max(rectangles[:,::2])]
		y1, y2 = [np.min(rectangles[:,1::2]), np.max(rectangles[:,1::2])]
	else:
		raise Exception("Input must be an array of rows of vertices => a 2d numpy array!")
	return np.array([x1,y1, x2,y1, x2,y2, x1,y2])

def intersect_rectangle(rectangles, verb=False):
	""" Returns the intersection rectangle, assuming edges of all input
		rectangles are aligned with the x-y axes. """
	x1, x2 = [np.max(np.min(rectangles[:,::2],0)), np.min(np.max(rectangles[:,::2],0))]
	y1, y2 = [np.max(np.min(rectangles[:,1::2],0)), np.min(np.max(rectangles[:,1::2],0))]
	intersect = np.array([x1,y1, x2,y1, x2,y2, x1,y2])
	if verb>2: print "Intersect: " + _str(intersect)
	if len(overlaps(rectangles,intersect, verb=verb))==len(rectangles):
		return intersect
	else:
		return None

def full_coverage(rectangles, verb=False):
	""" Return the 4 vertices of the maximum-coverage area. """

	# Divide the rectangles into N_dith dithers, assuming they are ordered already
	N_dith = 4
	nrec = len(rectangles)
	assert nrec%N_dith==0, "Your input survey mask seems not to have "+str(N_dith)+" dithers!"
	en = int(nrec/N_dith)

	# Find the extrema of each survey dither
	verts = np.zeros((N_dith,8))
	for i in range(N_dith):
		one_survey_dither = rectangles[(i*en):((i+1)*en)]
		verts[i] = union_rectangle(one_survey_dither) + np.array([0.0,0.0, XGAP,0.0, XGAP,YGAP, 0.0,YGAP])/SCALE

	# Now return the intersection
	return intersect_rectangle(verts)

def central_pointing(origin, COORD_SCALE=SCALE):
	""" Given the survey parametres, return the 4 vertices of the "square of influence"
		of the central pointing: i.e. pointing + gap. """
	xpoint = np.array([NDETX*(XGAP+XDET), 0.0])/COORD_SCALE
	ypoint = np.array([0.0, NDETY*(YGAP+YDET)])/COORD_SCALE
	origin = np.array(origin) + xpoint + ypoint
	return np.array([origin, origin+xpoint, origin+xpoint+ypoint, origin+ypoint]).reshape(8)

def inside(rec1, rec2, P=PRECISION/SCALE):
	""" Return true if rec1 is fully inside rec2. """

	# x-values
	# c.f. minima
	xmin = min(rec1[::2]) < min(rec2[::2]) - P
	# c.f. maxima
	xmax = max(rec1[::2]) > max(rec2[::2]) + P

	# y-values
	# c.f. minima
	ymin = min(rec1[1::2]) < min(rec2[1::2]) - P
	# c.f. maxima
	ymax = max(rec1[1::2]) > max(rec2[1::2]) + P

	return not (xmin or xmax or ymin or ymax)

def _str(r):

	try:
		dims = len(r.shape)
	except AttributeError:
		return str(r)

	if r.shape[-1]==10:
		meta=True
	elif r.shape[-1]==8:
		meta=False
	else:
		raise Exception("Rectangle should either have 8 or 8+2 elements, instead it has " + str(r.shape[-1]))

	if dims==1 and meta:
		return "#%+03d, w=%+d at (%+03.3f,%+03.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f)" % tuple(r)
	elif dims==1:
		return "(%+03.3f,%+03.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f)" % tuple(r)
	elif meta:
		out = ''
		for i in r:
			out += "#%+03d, w=%+d at (%+03.3f,%+03.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f)\n" % tuple(i)
		return out
	else:
		out = ''
		for j,i in zip(range(r.shape[1]),r):
			out += "#%+03d at (%+03.3f,%+03.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f), (%+03.3f,%+02.3f)\n" % tuple([j]+list(i))
		return out

def crop_survey(rectangles, cropvec, verb=False, drop=False):
	""" Only keep the survey area of 4-passes as if there were no chip and pointing gaps.
		This will break if we are on a sphere! E.g. if RA ~ 0! """
	
	if verb>2:
		print "Crop at: " +str(cropvec)
		print "Output: " + str(rectangles)

	keep = overlaps(rectangles[:,2:], cropvec, verb=verb)
	if verb>2: print '\n\t' + _str(np.array([-1,-1]+list(cropvec))) + '\n'

	for i, rec in zip(range(len(rectangles)),rectangles):
		if inside(rec[2:],cropvec): 
			continue
		elif i in keep:
			rec[2:] = intersect_rectangle(np.vstack([rec[2:],cropvec]))
		else:
			rec[1] = 0
	if verb>2: print _str(rectangles)

	if drop:
		rectangles = rectangles[rectangles[:,1]!=0]

	return rectangles

def coverage(rectangles, tot_lims=None):
	""" Add the areas of rectangles with non-zero weights. """

	no_w = max(rectangles[:,1])+1
	covperpass = np.zeros(no_w)

	for rec in rectangles:
		width = max(rec[2::2]) - min(rec[2::2])
		height = max(rec[3::2]) - min(rec[3::2])
		covperpass[rec[1]] += width*height

	# Now get total if total limits are given
	if tot_lims is not None:
		width = max(tot_lims[::2]) - min(tot_lims[::2])
		height = max(tot_lims[1::2]) - min(tot_lims[1::2])
		total = width*height
		# The rectangles that were intentionally set to have 0 weight must be ignored!
		covperpass[0]=total-sum(covperpass[1:])
	else:
		total = 1.0
		# The rectangles that were intentionally set to have 0 weight must be ignored!
		covperpass[0]=0.0

	return covperpass, total

def survey_coverage(dithvec=PATTERNS['J'], binary='./bin/create-euclid-patch', outpath='./', verb=0):
	""" Do the full calculation from input parameters to coverage. """

	# Always get a 3x3 survey so that can examine the central part
	npoint = 3
	# Note that this puts a restriction on how big the size parameter can be

	# Run the create-euclid-patch.c ubercal code to construct the Mangle mask
	inv, outv = create_euclid_patch(dithvec, outpath, nra=npoint, ndec=npoint, ndetx=NDETX, binary=binary, verb=verb)

	# Read the Mangle in- & output file
	inv = read_vrt(inv, verb=verb)
	outv = read_vrt(outv, verb=verb)
	origin = (np.min(inv[:,2::2]), np.min(inv[:,3::2]))

	# Find patch of full coverage
	#lim_rec = full_coverage(inv[:,2:], verb=arg.verb)
	lim_rec = central_pointing(origin)

	# Identify the edges that are close
	tmp = snap(np.vstack([lim_rec, outv[:,2:]]))
	lim_rec = tmp[0]; outv[:,2:] = tmp[1:]; tmp = None

	# Now throw away the polygons that not in the full coverage patch of the survey
	outv = crop_survey(outv, lim_rec, verb=verb)

	return outv, lim_rec

def _test(outv, lim_rec, dithv, verts=False):

	l = lim_rec.reshape(4,2)
	dithv = np.cumsum(np.array([0,0] + list(dithv)).reshape(4,2),0) + l[0]
	l = zip(*l)
	dithv = zip(*dithv)
	
	nout = outv.shape[0]
	polys = [0]*nout
	for i,poly in zip(range(nout), outv):
		polys[i] = pp.Poly(outv[i][0:2], 1.0)
		polys[i].vertices(outv[i,2:])
	limp = pp.Poly([0,1], 1.0)
	limp.vertices(lim_rec)

	fig, ax = plt.subplots()
	pp.add_Polys(ax, polys, maxa=0.99, lw=1., verts=verts)
	pp.add_Polys(ax, [limp], maxa=0.0, lw=5., verts=verts, col='m')
	plt.plot(dithv[0], dithv[1], 'k.-')
	plt.xlim([min(l[0])*1.1-max(l[0])*0.1, max(l[0])*1.1 - min(l[0])*0.1]); 
	plt.ylim([min(l[1])*1.1-max(l[1])*0.1, max(l[1])*1.1]- min(l[1])*0.1)
	ax.set_aspect('equal', adjustable='box')
	plt.ticklabel_format(useOffset=False)
	plt.xlabel("RA [deg]")
	plt.ylabel("dec [deg]")
	
if __name__=='__main__':

	import argparse
	p = argparse.ArgumentParser(description="Generate the Euclid survey geometry and calculate n-pass coverage.")
	p.add_argument("--pattern", "-p", default="S", choices=["J", "step", "S", "N", "X", "box"], help="Name of dither pattern.")
	p.add_argument("--size", "-s", default=1.0, type=float, help="Scale of dither, SIZE=1 means d_x=50'', d_y=100''.")
	p.add_argument("-o", "--outpath", default=op.join(HERE,"outputs/"), help="where you want your output files")
	p.add_argument("-c", "--srcpath", default=op.join(HERE,"ubercal/"), help="directory containing test-calibration source code")
	p.add_argument("-v", "--verb", default=0, type=int, help="verbosit y level. 0 for silent.")
	p.add_argument("-b", "--bin", default=None, help="Location of binaries if no need for compilation.")
	arg = p.parse_args()

	# Since only calculating a 3x3 survey, need to warn if size of dither is too big
	if arg.size > 2.0:
		warn("The coverage calculation is not guaranteed to work for very large dithers!")

	# Make sure we have built the latest version of the c-code
	target = 'create-euclid-patch'
	if not arg.bin:
		binary = build(path=arg.srcpath, target=target, verb=arg.verb)
	else:
		binary = op.join(arg.bin,target)

	# Survey parameters
	dithvec = arg.size*PATTERNS[arg.pattern]

	# Run code
	outv, lim_rec = survey_coverage(dithvec, binary, arg.outpath, arg.verb)

	# Calculate coverage
	cov, tot = coverage(outv, lim_rec)
	
	for p,f in zip(range(len(cov)),cov/tot*100):
		print "%d-pass coverage area fraction: %02.2f%%" % (p,f)

	# Test this module if verbosity is at debug level (i.e. 2)
	if arg.verb==2:
		_test(outv, lim_rec, dithvec/SCALE)
		plt.title(arg.pattern+'-pattern' + ', dx = '+str(round(dithvec[0],1)) + "'' & dy = "+str(round(dithvec[1],1))+"''")
		plt.show()	

