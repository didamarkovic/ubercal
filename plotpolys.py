""" Plot .vrt file from the ubercal code by Will Percival.

	29/07/2015, Dida Markovic
"""

import numpy as np
import matplotlib
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, colorConverter
import matplotlib.pyplot as plt
from itertools import izip

DEF_SCALE = 3600.0 # Mangle outputs will likely be in arcseconds, but want to convert to degrees
PRECISION = 1e-15

class Poly(object):
	""" Contains the properties of each set of vertices. """

	def __init__(self, inlist=None, COORD_SCALE=DEF_SCALE, XZERO=0.0, YZERO=0.0):
		""" Reading in from the first lines of the .vrt file. """
		self.COORD_SCALE = COORD_SCALE
		self.XZERO = XZERO
		self.YZERO = YZERO
		if len(inlist) == 5:
			self.ind = int(inlist[0])
			self.nverts = int(inlist[1])
			self.weight = float(inlist[2])
			self.centre = [float(inlist[3])/self.COORD_SCALE-self.XZERO, float(inlist[4])/self.COORD_SCALE-self.YZERO]
		elif len(inlist) == 2:
			self.ind = int(inlist[0])
			self.weight = float(inlist[1])
			self.nverts = 0
			self.centre = None
		else:
			self.ind = 0
			self.nverts = 0
			self.weight = 1
			self.centre = None

	def vertices(self, coolist):
		""" Get x,y coordinates of vertices from a list of even length. """
		if self.nverts==0:
			self.nverts = len(coolist)/2
		self.vrt = np.array(coolist).astype(float).reshape(self.nverts,2)/self.COORD_SCALE - [self.XZERO,self.YZERO]
		self.xvrt, self.yvrt = zip(*self.vrt)
		# Calculate these if they were not input
		# Assume polygons are simple enough that the average of the vertices gives a good estimate of the centroid
		if self.centre is None:
			self.centre = [sum(self.xvrt)/self.nverts, sum(self.yvrt)/self.nverts]

	def borderx(self):
		""" Repeat the first entry for a plotted line to go all the way around. """
		return self.xvrt + tuple([self.xvrt[0]])

	def bordery(self):
		""" Repeat the first entry for a plotted line to go all the way around. """
		return self.yvrt + tuple([self.yvrt[0]])

	def includein(self, minx, maxx, miny, maxy):
		""" Return the outermost edges, input or polygon. """
		if min(self.xvrt) < minx: minx = min(self.xvrt)
		if max(self.xvrt) > maxx: maxx = max(self.xvrt)
		if min(self.yvrt) < miny: miny = min(self.yvrt)
		if max(self.yvrt) > maxy: maxy = max(self.yvrt)
		return minx, maxx, miny, maxy

	def origin(self):
		""" Return the bottom left vertex. """
		return [min(self.xvrt), min(self.yvet)]

	def width(self):
		""" Return width. """
		return max(self.xvrt) - min(self.xvrt)

	def height(self):
		""" Return height. """
		return max(self.yvrt) - min(self.yvrt)

	def area(self):
		assert hasattr(self, 'vrt')
		return self.height()*self.width()

	def __str__(self):
		""" How to print out. """
		info = "Poly object no. " + str(self.ind)
		info += "\n\t number of vertices = " + str(self.nverts)
		info += "\n\t weight = " + str(self.weight)
		info += "\n\t centered at " + str(self.centre)
		if hasattr(self, 'vrt'):
			info += "\n\t polygon area = " + str(self.area())
		return info

	def reweight(self, weight):
		self.weight = weight

	def parents(self, parentlist):
		pass

def add_Polys(axes, polygons, col='g', maxa=1.0, weights=[], lw=0., verts=False, lc=None):
	""" Add a list of Poly objects to axes and return the unique weights for the legend """

	# Initial checks of inputs
	assert len(weights) <= len(polygons)
	assert maxa>=0 and maxa<=1

	# Which polygon weights to use
	if len(weights)==len(polygons):
		replaceweights = True
		all_weights = weights
	else:
		replaceweights = False
		all_weights = [polygons[i].weight for i in range(len(polygons))]
		
	# Whether to indicate weight by colour or alpha level
	if isinstance(weights, dict): 
		colors = True
	else:
		colors = False
		pass_range = float(max(all_weights)) + PRECISION # Normalise the alpha range
	
	forlegp = []
	forlegc = []
	donew = []
	i = 0
	for p in polygons:

		# Set colour depending on options
		if replaceweights:
			weight = weights[i]
			i+=1
		else:
			weight = p.weight
		if colors:
			try:
				fcol = weights[weight]
			except KeyError:
				weight = -1
		else:
			fcol = colorConverter.to_rgba(col, weight/pass_range*maxa)

		if weight == -1: continue

		if lc is None: 
			ec = fcol
		else:
			ec = lc
			
		# Generate polygons and add it to axes
		pp = Polygon(p.vrt, linewidth=lw, fc=fcol, ec=ec)
		axes.add_patch(pp)

		# Plot the vertices
		if verts is not False:
			verts = pp.get_xy().T
			axes.plot(verts[0],verts[1],col+'.')

		# Keep track of which weight values were added for the legend
		if p.weight not in donew:
			forlegp.append(pp)
			forlegc.append(int(p.weight))
			donew.append(p.weight)
	
	# Return sorted legend entries
	sort = np.argsort(forlegc)
	return list(np.array(forlegp)[sort]), list(np.array(forlegc)[sort])

def outline_strange(axes, polygons, cols='bry'):
	""" Outline polygons with other than 4 vertices """
	for p in polygons:
		if p.nverts!=4:
			axes.plot(p.borderx(),p.bordery(), cols[0]+'--', ms=3)
		if p.nverts==3:
			axes.plot(p.centre[0], p.centre[1], cols[1]+'.', ms=10)
		elif p.nverts==5:
			axes.plot(p.centre[0], p.centre[1], cols[2]+'.', ms=10)

def read_vrt(filename, verb=True, silent=False, SCALE=DEF_SCALE):
	""" Open the vertex file, read line by line
		Expecting a standard format for these .vrt files => remove the unneeded info
		Read the polygons into vertex arrays """	
	polygons = []
	minx = float('inf'); maxx = 0.0
	miny = float('inf'); maxy = 0.0
	with open(filename, 'r') as f:			
		""" Read one line of the vertex file. """
		l = 0
		for line in f:
			line = line.strip().split()
			if "polygons" in line:
				npolys = int(line[0])
			elif l==0 and "vertices" not in line:
				if verb and not silent: print ' '.join(line)
			elif 'vertices' in line:
				for word in ['vertices', 'vertices,', '(', 'weight,', 'mid):']: line.remove(word);
				polygons.append(Poly(line, SCALE))
				l+=1
			elif 'vertices' not in line:
				polygons[-1].vertices(line)
				minx, maxx, miny, maxy = polygons[-1].includein(minx, maxx, miny, maxy)
				l+=1
	assert len(polygons)>0, "This file does not contain Mangle vertices: " + filename + "!"
	return minx, maxx, miny, maxy, npolys, polygons

def test_fullsample(vrtfile='./test-calibration/full-pointing.vrt', txtfile='./test-calibration/full-survey-overlaps.txt', SCALE=DEF_SCALE):
	""" Test the vrt and txt out put files of creat-euclid-fullsample.c. Makes many assumptions
		about the format. Plots the polygons that are not discarded in green. In pink it plots those
		that are discarded by the code. The dashed outlines (dots in centre) are those that seem to have 
		3 or 5 vertices (not 4), which is unexpected. """

	minx, maxx, miny, maxy, npolys, polygons = read_vrt(vrtfile, SCALE=SCALE)

	""" Open the file containing areas and parent polygons. 
	Get the indices of the final polygons to only plot those. """
	with open(txtfile, 'r') as f:
		f.readline(); f.readline();
		npolys_out = int(f.readline().strip().split()[1])
		inds_out = []*npolys_out
		weights_out = []*npolys_out
		for line in f:
			line = line.split()
			inds_out.append(int(line[0]))
			weights_out.append(int(line[2]))

	discard = npolys!=len(inds_out)
	discarded = None
	if discard:
		assert len(polygons)>=max(inds_out), "There are more kept polygons than the input."

		""" Discard the discarded """
		discarded = [p for p in polygons if p.ind not in inds_out] 
		polygons = [polygons[i] for i in inds_out]
		for i in range(len(weights_out)):
			polygons[i].reweight(int(weights_out[i]))
	else:
		for i in range(len(weights_out)):
			polygons[inds_out[i]].reweight(int(weights_out[i]))

	return polygons, discarded, minx, maxx, miny, maxy

def plotpolys(polygons, discard=False, minx=np.NAN, maxx=np.NAN, miny=np.NAN, maxy=np.NAN):
	""" Make a nice plot of the polygons. """

	fig, ax = plt.subplots()
	forlegp, forlegc = add_Polys(ax, polygons, maxa=0.99, lw=0.)
	print "Plotted " + str(len(polygons)),
	
	if discard:
		flpd, flcd = add_Polys(ax, discarded, 'm', 0.2, lw=0.)
		forlegp = forlegp + flpd
		forlegc = forlegc + flcd
		polygons = polygons + discarded
		print "kept polygons and " +str(len(discarded)) + " discarded ones."
	else:
		print "polygons."
	
	plt.legend(forlegp, forlegc, bbox_to_anchor=(1.3, 1.0), title='weight')
	outline_strange(ax, polygons, 'mmc')

	""" Make plot pretty """
	ax.set_aspect('equal', adjustable='box')
	if not np.isnan(minx+maxx+miny+maxy): plt.xlim([minx*1.1-maxx*0.1, maxx*1.1 - minx*0.1]); plt.ylim([miny*1.1-maxy*0.1, maxy*1.1]-miny*0.1)
	plt.ticklabel_format(useOffset=False)
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')

if __name__=="__main__":

	import argparse, time, os.path

	parser = argparse.ArgumentParser(description="Plot Mangle output, optionally compare to ubercal output.")
	parser.add_argument("vrt", help="file containing Mangle vertices: the 'vertices # ( # vertices... )' format")
	parser.add_argument("-t", "--txt", default=None, help="full-survey-overlaps.txt file from ubercal")
	parser.add_argument("-s", "--save", default=None, help="folder or file to save the plot to, default shows the figure instead of saving")
	args = parser.parse_args()

	# Read from file and create polygons
	if args.txt is not None: 
		polygons, discarded, minx, maxx, miny, maxy = test_fullsample(vrtfile=args.vrt, txtfile=args.txt, SCALE=DEF_SCALE)
	elif args.vrt is "default":
		polygons, discarded, minx, maxx, miny, maxy = test_fullsample(SCALE=DEF_SCALE)
	else:
		minx, maxx, miny, maxy, npolys, polygons = read_vrt(args.vrt, SCALE=DEF_SCALE ,verb=False)
		discarded = None

	# Make the plot
	plotpolys(polygons, discarded is not None, minx, maxx, miny, maxy)

	# Show or save
	if args.save is None:
		plt.show();
	elif os.path.isdir(args.save):
		ts = time.strftime('%s')
		plt.savefig(args.save+"plotpolys-" + ts +'.pdf',dpi=400,bbox_inches='tight') 
		plt.close()
	else:
		plt.savefig(args.save+'.pdf',dpi=400,bbox_inches='tight') 
		plt.close()
