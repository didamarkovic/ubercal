""" plotsetup.py
	Make a nice plot of the survey, dither and coverage setup.
	KMarkovic, 09/12/2015
"""

import numpy as np
import time as t
import os
from matplotlib import rc
import matplotlib.pyplot as plt
import plotpolys as pp
import coverage
import os.path as op

HERE = op.dirname(op.abspath(__file__)) + '/'
NPOINT = 3
NDITH = 4
NDETX = 4

def _plot_survey(inv, outv, oned, NRA=3, NDEC=3, outfile=None, colours=False, ts=None, scale=coverage.SCALE, unit=r'$\degree$', fs=20):
	""" Make a nice plot of the dither configuration and coverage. """

	fig, ax = plt.subplots(figsize=(10,10))

	# Create and plot input 4-dither polygons
	nin = inv.shape[0]
	ipolys = [0]*nin
	for i,poly in zip(range(nin), inv):
		ipolys[i] = pp.Poly(inv[i][0:2], scale)
		ipolys[i].vertices(inv[i,2:])
	pp.add_Polys(ax, ipolys, maxa=0.0, lw=1., verts=False, col='w', lc='0.8')

	# Create and overplot a 1-dither survey
	oin = oned.shape[0]
	opolys = [0]*(oin-16)
	ii=0
	for i,poly in zip(range(oin), oned):
		if int(i/16)==4: continue
		opolys[ii] = pp.Poly(oned[i][0:2], scale)
		opolys[ii].vertices(oned[i,2:])
		ii+=1
	pp.add_Polys(ax, opolys, maxa=0, lw=1., verts=False, col='w', lc='0')
	
	# Create and plot the central pointing coverage polygons
	nout = outv.shape[0]
	polys = [0]*nout
	for i,poly in zip(range(nout), outv):
		polys[i] = pp.Poly(outv[i][0:2], scale)
		polys[i].vertices(outv[i,2:])
	if colours:
		leg = {1:'r',2:'y',3:'b',4:'g'}
		forlegp, forlegc = pp.add_Polys(ax, polys, weights=leg, maxa=0.99, lw=0., verts=False)
	else:
		forlegp, forlegc = pp.add_Polys(ax, polys, maxa=0.99, lw=0., verts=False)

	# Set plot limits that contain all these polygons (Python fails doing it automatically)
	l = zip(*coverage.union_rectangle(inv[:,2:]/scale).reshape(4,2))
	xlim = [min(l[0])*1.01-max(l[0])*0.01, max(l[0])*1.01 - min(l[0])*0.01]
	ylim = [min(l[1])*1.01-max(l[1])*0.01, max(l[1])*1.01]- min(l[1])*0.01
	plt.xlim(xlim); 
	plt.ylim(ylim)

	# Put pointing numbers in
	for i in range(9):
		x = (oned[i*16+15,4] + 125.)/scale
		y = (oned[i*16+15,5] + 50.)/scale
		if i==4:
			c='w'
		else:
			c='k'
		plt.text(x, y, str(i), family='serif', color=c, fontsize='45')

	# Plot the first detector dither corners
	fulldith = inv.shape[0]/NDITH
	fullpoint = inv.shape[0]/NDITH/NPOINT**2
	for Npoint in range(NPOINT**2):
		cornerx = [0]*NDITH
		cornery = [0]*NDITH
		for i in range(NDITH):
			pointind = Npoint*fullpoint + i*fulldith
			firstdet = inv[pointind:pointind+fulldith][0]/scale
			cornerx[i] = firstdet[4]
			cornery[i] = firstdet[5]
		if Npoint==4:
			style = 'w.-'
		else:
			style = 'k.-'
		plt.plot(cornerx, cornery, style)


	ax.set_aspect('equal', adjustable='box')
	plt.ticklabel_format(useOffset=False)
	plt.ylabel(r'$\delta$ [' + unit + ']')
	plt.xlabel(r'RA [' + unit + ']')
	#plt.legend(forlegp, forlegc, title='# passes', bbox_to_anchor=(1.13, 1.011), framealpha=1.0, frameon=False)
	forlegc = [str(forlegc[i])+'-pass' for i in range(len(forlegc))]
	leg = plt.legend(forlegp, forlegc, bbox_to_anchor=(0.0, 1.001, 1.0, 1.001), mode="expand", loc=3, borderaxespad=0.0, columnspacing=0.7, handletextpad=0.3, framealpha=1.0, frameon=False, ncol=4)
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels() + leg.get_texts()):
	    item.set_fontsize(fs)

	if outfile is not None:
		if ts is None: ts = t.strftime('%s')
		plt.savefig(outfile+"plotsetup-survey-" + ts +'.pdf', dpi=400, bbox_inches='tight') 
		plt.close()

	return fig, ax

def get_det_size(detectors, ndetx=NDETX, ndety=None):
	""" Get the area of the detector and the gap around it. """
	if ndety is None: ndety = ndetx
	# Assuming they move horizontally first
	dx, dy = (detectors[-1] - detectors[0])[2:4]
	return dx/(ndetx-1), dy/(ndety-1)

def _plot_pattern(origin, dx, dy, outv, ax=None, dithervec=None):
	""" Plot just one repetition of a pattern """

	lim_rec = list(origin) + [origin[0]+dx, origin[1]] + [origin[0]+dx, origin[1]+dy] + [origin[0], origin[1]+dy]
	outv = coverage.crop_survey(outv, lim_rec, verb=VERB, drop=True)

	n = outv.shape[0]
	polys = []
	for i,poly in zip(range(n), outv):
		p = pp.Poly(poly[0:2], 1.0)
		p.vertices(poly[2:])
		if p.weight>0 and p.area() > 1e-14:
			polys.append(p)
			#ax.plot(p.centre[0],p.centre[1],'r.')

	forlegp, forlegc = pp.add_Polys(ax, polys, maxa=0.99, lw=0.1, verts=False)
	forlegc = [str(forlegc[i])+'-pass' for i in range(len(forlegc))]

	x = [origin[0]]
	y = [origin[1]]
	if dithervec is not None:
		xdith = np.cumsum(dithervec[::2])
		ydith = np.cumsum(dithervec[1::2])
		#if np.sum(np.sqrt(xdith**2+ydith**2)) <= np.sqrt(dx**2+dy**2):
		x.extend(origin[0]+xdith)
		y.extend(origin[1]+ydith)
		ax.plot(x,y,'ko-')

	return forlegp, forlegc

def _create_survey(binary, outpath = './', pattern = 'J', size = 1.0, verb = 0, rewrite = False):
	""" Run the C-code and output the filename of the vertices created by Mangle. """

	inv = outpath + 'full-survey-' + pattern + '.vrt'
	outv = outpath + 'full-survey-overlaps-' + pattern + '.vrt'

	# Create the dither vector corresponding to the given pattern
	dithvec = size*coverage.PATTERNS[pattern]
	assert NDITH==int(len(dithvec)/2+1)

	# Run the create-euclid-patch.c ubercal code to construct the Mangle mask
	if not os.path.isfile(inv) or not os.path.isfile(outv) or rewrite: 
		runit = True
	else: 
		runit = False
		print "Warning, I'm using the old .vrt files in " + outpath + ' for my polygons!'
	if inv is None or outv is None or runit:
		inv2, outv2 = coverage.create_euclid_patch(dithvec, outpath, nra=NPOINT, ndec=NPOINT, ndetx=NDETX, binary=binary, verb=verb)
	if runit:
		os.rename(inv2, inv)
		os.rename(outv2, outv)
	else:
		inv2=inv
		outv2=outv

	# Read the Mangle in- & output file
	inv = coverage.read_vrt(inv, verb=verb, COORD_SCALE=1.0)
	outv = coverage.read_vrt(outv, verb=verb, COORD_SCALE=1.0)
	origin = (np.min(inv[:,2::2]), np.min(inv[:,3::2]))

	# Find patch of full coverage
	lim_rec = coverage.central_pointing(origin, COORD_SCALE=1.0)

	# Identify the edges that are close
	tmp = coverage.snap(np.vstack([lim_rec, outv[:,2:]]))
	lim_rec = tmp[0]; outv[:,2:] = tmp[1:]; tmp = None

	# Now throw away the polygons that not in the full coverage patch of the survey
	outv = coverage.crop_survey(outv, lim_rec, verb=VERB, drop=True)

	return inv, outv, lim_rec

if __name__=='__main__':

	BINARY = op.abspath(op.join(HERE,'../bin/create-euclid-patch'))
	OUTPATH = op.abspath(op.join(HERE,'../outputs/'))+'/'
	PATTERN = 'S'
	SIZE = 1.0
	VERB = 0
	PLOTPATH = OUTPATH
	PATTERNS = ['S', 'J', 'N', 'R', 'O', 'X']
	USEOLD = True

	# Record time stamp
	ts = t.strftime('%s')

	# Set up default plotting fonts
	fs = 15
	font = {'family' : 'serif',
    	    'weight' : 'normal',
        	'size'   : fs}
	rc('font', **font)	
	
	# PLOT 1:
	# Make the subplots of the different patterns
	fig, axes = plt.subplots(2,len(PATTERNS)/2, sharex=True, sharey=True, figsize=(10,10))
	plt.ticklabel_format(useOffset=False)

	# Plot each pattern
	i=0
	for ax, p in zip(np.hstack(axes), PATTERNS):
	
		# Get the vertices
		inv, outv, lim_rec = _create_survey(BINARY, OUTPATH, p, SIZE, VERB, not USEOLD)
		origin = lim_rec[:2]
		oned = inv[0:inv.shape[0]/NDITH]
		vector = SIZE*coverage.PATTERNS[p]
		if p == PATTERNS[0]:
			dx, dy = get_det_size(oned, ndetx=NDETX*NPOINT)

		# PRINT COVERAGE
		cov, tot = coverage.coverage(outv, lim_rec)
		print "For pattern " + p + ", size " + str(round(SIZE,1)) + ':'
		for c,f in zip(range(len(cov)),cov/tot*100):
			print "\t%d-pass coverage area fraction: %02.2f%%" % (c,f)

		# PLOT 1:
		if p == PATTERN:
			# Plot the full survey, highlight the full-coverage region
			fig2, axes2 = _plot_survey(inv, outv, oned, NRA=NPOINT, NDEC=NPOINT, outfile=PLOTPATH, ts=ts, scale=coverage.SCALE, fs=20)
		
		# SUBPLOT OF PLOT 2:
		if p == PATTERN:
			forlegp, forlegc = _plot_pattern(origin, dx, dy, outv, ax=ax, dithervec=vector/1.0)
		else:
			_plot_pattern(origin, dx, dy, outv, ax=ax, dithervec=vector/1.0)
		
		# Convert the pattern name if neccessary
		p = coverage.convpat(p)
		ax.text(origin[0]+dx*0.5, origin[1]+dy*0.475, p, fontsize=2*fs, family='serif', color='w', horizontalalignment='center')
		
		unit = "''"

		if i % axes.shape[1]==0: 
			ax.set_ylabel(r"$\delta$ [" + unit + "]")
			yticks = 1.0/4.0*np.arange(5)*dy
			ax.set_yticks(origin[1] + yticks)
			ax.set_yticklabels([str(round(t,1)) for t in yticks], font)
			#ax.yaxis.set_tick_params(labelsize=12)

		if i >= axes.shape[1]*axes.shape[0]-axes.shape[1]: 
			ax.set_xlabel(r"RA [" + unit + "]")
			xticks = 1.0/4.0*np.arange(5)*dx
			ax.set_xticks(origin[0] + xticks)
			ax.set_xticklabels([str(round(t,1)) for t in xticks], font)
			#ax.xaxis.set_tick_params(labelsize=12)

		ax.set_aspect('equal', adjustable='box-forced')
		ax.set_xlim([lim_rec[0], lim_rec[0]+dx])
		ax.set_ylim([lim_rec[1], lim_rec[1]+dy])
		
		i+=1
	
	#size = fig.get_size_inches()*fig.dpi
	fig.legend(forlegp, forlegc, fontsize=fs, bbox_to_anchor=(0.45, 0.47), loc="center", borderaxespad=0.0, framealpha=1.0, frameon=False, ncol=4)
	fig.set_size_inches(len(PATTERNS)/2 * 5.0, 9.0)
	#plt.tight_layout()

	if PLOTPATH is None:
		plt.show()
	else:
		plt.savefig(PLOTPATH+"plotsetup-patterns-" + ts +'.eps',dpi=400,bbox_inches='tight')
