#!/usr/bin/env python

# Import what we need to work with netcdf files
from netCDF4 import Dataset

# Used for arange
import numpy

# Copying files
from shutil import copy,copyfile

# Read user input
import getopt

# Constants stuff
import const

import os
import sys

############################################################
#
# Set up some config options

# End without a slash
const.upload_rel_dir="upload"

# Try to write cell bounding information
# I don't think this is properly implemented, as it seriously
# messes up the ability of ArcGIS to read in the NetCDF files
const.add_bound_info = False

# Add cell_measures data
const.add_cell_measures = True

# Add the Datum information.  Note, this seems to have no effect
const.set_datum = False

# Log file
const.do_log = True
const.log = "/tmp/pylog"

############################################################


def index(req):
	return "Hello"

def form():
	return """\
<html><body>
<form enctype="multipart/form-data" action="upload" method="post">
<p>File: <input type="file" name="file">:</p>
<p><label for="shift_cells">Shift Cells: <input type="checkbox" name="shift_cells" id="shift_cells" checked="checked" /></p>
<p><label for="debug">Debug: <input type="checkbox" name="debug" id="debug" /></p>
<p><input type="submit" value="Fix Spatial Information"></p>
</form>
</body></html>
"""

def upload(req):
	
	try: # Windows needs stdio set for binary mode.
		import msvcrt
		msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
		msvcrt.setmode (1, os.O_BINARY) # stdout = 1
	except ImportError:
		pass

	# A nested FieldStorage instance holds the file
	fileitem = req.form['file']

	# A nested FieldStorage instance holds the file
	shift_cells = True if 'shift_cells' in req.form else False

	# Output message
	debug = True if 'debug' in req.form else False
	

	# Test if the file was uploaded
	if fileitem.filename:

		# strip leading path from file name to avoid directory traversal attacks
		fname = os.path.basename(fileitem.filename)
		# build path to files directory
		dir_path = os.path.join(os.path.dirname(req.filename), const.upload_rel_dir)
		file_src = os.path.join(dir_path, fname)
		open(file_src, 'wb').write(fileitem.file.read())
		message = 'The file "%s" was uploaded successfully' % fname

	else:
		message = 'No file was uploaded'

	# Now that we have the file, lets pass it to our main function
	file_out=os.path.join(dir_path, "%s.%s"%('out', fname))
	rcode, stdout, stderr = fixIoapiSpatialInfo(file_src, file_out, shift_cells, debug)

	if const.do_log:
		f = open(const.log, 'a')
		f.write("Stdout: %s\n\nStderr: %s\n"%(stdout, stderr));
		f.close()
	
	if stderr == "" and not debug:
		req.headers_out.add('Content-Disposition', "attachment; filename=\"%s\""%os.path.basename(file_out))
		# Now, output the file
		f = open(file_out, 'r')
		return f.read();
	else:
		return """\
<html><body>
<p>%s</p>
<p><a href="./form">Upload another file</a></p>
<p>fname: %s</p>
<p>fileitem.filename: %s</p>
<p>stdout: %s</p>
<p>stderr: %s</p>
</body></html>
""" % (message, fname, fileitem.filename, stdout.replace("\n", "<br>"), stderr)


def getRectBounds(cols, rows, dx):

	# Implementing: http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/cf-conventions.html#cell-boundaries#{{{
	
	# used for arrays and loops, so adding one (because we don't have a zero)
	nc = len(cols)
	nr = len(rows)
	
	#colbnds = [x*2 for x in cols]
	colbnds=range(0, nr)
	
	# Initialize
	for j in range(0, nr):
		#print "j=%d"%j
		colbnds[j] = range(0, nc)
	
		for i in range(0, nc):
			#print "(j, i): (%d, %d)"%(j,i)
			colbnds[j][i] = range(0, 4)
	
	rowbnds=colbnds
	
	# Rows
	for j in range(0, nr):
		# Cols
		for i in range(0, nc):
			x = cols[i]
			y = rows[j]
	
			# Lower left
			colbnds[j][i][0] = x - dx/2
			rowbnds[j][i][0] = y - dx/2
	
			# Lower right
			colbnds[j][i][1] = x + dx/2
			rowbnds[j][i][1] = y - dx/2
	
			# Upper right
			colbnds[j][i][2] = x + dx/2
			rowbnds[j][i][2] = y + dx/2
	
			# Upper left
			colbnds[j][i][3] = x - dx/2
			rowbnds[j][i][3] = y + dx/2
	
	# Test the coordinates.  The lower right point of one box should be the
	# same as the lower left of the box to it's right
	#
	# Use this to test:
	# For 0 < j < n and 0 < i < m,
	#	If cells (j,i) and (j,i+1) are contiguous, then
	#		bnd(j,i,1)=bnd(j,i+1,0) 
	#		bnd(j,i,2)=bnd(j,i+1,3)
	#	If cells (j,i) and (j+1,i) are contiguous, then	
	#		bnd(j,i,3)=bnd(j+1,i,0) and bnd(j,i,2)=bnd(j+1,i,1)
	
	if True:
	
		# Rows
		for j in range(1, nr-1):
			# Cols
			for i in range(1, nc-1):
				print "(i,j)=(%d,%d), or [j][i]=[%d][%d]"%(i,j,j,i)
				if not colbnds[j][i][1] == colbnds[j][i+1][0]:
					print "Error! colbnds(%d,%d,%d) = %2.0f != colbnds(%d,%d,%d) = %2.0f"%(j,i,1,colbnds[j][i][1], \
						j,i+1,0,colbnds[j][i+1][0])
	
				if not colbnds[j][i][1] == colbnds[j][i+1][0]:
					print "Error! colbnds(%d,%d,%d) = %2.0f != colbnds(%d,%d,%d) = %2.0f"%(j,i,2,colbnds[j][i][2], \
						j,i+1,3,colbnds[j][i+1][3])
	
				if not colbnds[j][i][3] == colbnds[j+1][i][0]:
					print "Error! colbnds(%d,%d,%d) = %2.0f != colbnds(%d+1,%d,%d) = %2.0f"%(j,i,3,colbnds[j][i][3], \
						j,i,0,colbnds[j+1][i][0])
	
				if not colbnds[j][i][2] == colbnds[j+1][i][1]:
					print "Error! colbnds(%d,%d,%d) = %2.0f != colbnds(%d+1,%d,%d) = %2.0f"%(j,i,2,colbnds[j][i][2], \
						j,i,1,colbnds[j+1][i][1])
				print "--"#}}}
	
	return (colbnds,rowbnds)

def fixIoapiSpatialInfo(file_src, file_out, shift_cells = False, verbose = False):

	rcode=0#{{{
	stdout=""
	stderr=""

	if not os.access(file_src, os.R_OK):
		stderr += ("File %s is not readable!\n" % file_src)
		rcode=1
		return (rcode, stdout,stderr)
	
	if verbose:
		stdout += "Source: %s\n"%file_src
		stdout += "Output: %s\n"%file_out
		stdout += "Shifting Cells?: %s\n"%str(shift_cells)

	# Open the file
	## Check to see if the file exists
	nc = Dataset(file_src, 'r', format='NETCDF4')
	# Should include some sort of check here to make sure
	# the file is a proper Netcdf file, otherwise a bunch
	# of errors will be thrown below whenever we attempt
	# to access some properties on this.

	# Copy the source file to the output file, then we modify
	copyfile(file_src, file_out)
	nnc = Dataset(file_out, 'r+', clobber=True, format='NETCDF4')#}}}

	# Put the grid mapping on all the variables
	try:
		for var in nc.variables:
			varname = str(var)
			nnc.variables[varname].grid_mapping = 'Lambert_Conformal'
			if const.add_cell_measures:
				nnc.variables[varname].cell_measures = 'area: cell_area'
	except RuntimeError, err:
		stderr += ('Ack!  RuntimeError: %s\n' % str(err))

	try:
		# Lets create the projection variable
		llc=nnc.createVariable('Lambert_Conformal', 'i4')
		llc.grid_mapping_name="lambert_conformal_conic"
		llc.standard_parallel=[nc.P_ALP, nc.P_BET]
		llc.longitude_of_central_meridian=nc.P_GAM
		llc.latitude_of_projection_origin=nc.YCENT
		if const.set_datum:
			llc.datum='D_North_America_1983'
		else:
			if verbose:
				stdout += "Not attempting to set datum.\n"
	
		# I think the *ORIG variables aren't for the projection, but grids
		#llc.false_easting=-nc.XORIG/1000
		#llc.false_northing=-nc.YORIG/1000
		if verbose:
			stdout += "Lambert_Conformal variable added.\n"
	except (ValueError, RuntimeError) as (errno, strerror):
		stderr += "Ack!  Ran into a value error."
		stderr += "Variables in %s are: "%file_out
		stderr +=  "Current dimensions:"
		for dim in nnc.dimensions:
			stderr += "%s\n"%dim
		stderr +=  "\n"
		stderr += "Current variables:"
		for dim in nnc.variables:
			stderr += dim

	if const.add_cell_measures:
		try:
			# Lets create the projection variable
			cm=nnc.createVariable('cell_area', 'f4')
			cm.long_name = "area of grid cell";
			cm.standard_name = "area";
			cm.units = "km2";
			
			if verbose:
				stdout += "Cell Measures variable added.\n"
		except (ValueError, RuntimeError) as (errno, strerror):
			stderr += "Ack!  Ran into a value error."
			stderr += "Variables in %s are: "%file_out
			stderr +=  "Current dimensions:"
			for dim in nnc.dimensions:
				stderr += "%s\n"%dim
			stderr +=  "\n"
			stderr += "Current variables:"
			for dim in nnc.variables:
				stderr += dim
	else:
		if verbose:
			stdout += "Not adding cell measures attributes.\n";

	if const.add_bound_info:
		# Add number of vertices for cell geometry
		try:
			nnc.createDimension('nv', 4);
		except (ValueError, RuntimeError) as (errno, strerror):
			stderr += "Error: %s"%strerror
	else:
		if verbose:
			stdout += "Not adding boundary information.\n"

	# Add col/row variables
	try:
		col_var=nnc.createVariable('COL', 'f4', ('COL',))
		col_var.units = "km"
		col_var.long_name = "x coordinate of projection"
		col_var.standard_name = "projection_x_coordinate"

		if const.add_bound_info:
			col_var.bounds = "col_bnds"
			colbnds=nnc.createVariable('COL_bnds', 'f4', ('ROW','COL','nv'))

		if verbose:
			stdout += "Added COL variable\n"
	except (ValueError, RuntimeError) as (errno, strerror):
		stderr += "Ack!  Error"
		stderr += "Variables in %s are: "%file_out
		stderr +=  "Current dimensions:"
		for dim in nnc.dimensions:
			stderr += "%s\n"%dim
		stderr +=  "\n"
		stderr += "Current variables:"
		for dim in nnc.variables:
			stderr += dim
	
	
	try:
		row_var=nnc.createVariable('ROW', 'f4', ('ROW',))
		row_var.units = "km"
		row_var.long_name = "y coordinate of projection"
		row_var.standard_name = "projection_y_coordinate"

		if const.add_bound_info:
			row_var.bounds = "ROW_bnds"
			rowbnds=nnc.createVariable('ROW_bnds', 'f4', ('ROW','COL','nv'))

		if verbose:
			stdout += "Added ROW variable\n"
	except (ValueError, RuntimeError) as (errno, strerror):
		stderr += "Ack!  Ran into a value error."
		stderr += "Variables in %s are: "%file_out
		stderr +=  "Current dimensions:"
		for dim in nnc.dimensions:
			stderr += "%s\n"%dim
		stderr +=  "\n"
		stderr += "Current variables:"
		for dim in nnc.variables:
			stderr += dim
	
	# Populate row/col vars
	ncols = len(nc.dimensions['COL'])
	nrows = len(nc.dimensions['ROW'])
	
	yorig=nc.YORIG
	xorig=nc.XORIG
	cellsize=nc.XCELL
	if nc.XCELL != nc.YCELL:
		stderr += "Error! X/Y cell sizes are different!  Should throw an error.."
		return (rcode, stdout,stderr)
	
	cols = numpy.arange(xorig, xorig + (cellsize*ncols), cellsize)
	rows = numpy.arange(yorig, yorig + (cellsize*nrows), cellsize)


	if shift_cells:
		# ArcGIS interprets the starting location of a cell
		# differently than CMAQ.  Therefore, we shift these
		# cells.  (Double check which direction they should
		# be pushed.)
		cols += cellsize/2
		rows += cellsize/2
	
	# Convert to KM:
	cols = cols/1000.0
	rows = rows/1000.0

	# Some debug
	if verbose:
		stdout += "Cellsize = %d\n"%(cellsize/1000)
		stdout += "Cols(%d): %d:%d, note xorig=%d, ncols=%d\n"%(len(cols), cols[0], cols[len(cols) - 1], xorig/1000, ncols)
		stdout += "Rows(%d): %d:%d, note yorig=%d, nrows=%s\n"%(len(rows), rows[0], rows[len(rows) - 1], yorig/1000, nrows)
	
	try:
		col_var[:] = cols;
		row_var[:] = rows;
	except IndexError, err:
		stderr += "IndexError Exception: %s"%err

		stderr += "nc.COL = %d, len(cols) = "%(ncols,len(cols))
		stderr += "nc.ROW = %d, len(rows) = "%(nrows,len(rows))

	if const.add_bound_info:
		# Retrieve bounds
		cbvals, rbvals = getRectBounds(cols, rows, cellsize/1000)
		colbnds = cbvals
		rowbnds = rbvals
		if verbose:
			stdout += "Added bounds\n";
			stdout += str(colbnds)


	nnc.close();
	nc.close();

	return (rcode, stdout, stderr)

#!#############################################################
#!##
#!## Were we called from the command line?
#!##
#!#############################################################
#!#
#!#
#!## Define our usage help screen
#!#def usage():
#!#	print "This script is used to enter projection information into the Netcdf file using the IOAPI attributes in the file."
#!#	print "Usage:"
#!#	print "  -i <input file>"
#!#	print "  -o <output file>"
#!#	print "\nNote, input and output file cannot be the same."
#!#
#!#try:
#!#	opts, args = getopt.getopt(sys.argv[1:], 'i:o:hv', ['help'])
#!#except getopt.GetoptError, err:
#!#	verbose = False
#!#	file_src = "";
#!#	file_out = "";
#!#	for opt,val in opts:
#!#		if opt == "-i":
#!#			file_src = val
#!#		elif opt == "-o":
#!#			print "Error: ", str(err)
#!#			usage()
#!#			sys.exit(2)
#!#
#!#verbose = False
#!#file_src = "";
#!#file_out = "";
#!#for opt,val in opts:
#!#	if opt == "-i":
#!#		file_src = val
#!#	elif opt == "-o": 
#!#		file_out = val
#!#	elif opt in ("-h", "--help"):
#!#		usage();
#!#		sys.exit();
#!#	elif opt == "-v":
#!#		verbose = True
#!#	else:
#!#		print "uhh, do something..."
#!#
#!#if not os.access(file_src, os.R_OK):
#!#	sys.stderr.write("File %s is not readable!\n" % file_src)
#!#	usage()
#!#	sys.exit(2)
#!#
#!#if verbose:
#!#	print "Source:", file_src
#!#	print "Output:", file_out
#!#
#!#arr = fixIoapiSpatialInfo(file_src, file_out)
