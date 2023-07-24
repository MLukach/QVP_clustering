import sys
import os
import warnings
import numpy as np

from optparse import OptionParser
import copy
import pyart

import matplotlib.pyplot as plt

from matplotlib import rcParams


def read_options_and_input_variables(field_list,field_lables_list,argv):
	
	usage = """%prog [options] /input/directory/PPI_file.nc 
    Use the option -h or --help to get all possible options
    """

	# parsing options and arguments
	parser = OptionParser(usage=usage)
	parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
	parser.add_option("-l", "--list", dest="field_list", default=field_list, help="list of QVP fields to process")
	parser.add_option("-p", "--plot", action="store_true", dest="plot", default=False, help="plot option: default is False")
	parser.add_option("-e", "--elevation", dest="elev", default=3, help="elevation in degrees: default is 3 deg.")
	
	(options, files) = parser.parse_args()
	debug = options.debug
	verbose = options.verbose
	field_list = options.field_list
	elev=options.elev
	# if we what to generate plots on the run
	plot = options.plot
	
	# Geting input cariables
	if(len(argv)<3):
		print "The number of input values should be 3: the input PPI file name, clustering file and the output directory should be provided."
		sys.exit()
	
	ppi_file = argv[1]
	if verbose: print "First argument is ppi_file = {}".format(ppi_file)
	
	clustering_file = argv[2]
	if verbose: print "Second argument is {}".format(clustering_file)
	
	# with ./data plotds will be saved into ./data/20170517_QVP/plots/
	output_dir = argv[3]
	if verbose: print "Third argument is {}".format(output_dir)
	
	# create unique directory for this run
	directory = "{}/spectral_{}_".format(output_dir,"-".join(field_lables_list)) + "/"
	if not os.path.exists(directory):os.makedirs(directory)
	output_dir = directory
	
	
	return(debug,verbose,field_list,elev,plot,ppi_file,clustering_file,output_dir)



def read_ppi_data(file_):
	radar = pyart.io.read(file_)
	sweep_time = []
	dBZ_field = (copy.deepcopy(radar.fields['dBZ']['data']))
        
	return (radar,dBZ_field)


def field_fill_to_nan(radar, mask_field):
	
	
	if(np.ma.is_masked(radar.fields[mask_field]['data'])):
		mask = radar.fields[mask_field]['data'].mask
		radar.fields[mask_field]['data'] = radar.fields[mask_field]['data'].data
		radar.fields[mask_field]['data'][mask==True] = np.nan

	#radar[mask_field] = np.nan;
	return
	
def plot_ppi(display,sweep,output_dir,mask_field = 'dBZ',filename=None):
	# Create figure instance 
	fig = plt.figure()
	# Create axes in figure
	ax = fig.add_subplot(111)
	
	xlabel = 'Distance from Radar [km]'
	ylabel = 'Distance from Radar [km]'
	colorbar_label = 'Horizontal Reflectivity [dBZ]'
	
	#Define Colormap Lot of Options
	# We shoudl define a certain set
	plot_cmap='pyart_Theodore16'
	
	#sweep=2
	#Plot With Labels
	display.plot_ppi(mask_field, sweep=sweep, vmin=-15, vmax=15.,colorbar_orient='vertical',
	                 cmap=plot_cmap, axislabels=(xlabel, ylabel),colorbar_label=colorbar_label)
	
	#Add Grid to Plot
	display.plot_grid_lines(ax=None, col='k', ls=':')
	
	
	# Add Range Rings
	
	display.plot_range_ring(25, col='k', ls='--')
	display.plot_range_ring(50, col='k', ls='--')
	display.plot_range_ring(100., col='k', ls='-.')

	#Add Cross Hair
	#display.plot_cross_hair(5.)
	# Define paramters from package
	rcParams['axes.labelsize'] = 14
	rcParams['axes.titlesize'] = 20
	rcParams['xtick.labelsize'] = 12
	rcParams['ytick.labelsize'] = 11
	
	if(filename==None):filename = "test_ppi_sweep{}".format(sweep)
	savefile2 = output_dir + filename
	plt.savefig(savefile2 , format='png',dpi=200)
	plt.close()
	
	return (savefile2)
