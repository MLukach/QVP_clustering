# -*- coding: utf-8 -*-
# Plot cluster data 
# M. Lukach, May 2018
import sys
import os
from optparse import OptionParser

# from netCDF4 import Dataset
import warnings
import numpy as np
# import matplotlib.pyplot as plt
# import scipy as sp
import pandas as pd
#import netCDF4 as nc
#import h5py

from netCDF4 import Dataset
from netCDF4 import date2num
# import os
# import pyart
# import math
# import matplotlib.dates as mdates
import seaborn as sns
# import datetime as dt
import QVP_plot_functions as qvpp

# from sklearn.decomposition import PCA
import pytz
from astral import Astral
#from astral import LocationInfo
#from astral.sun import sun

#from astroplan import Observer
#import astropy.coordinates as coord
#import astropy.units as u
#import astropy.time as t
#import urllib

from clusters_lib import clusternode
from clusters_lib import read_clusters,plotting_tree,get_active_clusteres_data,get_all_children_full_name,get_night_hours_for_one_day,get_night_hours
import cluster_centroids_plot_lib as ccplib
import csv
import datetime as dt
import math


from scipy.stats import skew,kurtosis


import netCDF4 as nc

import matplotlib.dates as mdates

import matplotlib.pyplot as plt

from matplotlib import cm, colors

from math import pi
from math import radians
import copy

from pandas.plotting import parallel_coordinates

import xarray as xr
from collections import Counter

#from statistics import mode
 



plot_kwds = {'alpha' : 0.25, 's' : 80, 'linewidths':0}

max_height = 1000

def read_input_data_night(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case,elev,verbose=False,debug=True,sys_phi_list=[]):
	"""
	Reads input data from the qvp file and generates the panda dataset
	
	Parameters:
    -----------------------------------------
    qvp_file - file-object daataset from the read_file() function
    field_list - list of fields from the input file to be analysed
    upper_indx,count_thr,alts
    
    Returns:
    -----------------------------------------
    file_obj a dataset formed from the opend file
	"""
	print "in read_input_data_night"
	# get uPhiDP
	#f = qvp_file['uPhiDP/Means'][0:upper_indx]
	uPhiDP = np.reshape(np.where(qvp_file['uPhiDP/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['uPhiDP/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	#PhiDP = np.reshape(np.where(qvp_file['PhiDP/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['PhiDP/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	# get RhoHV for the non-meteo mask
	RhoHV = np.reshape(np.where(qvp_file['RhoHV/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['RhoHV/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
    
    #  Get height for range calculation to remove nearest ranges up to 400m
	Height = qvp_file['Height'][0:upper_indx]
	
    # Convert height to range assumint 20 deg elevation
	Range = np.cos(np.radians(int(elev))) * Height
	
	sys_phi_list.append(np.nanmedian(uPhiDP))		
	
	print field_list
	print "field_list[0:-3]"
	print field_list[0:-3]
	
	#for i,field in enumerate(field_list[0:-3]):
	for i,field in enumerate(field_list[0:-4]):
		if(debug):print "adding %d field (%s) to the matrix" %(i+1,field)
		field_count = field+'/Counts'
		field_mean = field+'/Means'
		
		f = qvp_file[field_mean][0:upper_indx]
		
		# remove bins near the radar up to 400 m range (due to influence by the side lobes) 
		#f[Range<=400,:] = np.nan
		
		if(debug):print "original shape is "
		if(debug):print f.shape
		X[:,i] = np.reshape(np.where(qvp_file[field_count[:]][0:upper_indx]>count_thr,f,np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
		#X[RhoHV < 0.8,i] = np.nan                
		#if (field is "uPhiDP"):
			#if(debug):
				#print "we are modifying the field"
				#print "np.nanmean(X[:,i]) is {}, i is {}, field is {}".format(np.nanmean(X[:,i]),i,field)
			#X[:,i] = X[:,i] - np.nanmedian(uPhiDP)
			#if(debug):
				#print "np.nanmean(PhiDP)"
				#print np.nanmedian(PhiDP)
				#print "np.nanmean(X[:,i])"
				#print np.nanmedian(X[:,i])
		if (field is "temperature_2"):
			X[:,i] = X[:,i] - 273.15
		if (field is "temperature"):
			X[:,i] = X[:,i] - 273.15
		if(debug):
			print "X[:,i].shape"
			print X[:,i].shape
			
		#if(debug):plot_hist(X,i,output_dir,field)
		
	
		if(i==0):
			# update Y time-vector
			X[:,-1] = np.reshape(np.where(qvp_file[field_count[:]][0:upper_indx]>count_thr,(qvp_file['Time'][:]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	
		
	
	# add Height field
	X[:,-3] = np.reshape(np.tile(alts[0:upper_indx], (len(qvp_file['Time'][:]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	
	
	# add Case field
	X[:,-2] = np.reshape(np.tile(case, (len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	if(debug):print "case index is np.unique(X[:,-2]) {}".format(np.unique(X[:,-2]))
	
	
	#get the night mask
	(hours_mask,night_hours_mask,dusk,dawn) = get_night_hours(1,np.asarray([[nc.num2date(x,qvp_file['Time'].units,qvp_file['Time'].calendar) for x in qvp_file['Time']]]),alts.shape)
	
	# add Night field, which is the 
	#X[:,-2] 
	night = np.reshape(hours_mask[0:upper_indx],(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	# if(debug):print ":::::::nightmask index is np.unique(night) {}".format(np.unique(night))
	# print ":::::::nightmask index is np.unique(night) {}".format(np.unique(night[night==False]))
	# print "np.where(hours_mask[0:upper_indx]==True)"
	# print np.where(hours_mask[0:upper_indx]==True)
	# print "np.where(night_hours_mask[0:upper_indx]==True)"
	# print np.where(night_hours_mask[0:upper_indx]==True)
	# print "night.shape"
	# print night.shape
	# print "X.shape"
	# print X.shape
	# print ":::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	# print "night_hours_mask.shape"
	# print night_hours_mask.shape
	# print "night_hours_mask[0:upper_indx].shape"
	# print night_hours_mask[0:upper_indx].shape
	X[night==False,0:-3] = np.nan
	
	
	
	
	indexes = ~np.isnan(X).any(axis=1)
	# update pd_X with the X data
	if(pd_X is None):pd_X = pd.DataFrame(data = X[indexes], columns = field_lables_list)
	else: 
		df_X = pd.DataFrame(data = X[indexes], columns = field_lables_list)
		pd_X = pd_X.append(df_X)

	if(debug):
		print "describe the data"
		print pd_X.describe()
		print pd_X.keys().tolist()
	return (pd_X,indexes,sys_phi_list,hours_mask)
	

def read_CIP(CIP_directory_name, CIP_filename, event = '20180214', CIP=15 ):
	print "In read_CIP and filename is:"
	#print "{}/{}_QVP/FAAM/{}".format(CIP_directory_name,event,CIP_filename)
	CIP_file = nc.Dataset("{}/{}_QVP/FAAM/{}".format(CIP_directory_name,event,CIP_filename),'r')
	conc_edge = None
	if('processed_group' in CIP_file.groups):#PICASSO_group
		processed_group = CIP_file.groups['processed_group']#['PICASSO_group']#['processed_group']#['PADS_group']#['raw_group']
		#print "!!!!! processed_group in place of PICASSO"#"PICASSO_group"	
		conc_ice = np.asarray(processed_group.variables['cip{}_iwc_psd'.format(CIP)])*1000#['cip{}_conc_psd_pads'.format(CIP)])#['cip{}_conc_psd'.format(CIP)])
		
		#np.asarray(processed_group.variables['cip{}_iwc_psd_ai'.format(CIP)])#['cip{}_conc_psd_pads'.format(CIP)])#['cip{}_conc_psd'.format(CIP)])
		conc_water = np.asarray(processed_group.variables['cip{}_lwc_psd'.format(CIP)])*1000#np.asarray(processed_group.variables['cip{}_lwc_psd_ai'.format(CIP)])
		if('PICASSO_group' in CIP_file.groups):#PICASSO_group
			conc_edge = np.asarray(processed_group.variables['cip{}_psd_edge'.format(CIP)])*1000#np.asarray(processed_group.variables['cip{}_lwc_psd_ai'.format(CIP)])
			conc_edge[conc_edge==-9999.0] = np.nan
			
	else:
		processed_group = CIP_file.groups['processed_group']#['PADS_group']#['raw_group']	
		#print "processed_group"	
		
		conc_ice = np.asarray(processed_group.variables['cip{}_iwc_psd'.format(CIP)])#['cip{}_conc_psd_pads'.format(CIP)])#['cip{}_conc_psd'.format(CIP)])
		conc_water = np.asarray(processed_group.variables['cip{}_lwc_psd'.format(CIP)])
		
	
	#print processed_group.variables
	
	#conc_ice = np.asarray(processed_group.variables['cip{}_iwc_psd_ai'.format(CIP)])#['cip{}_conc_psd_pads'.format(CIP)])#['cip{}_conc_psd'.format(CIP)])
	conc_ice[conc_ice==-9999.0] = np.nan
	#conc_water = np.asarray(processed_group.variables['cip{}_lwc_psd_ai'.format(CIP)])
	conc_water[conc_water==-9999.0] = np.nan
	
	
	PADS_group = CIP_file.groups['PADS_group']#['raw_group']	
	#print "PADS_group"	
	#print PADS_group.variables
	
	conc_pads = np.asarray(PADS_group.variables['cip{}_conc_psd_pads'.format(CIP)])#['cip{}_conc_psd'.format(CIP)])
	conc_pads[conc_pads==-9999.0] = np.nan
	
	bins = np.asarray(CIP_file.variables['cip{}_bin_centre'.format(CIP)])
	#print  "bins.shape"
	#print  bins.shape
	#print bins
	faam_times = CIP_file.variables['time']
	faam_altitude = CIP_file.variables['altitude']
	#print "len(faam_times)"
	##print len(faam_times)
	#print faam_times
	#jd = np.asarray(nc.num2date(faam_times[:],faam_times.units,faam_times.calendar), dtype=object)
	#print "jd"
	#print jd
		
			
	return 	(conc_ice,conc_water,conc_pads,conc_edge, bins, faam_times, faam_altitude)


def faam_time(zoom_interval,verbose):
	# transformes the zoom_interval into faam-related variables
	#Make aware datetime objects from the start and stop dates (with time zone info included in the object)
	utc = pytz.timezone("UTC")
	if(verbose):
		print "zoom_interval"
		print zoom_interval
		print len(zoom_interval)
		
	if((zoom_interval is not None)):
		if (len(zoom_interval)==2):
			if(verbose):print "zoom intterval is OK!!!"
			#Make aware datetime objects from the zoom_interval dates (with time zone info included in the object)
			zoom_start = dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
									int(zoom_interval[0][6:8]), int(zoom_interval[0][9:11]),
									int(zoom_interval[0][11:13]),int(zoom_interval[0][13:15]), tzinfo=utc)
			zoom_end = dt.datetime(int(zoom_interval[1][0:4]),int(zoom_interval[1][4:6]),
									int(zoom_interval[1][6:8]), int(zoom_interval[1][9:11]),
									int(zoom_interval[1][11:13]),int(zoom_interval[1][13:15]), tzinfo=utc)
			
			# Set the visibility intervals for each case
			# the visibility intervals were selected by hand from the kml files describing the FAAM flights
			if(zoom_interval[0][0:8] == '20170517'):
				# C013 flight
				faam_visib = [[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										06,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										16,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										36,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										48,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										59,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 14,
										14,00)]]
				cip15_file = "core-cloud-phy_faam_20170517_v002_r0_c013_cip15.nc"
				cip100_file = "core-cloud-phy_faam_20170517_v002_r0_c013_cip100.nc"
				#"{}/{}_QVP/FAAM/{}".format(CIP_directory_name,event,CIP_filename)
			elif(zoom_interval[0][0:8] == '20180124'):
				# C076 flight
				faam_visib = [[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 10,
										25,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 10,
										40,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 11,
										00,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 11,
										12,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 11,
										32,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 11,
										50,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 12,
										07,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 12,
										26,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 12,
										48,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										03,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										24,35),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 13,
										39,00)]]
				
				cip15_file = "core-cloud-phy_faam_20180124_v005_r0_c076_cip15.nc"
				cip100_file = "core-cloud-phy_faam_20180124_v005_r0_c076_cip100.nc"
				#cip15_file = "core-cloud-phy_faam_20180124_v002_r0_c076_cip15.nc"
				#cip100_file = "core-cloud-phy_faam_20180124_v002_r0_c076_cip100.nc"
				
			elif(zoom_interval[0][0:8] == '20180213'):
				# C081 flight
				faam_visib = [[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 5,
										38,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 5,
										48,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 6,
										18,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 6,
										28,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 7,
										02,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 7,
										14,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 7,
										41,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 7,
										55,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 8,
										13,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 8,
										25,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 8,
										42,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 8,
										54,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 9,
										06,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 9,
										17,00)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 9,
										29,00),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 9,
										37,00)]]
				cip15_file = "core-cloud-phy_faam_20180213_v005_r0_c081_cip15.nc"#"core-cloud-phy_faam_20180213_v002_r0_c081_cip15.nc"
				cip100_file = "core-cloud-phy_faam_20180213_v005_r0_c081_cip100.nc"#"core-cloud-phy_faam_20180213_v002_r0_c081_cip100.nc"
				#cip15_file = "core-cloud-phy_faam_20180213_v002_r0_c081_cip15.nc"#"core-cloud-phy_faam_20180213_v002_r0_c081_cip15.nc"
				#cip100_file = "core-cloud-phy_faam_20180213_v002_r0_c081_cip100.nc"#"core-cloud-phy_faam_20180213_v002_r0_c081_cip100.nc"
			elif(zoom_interval[0][0:8] == '20180214'):
				# C082 flight to be updated by Manchester team
 				faam_visib = [[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 16,
										47,54),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 16,
										56,14)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 17,
										14,34),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 17,
										26,14)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 17,
										44,34),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 17,
										56,14)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 18,
										16,14),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 18,
										28,44)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 18,
										50,24),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 19,
										03,44)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 19,
										16,14),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 19,
										28,44)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 19,
										47,54),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 20,
										01,14)],[dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 20,
										19,34),dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
										int(zoom_interval[0][6:8]), 20,
										24,34)]]
				cip15_file = "core-cloud-phy_faam_20180214_v005_r0_c082_cip15.nc"#"core-cloud-phy_faam_20180214_v002_r0_c082_cip15.nc"
				cip100_file = "core-cloud-phy_faam_20180214_v005_r0_c082_cip100.nc"#"core-cloud-phy_faam_20180214_v002_r0_c082_cip100.nc"
				#cip15_file = "core-cloud-phy_faam_20180214_v002_r0_c082_cip15.nc"#"core-cloud-phy_faam_20180214_v002_r0_c082_cip15.nc"
				#cip100_file = "core-cloud-phy_faam_20180214_v002_r0_c082_cip100.nc"#"core-cloud-phy_faam_20180214_v002_r0_c082_cip100.nc"
		
		else:
			if verbose: print "Zoom interval should be a list with two dates. The event's start and end are used instead."
	
	return(zoom_start,zoom_end,faam_visib,cip15_file,cip100_file)


def read_file(directory,event,elev,verbose=False,debug=False):
	"""
	Read file from the provided directory and send back the file-object.
	
	Parameters:
    -----------------------------------------
    directory a string with the directory name
    event a string with the date of the event in YYYYMMDD format (20170517)
    
    Returns:
    -----------------------------------------
    file_obj a dataset formed from the opend file
	"""
	
	file_name = "%s/%s_QVP/%s_QVP_%sdeg.nc" %(directory,event,event,elev)
	
	
	#Dataset('./data/20170517_QVP/20170517_QVP_20deg.nc','r')#../scripts/QVP_20160706_20deg.nc','r')
	#./data/20170517_QVP/20170517_QVP_90deg.nc
			
	file_obj = nc.Dataset(file_name,'r')#../scripts/QVP_20160706_20deg.nc')
	if verbose or debug: 
		print file_name
		# print "file_obj.dimensions.keys()"
		# print file_obj.dimensions.keys()
		# print file_obj.variables.keys()
	return(file_obj)


def trap_time(date, zoom_interval,verbose):
	""""
	Transforms the zoom_interval into datetime objects using the date. 
	Makes aware datetime objects from the start and stop dates (with time zone info included in the object)
	
	Parameters:
    -----------------------------------------
    zoom_interval is a vector with two strings: ['20170517T090000','20170517T145959']
    
    Returns:
    -----------------------------------------
    zoom_start,zoom_end : two datetime objects
	"""
	utc = pytz.timezone("UTC")
	if(verbose):
		print "in trap_time() function"
		print "trap_interval is {}".format(zoom_interval)
		
	date = dt.datetime(int(date[0:4]),int(date[4:6]),
									int(date[6:8]), 0,0,0)
	date1 = date + datetime.timedelta(days=1)
	date1 = date1.strftime("%Y%m%d")
									
	if((zoom_interval is not None)):
		if (len(zoom_interval)==2):
			#Make aware datetime objects from the zoom_interval dates (with time zone info included in the object)
			zoom_start = dt.datetime(int(date[0:4]),int(date[4:6]),
									int(date[6:8]), int(zoom_interval[0][9:11]),
									int(zoom_interval[0][11:13]),int(zoom_interval[0][13:15]), tzinfo=utc)
			zoom_end = dt.datetime(int(date1[0:4]),int(date1[4:6]),
									int(date1[6:8]), int(zoom_interval[1][9:11]),
									int(zoom_interval[1][11:13]),int(zoom_interval[1][13:15]), tzinfo=utc)
			
		else:
			if verbose: print "Zoom interval should be a list with two dates. The event's start and end are used instead."
	
	return(zoom_start,zoom_end)

def event_time(zoom_interval,verbose):
	""""
	Transforms the zoom_interval into datetime objects. 
	Makes aware datetime objects from the start and stop dates (with time zone info included in the object)
	
	Parameters:
    -----------------------------------------
    zoom_interval is a vector with two strings: ['20170517T090000','20170517T145959']
    
    Returns:
    -----------------------------------------
    zoom_start,zoom_end : two datetime objects
	"""
	utc = pytz.timezone("UTC")
	if(verbose):
		print "in event_time() function"
		print "zoom_interval is {}".format(zoom_interval)
	if((zoom_interval is not None)):
		if (len(zoom_interval)==2):
			#Make aware datetime objects from the zoom_interval dates (with time zone info included in the object)
			zoom_start = dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
									int(zoom_interval[0][6:8]), int(zoom_interval[0][9:11]),
									int(zoom_interval[0][11:13]),int(zoom_interval[0][13:15]), tzinfo=utc)
			zoom_end = dt.datetime(int(zoom_interval[1][0:4]),int(zoom_interval[1][4:6]),
									int(zoom_interval[1][6:8]), int(zoom_interval[1][9:11]),
									int(zoom_interval[1][11:13]),int(zoom_interval[1][13:15]), tzinfo=utc)
			
		else:
			if verbose: print "Zoom interval should be a list with two dates. The event's start and end are used instead."
	
	return(zoom_start,zoom_end)
	
	
def plot_cluster_field_with_faam(clusters_field,
               Cmap,
               Cnorm,
               x_lims,
               vmin, vmax,new_labels,timeticks,savefile,alts=None,
               count_threshold=None, ax1=None, hs=None, hvis=None,timeserie=None):
				   
    print "!!!!     ---------------        IN plot_cluster_field_with_faam"
    print savefile
    
    hmax = max_height
    date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    
    #palette = sns.color_palette('deep', len(new_labels)+1)
    #print "len(timeticks)"
    #print len(timeticks)
    
    #print "clusters_field.shape"
    #print clusters_field.shape
    
    #bounds=[0,1,2,3]
    bounds=np.linspace(0,len(new_labels),len(new_labels)+1)
    #print "bounds"
    #print bounds
    
    norm = colors.BoundaryNorm(bounds, Cmap.N)
    #print "!!!!!   Cmap.N"
    #print Cmap.N
    
    #print "vmin"
    #print vmin
    
    vmax = vmax +1
    #print "vmax"
    #print vmax
    
    #print "norm"
    #print norm
    
    #print "!!!!!   np.unique(clusters_field[clusters_field != np.nan])"
    #print np.unique(clusters_field[clusters_field <> np.nan])
    #print int(np.nanmax(np.unique(clusters_field)))
    
    
    fig, ax1 = plt.subplots(sharey=False,figsize=(10*(len(timeticks)-1)/4,5))
    
		
    this_case_clusters = []
    this_case_ts = []
    #this_cases_clusternames = []
    
    if( ax1 is None):
		#print "ax is None !!!"
		
		if(hs is not None):ax1 = hs.plot(color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
		if(hvis is not None):hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    
		#print "clusters_field.shape before printing"
		#print clusters_field.shape
		#print "x_lims[0], x_lims[1],alts[0],alts[-1]"
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    else:
		
		if(hs is not None):hs.plot(ax=ax1,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
		if(hvis is not None):
			hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
			timeserie = np.array([each_ts.replace(microsecond=0) for each_ts in timeserie])
			#print "timeserie.shape"
			#print timeserie.shape
			timeserie = pd.Series(np.arange(len(timeserie)),timeserie)#pd.DataFrame(timeserie)
			
			events = np.split(hvis, np.where(np.isnan(hvis))[0])
			
			 
			
			# removing NaN entries
			events = [ev[~np.isnan(ev)] for ev in events if not isinstance(ev, np.ndarray)]
			# removing empty DataFrames
			events = [ev for ev in events if not ev.empty]
			
			for each_event in events:
				#print "each_event"
				#print each_event
				if(each_event.empty == False):print "		event is {} it starts {} and ends with {} the height is {}".format(each_event,min(each_event.keys()),max(each_event.keys()),each_event[min(each_event.keys())])
				for i,dt in enumerate(each_event.keys()):
					#print "dt is {}".format(pd.to_datetime(dt))
					dt = pd.to_datetime(dt)
					
					t = timeserie.index.get_loc(dt, method='nearest')#timeserie.iloc[timeserie.index.get_loc(dt, method='nearest')]
					
					
					#print "t is {}".format(t)
					# print " timeserie[t] = {}".format(timeserie.index[t])
					# print "each_event.keys() = {}".format(each_event.keys())
					#print "each_event.keys()[{}] = {}".format(i,each_event.keys()[i])
					#print "each_event[each_event.keys()[{}]] = {}".format(i,each_event[each_event.keys()[i]])
					# print "alts = {}".format(alts)
					# print "np.abs(alts - each_event[min(each_event.keys())]) = {}".format(np.abs(alts - each_event[each_event.keys()[i]]))
					# print "np.abs(alts - each_event[min(each_event.keys())]).argmin() = {}".format(np.abs(alts - each_event[each_event.keys()[i]]).argmin())
					# print "np.abs(alts - each_event[each_event.keys()[i]])[np.abs(alts - each_event[each_event.keys()[i]]).argmin()] = {}".format(np.abs(alts - each_event[each_event.keys()[i]])[np.abs(alts - each_event[each_event.keys()[i]]).argmin()])
					
					h = np.abs(alts - each_event[each_event.keys()[i]]).argmin()
					
					# print "clusters_field[:,t] = {}".format(clusters_field[:,t])
					# print "clusters_field[h,t] = {}".format(clusters_field[h,t])
					
					
					
					this_case_clusters.append(clusters_field[h,t])
					this_case_ts.append(dt)
					#if(np.isnan(clusters_field[h,t])):this_cases_clusternames.append("--*--")
					#else:this_cases_clusternames.append(new_labels[int(clusters_field[h,t])])
					#print this_cases_clusternames
	
	    
		#im = ax1.imshow(clusters_field,origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap, norm=Cnorm)
		#print "clusters_field.shape before printing"
		#print clusters_field.shape
		#print "alts[0],alts[-1]"
		#print alts[0],alts[-1]
		#faam_times.units,faam_times.calendar) 
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
		
		
    qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax1, im=im)
    #qvpp.t_contourplt(T[0:clusters_field.shape[0]], x_lims, alts[0], alts[-1], ax=ax1)
    
    
    
    plt.tight_layout()
    plt.savefig(savefile , format='png',dpi=100)
    
    #this_case_clusters = np.array(this_case_clusters).astype(int)+1
    
    
    if (len(this_case_clusters) > 0):
	    fig, ax2 = plt.subplots(sharey=False,figsize=(10,2))
	    a = np.array(this_case_clusters).reshape(1, -1)#np.linspace(0, 1, 256).reshape(1, -1)
	    #print ('Cluters observed by FAAM instruments on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    
	    #print "this_case_clusters = {}".format(this_case_clusters)
	    #print a
	    
	   # df_FAAM = pd.DataFrame(data={"data":np.array(this_case_clusters), "clusters":np.array(this_cases_clusternames)},index=np.array(this_case_ts))
	    ax2.imshow(a, aspect='auto', cmap=Cmap, origin='lower',vmin=vmin, vmax=vmax)
	    start, end = ax2.get_xlim()
	    stepsize = np.abs(end - start)/len(new_labels)
	    
	    for x,s in zip(np.arange(start+stepsize/4,end+stepsize,stepsize),new_labels):
			
			plt.text(x, 0, s, fontsize=12)
	    
	    
	    ax2.xaxis.set_ticks(np.arange(start, end, stepsize))
	    ax2.grid(color='k', linestyle='-', linewidth=1, axis ='x')
	    ts_time = [el.strftime('%H:%M:%S') for el in this_case_ts]
	    ax2.xaxis.set_ticklabels(ts_time)
	    #pos = list(ax2.get_position().bounds)
	    
	    plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30 )
	    ax2.set_yticks([])
	    plt.title('Cluters observed by FAAM instruments on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    plt.tight_layout()
	    
	    plt.savefig(savefile[:-4]+"_FAAM.png" , format='png',dpi=100)
	
    plt.close(savefile)
    
    
    print "figure %s is ok" %savefile
	
    return im



def plot_cluster_field_night_hours(clusters_field,
               Cmap,
               Cnorm,
               x_lims,
               vmin, vmax,new_labels,timeticks,savefile,alts=None,
               count_threshold=None, ax1=None, hs=None, hvis=None,timeserie=None):
				   
    print "!!!!     ---------------        IN plot_cluster_field_with_faam"
    print savefile
    
    hmax = max_height
    date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    
    #palette = sns.color_palette('deep', len(new_labels)+1)
    print "len(timeticks)"
    print len(timeticks)
    
    print "clusters_field.shape"
    print clusters_field.shape
    
    #bounds=[0,1,2,3]
    bounds=np.linspace(0,len(new_labels),len(new_labels)+1)
    #print "bounds"
    #print bounds
    
    norm = colors.BoundaryNorm(bounds, Cmap.N)
    #print "!!!!!   Cmap.N"
    #print Cmap.N
    
    print "vmin"
    print vmin
    
    vmax = vmax +1
    print "vmax"
    print vmax
    
    #print "norm"
    #print norm
    
    #print "!!!!!   np.unique(clusters_field[clusters_field != np.nan])"
    #print np.unique(clusters_field[clusters_field <> np.nan])
    #print int(np.nanmax(np.unique(clusters_field)))
    
    
    fig, ax1 = plt.subplots(sharey=False,figsize=(10*(len(timeticks)-1)/4,5))
    
		
    this_case_clusters = []
    this_case_ts = []
    #this_cases_clusternames = []
    
    if( ax1 is None):
		print "ax is None !!!"
		
		if(hs is not None):ax1 = hs.plot(color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
		if(hvis is not None):hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    
		print "clusters_field.shape before printing"
		print clusters_field.shape
		print "x_lims[0], x_lims[1],alts[0],alts[-1]"
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    else:
		
		if(hs is not None):hs.plot(ax=ax1,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
		if(hvis is not None):
			hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
			timeserie = np.array([each_ts.replace(microsecond=0) for each_ts in timeserie])
			print "timeserie.shape"
			print timeserie.shape
			timeserie = pd.Series(np.arange(len(timeserie)),timeserie)#pd.DataFrame(timeserie)
			
			events = np.split(hvis, np.where(np.isnan(hvis))[0])
			
			 
			
			# removing NaN entries
			events = [ev[~np.isnan(ev)] for ev in events if not isinstance(ev, np.ndarray)]
			# removing empty DataFrames
			events = [ev for ev in events if not ev.empty]
			
			for each_event in events:
				#print "each_event"
				#print each_event
				if(each_event.empty == False):print "		event is {} it starts {} and ends with {} the height is {}".format(each_event,min(each_event.keys()),max(each_event.keys()),each_event[min(each_event.keys())])
				for i,dt in enumerate(each_event.keys()):
					#print "dt is {}".format(pd.to_datetime(dt))
					dt = pd.to_datetime(dt)
					
					t = timeserie.index.get_loc(dt, method='nearest')#timeserie.iloc[timeserie.index.get_loc(dt, method='nearest')]
					
					
					#print "t is {}".format(t)
					# print " timeserie[t] = {}".format(timeserie.index[t])
					# print "each_event.keys() = {}".format(each_event.keys())
					#print "each_event.keys()[{}] = {}".format(i,each_event.keys()[i])
					#print "each_event[each_event.keys()[{}]] = {}".format(i,each_event[each_event.keys()[i]])
					# print "alts = {}".format(alts)
					# print "np.abs(alts - each_event[min(each_event.keys())]) = {}".format(np.abs(alts - each_event[each_event.keys()[i]]))
					# print "np.abs(alts - each_event[min(each_event.keys())]).argmin() = {}".format(np.abs(alts - each_event[each_event.keys()[i]]).argmin())
					# print "np.abs(alts - each_event[each_event.keys()[i]])[np.abs(alts - each_event[each_event.keys()[i]]).argmin()] = {}".format(np.abs(alts - each_event[each_event.keys()[i]])[np.abs(alts - each_event[each_event.keys()[i]]).argmin()])
					
					h = np.abs(alts - each_event[each_event.keys()[i]]).argmin()
					
					# print "clusters_field[:,t] = {}".format(clusters_field[:,t])
					# print "clusters_field[h,t] = {}".format(clusters_field[h,t])
					
					
					
					this_case_clusters.append(clusters_field[h,t])
					this_case_ts.append(dt)
					#if(np.isnan(clusters_field[h,t])):this_cases_clusternames.append("--*--")
					#else:this_cases_clusternames.append(new_labels[int(clusters_field[h,t])])
					#print this_cases_clusternames
	
	    
		#im = ax1.imshow(clusters_field,origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap, norm=Cnorm)
		print "clusters_field.shape before printing"
		print clusters_field.shape
		print "alts[0],alts[-1]"
		print alts[0],alts[-1]
		#faam_times.units,faam_times.calendar) 
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
		
		
    qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax1, im=im)
    #qvpp.t_contourplt(T[0:clusters_field.shape[0]], x_lims, alts[0], alts[-1], ax=ax1)
    
    
    
    plt.tight_layout()
    plt.savefig(savefile , format='png',dpi=100)
    
    #this_case_clusters = np.array(this_case_clusters).astype(int)+1
    
    
    if (len(this_case_clusters) > 0):
	    fig, ax2 = plt.subplots(sharey=False,figsize=(10,2))
	    a = np.array(this_case_clusters).reshape(1, -1)#np.linspace(0, 1, 256).reshape(1, -1)
	    #print ('Cluters observed by FAAM instruments on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    
	    #print "this_case_clusters = {}".format(this_case_clusters)
	    #print a
	    
	   # df_FAAM = pd.DataFrame(data={"data":np.array(this_case_clusters), "clusters":np.array(this_cases_clusternames)},index=np.array(this_case_ts))
	    ax2.imshow(a, aspect='auto', cmap=Cmap, origin='lower',vmin=vmin, vmax=vmax)
	    start, end = ax2.get_xlim()
	    stepsize = np.abs(end - start)/len(new_labels)
	    
	    for x,s in zip(np.arange(start+stepsize/4,end+stepsize,stepsize),new_labels):
			
			plt.text(x, 0, s, fontsize=12)
	    
	    
	    ax2.xaxis.set_ticks(np.arange(start, end, stepsize))
	    ax2.grid(color='k', linestyle='-', linewidth=1, axis ='x')
	    ts_time = [el.strftime('%H:%M:%S') for el in this_case_ts]
	    ax2.xaxis.set_ticklabels(ts_time)
	    #pos = list(ax2.get_position().bounds)
	    
	    plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30 )
	    ax2.set_yticks([])
	    plt.title('Cluters observed by FAAM instruments on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    plt.tight_layout()
	    
	    plt.savefig(savefile[:-4]+"_FAAM.png" , format='png',dpi=100)
	
    plt.close(savefile)
    
    
    print "figure %s is ok" %savefile
	
    return im


#--------------
def plot_faam_with_centroids(clusters_field,
               Cmap,
               Cnorm,
               x_lims,
               vmin, vmax,new_labels,timeticks,savefile,alts=None,
               count_threshold=None, ax1=None, hs=None, hvis=None,timeserie=None, faam_h=None,faam_times=None):
#plot_faam_with_centroids(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,savefile_faam,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1],faam_h=h_list[p-1],faam_times=faam_times_list[p-1])
		
				   
    print "!!!!     ---------------        IN plot_faam_with_centroids"
    
    print "clusters_field"
    print clusters_field.shape
    
    #print "temperature_field"
    #print temperature_field.shape
    
    # if (faam_h is not None):
		# #print len(faam_h)
		# #print "faam_h.keys()"
		# print faam_h.keys()
		
		# print "faam_h[:][0]"
		# print faam_h[:][0]
		
		# print "faam_h[0][0]"
		# print faam_h[0][0]
		
		
		
		# print "faam_h[0][len(faam_h)-4]"
		# print faam_h[0][len(faam_h)-4]
		# print "faam_h[0]len(faam_h)-1]"
		# print faam_h[0][len(faam_h)-1]
    
		# print "faam_times"
		# print faam_times
		# print len(faam_times)
   	
    
    # transform faam dates to integers
    #jd = np.asarray(nc.num2date(faam_times[:],faam_times.units,faam_times.calendar), dtype=object)
    ##print "jd"
    #print jd
    #faam_height = pd.DataFrame(h[:],index=jd)
    #print "faam_height"
    #print faam_height
    #faam_dates_in_times = np.asarray(nc.date2num(jd,times.units,calendar=times.calendar)).astype(int)
    #print "faam_dates_in_times"
    #print faam_dates_in_times
				
    
			
    
    hmax = max_height
    date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    bounds=np.linspace(0,len(new_labels),len(new_labels)+1)
    norm = colors.BoundaryNorm(bounds, Cmap.N)
    vmax = vmax +1
   	
    this_case_clusters = []
    this_case_ts = []
    #this_cases_clusternames = []
    
   # if(hs is not None):hs.plot(ax=ax1,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
    if(hvis is not None):
		
		pd_faam_times = pd.Series(np.arange(len(faam_times)),faam_times)#pd.DataFrame(timeserie)
	
		
		#hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
		timeserie = np.array([each_ts.replace(microsecond=0) for each_ts in timeserie])
		#print "timeserie"
		#print timeserie
		timeserie = pd.Series(np.arange(len(timeserie)),timeserie)#pd.DataFrame(timeserie)
		
		events = np.split(hvis, np.where(np.isnan(hvis))[0])
		# removing NaN entries
		events = [ev[~np.isnan(ev)] for ev in events if not isinstance(ev, np.ndarray)]
		# removing empty DataFrames
		events = [ev for ev in events if not ev.empty]
		
		#print "events"
		#print events
		#print len(events)
		
		
		for each_event in events:
			print "each_event"
			print each_event
			print "each_event.keys()"
			print each_event.keys()
			print "enumerate(each_event.keys() = {} ".format( enumerate(each_event.keys()))
			print "each_event.empty "
			print each_event.empty 
			print "		event is {} it starts {} and ends with {} the height is {}".format(each_event,min(each_event.keys()),max(each_event.keys()),each_event[min(each_event.keys())])
			
			
			#if(each_event.empty == False):print "		event is {} it starts {} and ends with {} the height is {}".format(each_event,min(each_event.keys()),max(each_event.keys()),each_event[min(each_event.keys())])
			
			for i,dt in enumerate(each_event.keys()):
				print "i={}".format(i)
				print "dt is {}".format(pd.to_datetime(dt))
				dt = pd.to_datetime(dt)
				t = timeserie.index.get_loc(dt, method='nearest')#timeserie.iloc[timeserie.index.get_loc(dt, method='nearest')]
				
				
				t_faam = pd_faam_times.index.get_loc(dt, method='nearest')
				#print "t_faam = {}".format(t_faam)
				#print pd_faam_times[t_faam]
				#print "pd_faam_times[t_faam-10] = {}".format(pd_faam_times[t_faam-10])
				#print "pd_faam_times[t_faam+10] = {}".format(pd_faam_times[t_faam+10])
				#print faam_times[t_faam-10:t_faam+10]
				#print "faam_h[t_faam-10:t_faam+10]"
				#print faam_h[t_faam-10:t_faam+10]
				
				
				
				
				#print "each_event[each_event.keys()[i]]"
				#print each_event[each_event.keys()[i]]
				
				h = np.abs(alts - each_event[each_event.keys()[i]]).argmin()
				#print "np.abs(alts - each_event[each_event.keys()[i]])"
				#print np.abs(alts - each_event[each_event.keys()[i]])
				
				print "h,t"
				print h,t
				
				print "temperature_field[h,t]"
				print temperature_field[h,t]
				
				
				print "clusters_field[h,t]"
				print clusters_field[h,t]
				
				#print "dt"
				#print dt
				
				this_case_clusters.append(clusters_field[h,t])
				this_case_ts.append(dt)
				#if(np.isnan(clusters_field[h,t])):this_cases_clusternames.append("--*--")
				#else:this_cases_clusternames.append(new_labels[int(clusters_field[h,t])])
				

    
    if (len(this_case_clusters) > 0):
		
	    fig, ax2 = plt.subplots(sharey=False,figsize=(10,2))
	    a = np.array(this_case_clusters).reshape(1, -1)#np.linspace(0, 1, 256).reshape(1, -1)
	    print "a"
	    print a
	    
	   # df_FAAM = pd.DataFrame(data={"data":np.array(this_case_clusters), "clusters":np.array(this_cases_clusternames)},index=np.array(this_case_ts))
	    ax2.imshow(a, aspect='auto', cmap=Cmap, origin='lower',vmin=vmin, vmax=vmax, norm=Cnorm,interpolation='none')
	    
	    #im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap, alpha=alpha)#, norm=Cnorm)
		
	    start, end = ax2.get_xlim()
	    stepsize = np.abs(end - start)/len(new_labels)
	    
	    for x,s in zip(np.arange(start+stepsize/4,end+stepsize,stepsize),new_labels):
			
			plt.text(x, 0, s, fontsize=12)
	    
	    
	    ax2.xaxis.set_ticks(np.arange(start, end, stepsize))
	    ax2.grid(color='k', linestyle='-', linewidth=1, axis ='x')
	    ts_time = [el.strftime('%H:%M:%S') for el in this_case_ts]
	    ax2.xaxis.set_ticklabels(ts_time)
	    #pos = list(ax2.get_position().bounds)
	    
	    plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30 )
	    ax2.set_yticks([])
	    plt.title('Cluters corresponding to the FAAM altitudes on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    plt.tight_layout()
	    
	    plt.savefig(savefile[:-4]+"_FAAM.png" , format='png',dpi=300)
	
    plt.close(savefile)
    
    
    print "figure %s is ok" %savefile
	
    #return im

#--------------


def plot_cluster_field_with_faam_and_centroids(clusters_field,
               Cmap,
               Cnorm,
               x_lims,
               vmin, vmax,new_labels,timeticks,savefile,alts=None,
               count_threshold=None, ax1=None, hs=None, hvis=None,timeserie=None,alpha=0.5):
				   
    print "!!!!     ---------------        IN plot_cluster_field_with_faam_and_centroids"
    
    hmax = max_height
    date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
 
    bounds=np.linspace(0,len(new_labels),len(new_labels)+1)
    norm = colors.BoundaryNorm(bounds, Cmap.N)
    vmax = vmax + 1
    
    
    fig, ax1 = plt.subplots(sharey=False,figsize=(10,5))
    
		
    this_case_clusters = []
    this_case_ts = []
    #this_cases_clusternames = []
    
    if( ax1 is None):
		#print "ax is None !!!"
		
		
		
		#print "hs"
		#print hs
		#print len(hs)
		#print "hvis"
		#print hvis
		#print len(hvis)
		
		if(hs is not None):ax1 = hs.plot(color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
		if(hvis is not None):hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    else:
		
		
		
		
		if(hs is not None):
			hs.plot(ax=ax1,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True) # alpha=0.5,
			#print "hs"
			#print hs
			#print len(hs)
		if(hvis is not None):
			#print "hvis"
			#print hvis
			#print len(hvis)
			hvis.plot(ax=ax1,color = 'black', alpha=0.5, style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
			timeserie = np.array([each_ts.replace(microsecond=0) for each_ts in timeserie])
			timeserie = pd.Series(np.arange(len(timeserie)),timeserie)#pd.DataFrame(timeserie)
			
			events = np.split(hvis, np.where(np.isnan(hvis))[0])
			# removing NaN entries
			events = [ev[~np.isnan(ev)] for ev in events if not isinstance(ev, np.ndarray)]
			# removing empty DataFrames
			events = [ev for ev in events if not ev.empty]
			#print "events are {}".format(events)
			for each_event in events:
				if(each_event.empty == False):print "		event is {} it starts {} and ends with {} the height is {}".format(each_event,min(each_event.keys()),max(each_event.keys()),each_event[min(each_event.keys())])
				for i,dt in enumerate(each_event.keys()):
					#print "dt is {}".format(pd.to_datetime(dt))
					dt = pd.to_datetime(dt)
					
					t = timeserie.index.get_loc(dt, method='nearest')#timeserie.iloc[timeserie.index.get_loc(dt, method='nearest')]
					h = np.abs(alts - each_event[each_event.keys()[i]]).argmin()
					
					this_case_clusters.append(clusters_field[h,t])
					this_case_ts.append(dt)
					#if(np.isnan(clusters_field[h,t])):this_cases_clusternames.append("--*--")
					#else:this_cases_clusternames.append(new_labels[int(clusters_field[h,t])])
	    
		#im = ax1.imshow(clusters_field,origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap, norm=Cnorm)
		im = ax1.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap, alpha=alpha)#, norm=Cnorm)
		
		# field that will be field by cluster indexes at teh altitudes, where the centroid is located
		mean_cluster_height = np.full_like(clusters_field, np.nan)
		
		for cluster in np.unique(clusters_field[~np.isnan(clusters_field)]): 
			
			this_clusters_field = np.full_like(clusters_field, np.nan)
			
			alts_tiled = np.tile(alts.T,len(this_clusters_field[0])).reshape([clusters_field.shape[1],clusters_field.shape[0]])
			
			alts_tiled = alts_tiled.T
			
			this_clusters_field[clusters_field==cluster] = alts_tiled[clusters_field==cluster]
			
			if(len(np.where(clusters_field==cluster)[1])>0):
				for j in np.unique(np.where(clusters_field==cluster)[1]):
					
					if (np.sum(np.isfinite((this_clusters_field[:,j])))>2):
					
						mean_height = np.nanmean(this_clusters_field[:,j],axis=0)
						
						diff = (np.abs(this_clusters_field[:,j] - mean_height))
						
						idx = np.where( diff == np.nanmin(diff) )
						
						mean_cluster_height[idx[0][0],j]=cluster
		
		im = ax1.imshow(mean_cluster_height,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='nearest',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
		
		
    qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax1, im=im)
    #qvpp.t_contourplt(T[0:clusters_field.shape[0]], x_lims, alts[0], alts[-1], ax=ax1)
    
    
    
    plt.tight_layout()
    plt.savefig(savefile , format='png',dpi=300)
    print savefile
    print "!!!!!!!!"
    
    #this_case_clusters = np.array(this_case_clusters).astype(int)+1
    
    
    # if (len(this_case_clusters) > 0):
	    # fig, ax2 = plt.subplots(sharey=False,figsize=(10,2))
	    # a = np.array(this_case_clusters).reshape(1, -1)#np.linspace(0, 1, 256).reshape(1, -1)
	    
	   # # df_FAAM = pd.DataFrame(data={"data":np.array(this_case_clusters), "clusters":np.array(this_cases_clusternames)},index=np.array(this_case_ts))
	    # ax2.imshow(a, aspect='auto', cmap=Cmap, origin='lower',vmin=vmin, vmax=vmax)
	    # start, end = ax2.get_xlim()
	    # stepsize = np.abs(end - start)/len(this_cases_clusternames)
	    
	    # for x,s in zip(np.arange(start+stepsize/4,end+stepsize,stepsize),this_cases_clusternames):
			
			# plt.text(x, 0, s, fontsize=12)
	    
	    
	    # ax2.xaxis.set_ticks(np.arange(start, end, stepsize))
	    # ax2.grid(color='k', linestyle='-', linewidth=1, axis ='x')
	    # ts_time = [el.strftime('%H:%M:%S') for el in this_case_ts]
	    # ax2.xaxis.set_ticklabels(ts_time)
	    # #pos = list(ax2.get_position().bounds)
	    
	    # plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30 )
	    # ax2.set_yticks([])
	    # plt.title('Cluters observed by FAAM instruments on {}'.format(this_case_ts[0].strftime("%Y - %m - %d")))
	    # plt.tight_layout()
	    
	    # plt.savefig(savefile[:-4]+"_FAAM.png" , format='png',dpi=300)
	
    # plt.close(savefile)
    
    
	
    return im



def plot_clusters_centroids_in_data(root, timeticks_cases, alts_cases, alt_shape, time_shape,h_list,faam_times_list, colors_list=[], file_name="clusters_in_data.png",output_dir="",hs=None,hvis=None,timeserie_cases=[],verbose=False):
	"""
	input:
		root - 				cluster id
		timeticks_cases - 	list of timeticks per case to use in the plot
		alt_shape -			list of altitude shapes to use in reshaping the data
		time_shape -		list of time shapes to use in reshaping the data
		indexes_cases -		lest of indexes for the right positioning of the labels in the plot
	
	"""
	#if(verbose):
		
	print "we are in plot_clusters_centroids_in_data"
	
	
	(all_childrens_list,all_data,all_labels,all_medoids,all_indxs) = get_active_clusteres_data(root,childrens_list=[])
	cluster_labels = get_all_children_full_name(root,childrens_list=[])
	
	cluster_labels = ([s.replace('cl', '') for s in cluster_labels])
	cluster_labels = ([s.replace('root', '0') for s in cluster_labels])
	new_labels = (['f_cl{}'.format(s+1) for s in np.arange(all_medoids.shape[0])])
	
	two_levels_list = []
	for a_label in cluster_labels:
		a_label = a_label.replace("cl", "")
		if(len(a_label.split('.',2))==3): 
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append('.'.join(a_label.split('.',2)[:-1]))
		else:
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append(a_label)
	first_part_labels = ['.'.join(x.split('.',2)[:-1]) for x in cluster_labels]
	
	clusters_number = len(all_childrens_list)
	
	# number of cases?
	plots_number = np.unique(all_data[:,-2])
	print "plots_number"
	print plots_number
	print "all_data[0,:]"
	print all_data[0,:]
	
	
	
	cluster_colors = []
	
	for p in plots_number.astype(int):
		
		this_labels = all_labels[all_data[:,-2]==p]
		print "this_labels"
		print this_labels
		
		# get the original indxs for this case. The list of indxs is appended list of all indexes for all cases
		if p>1: prev_max_index = sum([(x*y) for (x,y) in zip(alt_shape[:(p-1)],time_shape[:(p-1)])])
		else:prev_max_index = 0
		this_indxs = all_indxs[all_data[:,-2]==p] - prev_max_index
		#print "timeticks_cases[p-1][0] = {}".format(timeticks_cases[p-1][0])
		event = timeticks_cases[p-1][0].strftime('%Y%m%d')
		#print "timeticks_cases[p-1][0] = {}; event is {}, directory is {}{}/plot".format(timeticks_cases[p-1][0],event,output_dir,event)
		
		
		# create directory if it doesn't exist and create a filename with the full path
		if not os.path.exists("{}{}".format(output_dir,event)):os.makedirs("{}{}".format(output_dir,event))
		savefile = "{}{}/case{}_{}_{}".format(output_dir,event,p,clusters_number,file_name)
		savefile_centroids = "{}{}/centroids_case{}_{}_{}".format(output_dir,event,p,clusters_number,file_name)
		savefile_faam = "{}{}/faam_case{}_{}_{}".format(output_dir,event,p,clusters_number,file_name)
		print "savefile_centroids"
		print savefile_centroids
		print "savefile_faam"
		print savefile_faam
		
		if(verbose):print "the plot will be saved in the file {}".format(savefile)
		
		#fig = plt.figure(figsize=(10,5))#True)
		
		#ax = fig.add_subplot(111)
		
		altitudes = alts_cases[p-1]
		print "altitudes"
		print altitudes
		if(colors_list == []):
			cluster_colors = sns.color_palette("hls", clusters_number)
		else:
			white = Color("white")
			main_clusters = np.unique(two_levels_list)
			
			if (len(colors_list) <> len(two_levels_list)):
				for c,a_color in enumerate(colors_list):
					number_of_subclusters = len(re.findall(main_clusters[c], ':'.join(two_levels_list)))
					if(number_of_subclusters > 1):cluster_colors.extend([x.hsl for x in list(white.range_to(Color(hsl=a_color),number_of_subclusters))][:])
					else: cluster_colors.append(a_color)
			else:cluster_colors = colors_list
				
		cluster_map,cluster_norm = qvpp.leveled_discreet_cmap(np.arange(clusters_number),cluster_colors)
		
		#z_colors = sns.color_palette("hls", 12)
		
		# print "!!! cluster_colors length"
		# print cluster_colors
		
		# print "np.arange(clusters_number)"
		# print np.arange(clusters_number)
		
		# print "cluster_map"
		# print cluster_map
		
		# print "cluster_norm"
		# print cluster_norm
		# print "clusters_number"
		
		# print clusters_number
		
		cm = colors.LinearSegmentedColormap.from_list('my_list', cluster_colors, N=clusters_number)
			
		timeticks = timeticks_cases[p-1]
		
		x_lims = mdates.date2num([timeticks[0],timeticks[-1]])
		#[nc.num2date((this_data[:,-1][0]-this_data[:,-1][0]%3600)+x) for x in range(0,int(this_data[:,-1][-1]-this_data[:,-1][0]),3600)]
		#[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
		#print "timeticks"
		#print timeticks
		
		Z = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=int)*np.nan#np.empty_like(this_data[:,0], dtype=int)*np.nan
		T = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=float)*np.nan
		
		#print T.shape
		#print this_labels.shape
		
		Z[this_indxs] = this_labels
		
		# Extract temperature information for this case (p)
		T[this_indxs] = all_data[all_data[:,-2]==p,-4]
		print "T.shape"
		print T.shape
		
		#print all_data[all_data[:,-2]==p,-4]   
		#print all_data[all_data[:,-2]==p,-4].shape 
		#T = all_data[all_data[:,-2]==p,-4]
		
		
		print  "!!!        Temperature !!!!!!! "
		
		#print  T
		print "in Fahrenheit min(T) = {}, max(T) = {}".format(np.nanmin(T), np.nanmax(T))		
		
		#T = (T - 273.15)
		#print T.shape
		Z = Z.reshape(alt_shape[p-1],time_shape[p-1])
		T = T.reshape(alt_shape[p-1],time_shape[p-1])
		
		
		
		#fig, ax = plt.subplots(figsize=(10,5))
	
		#im = qvpp.plot_cluster_field(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),alts=altitudes,ax=ax)
		
		#im = plot_cluster_field_with_faam(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,T_cases[p-1],savefile,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1])
		
		im_centroids = plot_cluster_field_with_faam_and_centroids(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,savefile_centroids,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1],alpha=0.3)
		#plot_cluster_field_with_faam(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,T_cases[p-1],savefile,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1])


		im_faam = plot_faam_with_centroids(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,savefile_faam,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1],faam_h=h_list[p-1],faam_times=faam_times_list[p-1])
		
		
		#if(hs is not None):hs.plot(ax=ax,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
		
		#qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax, im=im)
		
		#print "hs !!!!"
		#print hs
		#if(hvis is not None):hvis.plot(ax=ax,color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    
		#print "qvpp.t_contourplt(T, {}, {}, ax=None)".format(x_lims,np.unique(all_data[:,-3]))
		
		#
		#qvpp.t_contourplt(T, x_lims, altitudes[0], altitudes[-1], ax=ax)
		
		#plt.colorbar(ticks=np.unique(this_labels), label='clusters')
		#plt.tight_layout()
		
		#savefile2 = output_dir + "{}{}_{}".format(output_dir,clusters_number,file_name)
		#plt.savefig(savefile , format='png',dpi=300)
		
		#print ""
		#plt.close(savefile)
		#if(verbose):print "figure %s is ok" %savefile
		#print "function is done"
		
		
	return cluster_colors
	



def plot_hist(bins, values, faam_altitude, faam_temp, cip=15, title="Histogram", filename="./hist.png"):
	#cip_bins_15,hist_cip15_pads,np.mean(height_vis),np.mean(tempr_vis),cip=15,
	
	fig, ax = plt.subplots(sharey=False,figsize=(5,3))
	if(cip==15):width = 10
	elif(cip==100):width = 60
	values = values.astype(int) 
	p = plt.bar(bins,values,width,color='tab:gray')
	
	zip_values = np.array(zip(bins,values))
	#print zip_values
	stats=Stats(zip_values)
	#print vars(stats)
	hist_mean = vars(stats)["mean"]
	hist_std = vars(stats)["stdv"]
	hist_skew = vars(stats)["skew"]
	hist_kurtosis = vars(stats)["kurtosis"]
	
	#plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp-273.15)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f},\ Kurtosis={:.2f}$'.format(hist_mean ,hist_std,hist_skew,hist_kurtosis),{'fontsize': 9})
	plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp-273.15)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f}$'.format(hist_mean ,hist_std,hist_skew),{'fontsize': 9})
	#if(cip==15):
	ax.set_yscale('log')
	##ax.set_ylim(0.,np.max(values))
	plt.ylabel('Mean particles number ($\#\ l^{-1}$)',{'fontsize': 9})
	plt.xlabel('Particles size ($\mu$m)',{'fontsize': 9})
	#plt.xticks(bins)#, ('Bill', 'Fred', 'Mary', 'Sue'))
	#plt.legend((p[0]), ('PADS'),fontsize='small')

	#plt.xticks(bins)#, ('Bill', 'Fred', 'Mary', 'Sue'))
	plt.tight_layout()
	plt.savefig(filename , format='png',dpi=300)
    
	return 
	
	
	
	
class Stats:
    """ Calculate Histogram statistics see https://math.stackexchange.com/questions/857566/how-to-get-the-standard-deviation-of-a-given-histogram-image """
    def __init__(self,histindexed):
        self.n=len(histindexed)
        self.sum=0
        self.prod=0
        self.sqsum=0
        self.cubsum=0
        self.forthsum=0

        for x,y in histindexed:
            self.sum+=y
            self.prod+=x*y
        self.mean=self.prod/self.sum
        for x,y in histindexed:
            dx=x-self.mean
            self.sqsum+=y*dx*dx
            self.cubsum+=y*dx*dx*dx
            self.forthsum+=y*dx*dx*dx
        # sigma squared
        self.variance=self.sqsum/self.sum
        
        self.thirdmoment=self.cubsum/self.sum
        self.forthmoment=self.forthsum/self.sum
        self.stdv=math.sqrt(self.variance)
        self.skew=self.thirdmoment/(np.power(self.stdv,3))
        self.kurtosis=self.forthmoment/(np.power(self.stdv,4))
        
	
def plot_hist_ice_water(bins,values_ice,values_water,faam_altitude,faam_temp,cip=15,title="Histogram",filename="./hist.png"):
	
	fig, ax = plt.subplots(sharey=False,figsize=(5,3))
	if(cip==15):width = 10
	elif(cip==100):width = 60
	values_ice = values_ice.astype(int) 
	values_water = values_water.astype(int)
	p1 = plt.bar(bins, values_water , width, color='y')#, yerr=menStd)
	p2 = plt.bar(bins, values_ice, width, bottom=values_water)#, color='b')# , yerr=womenStd)

	values = values_water + values_ice
	#print "values"
	#print values
	#print('Mean:', np.nansum(values*bins)/np.nansum(values))
	#print('Standard Deviation:', np.std(values))
	#ax.text(0.8, .3, r'$\mu={},\ \sigma={},\ kurtosis={}$'.format(np.mean(values),np.std(values),kurtosis(values)),{'fontsize': 10})
	#ax.text(2.5, 2.5,'at {:.2f} m'.format(faam_altitude),{'fontsize': 8})
	#plt.bar(bins,values,width)
	#hist_mean = np.nansum(values*bins)/np.nansum(values)
	#hist_std = math.sqrt(1/np.nansum(values) * np.sum(values*np.square(bins-hist_mean)))
	zip_values = np.array(zip(bins,values))
	#print zip_values
	stats=Stats(zip_values)
	#print vars(stats)
	hist_mean = vars(stats)["mean"]
	hist_std = vars(stats)["stdv"]
	hist_skew = vars(stats)["skew"]
	hist_kurtosis = vars(stats)["kurtosis"]
	
	#plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp-273.15)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f},\ Kurtosis={:.2f}$'.format(hist_mean ,hist_std,hist_skew,hist_kurtosis),{'fontsize': 9})
	plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp-273.15)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f}$'.format(hist_mean ,hist_std,hist_skew),{'fontsize': 9})
	
	#if(cip==15):
	ax.set_yscale('log')
	##ax.set_ylim(0.,np.max(values))
	plt.ylabel('Mean particles number ($\#\ l^{-1}$)',{'fontsize': 9})
	plt.xlabel('Particles size ($\mu$m)',{'fontsize': 9})
	#plt.xticks(bins)#, ('Bill', 'Fred', 'Mary', 'Sue'))
	plt.legend((p1[0], p2[0]), ('Water content', 'Ice content'),fontsize='small')
	plt.tight_layout()
	plt.savefig(filename , format='png',dpi=300)
	print filename
	#print faam_altitude
    
	return 	
#--------------

def plot_hist_ice_water_edge(bins,values_ice,values_water,values_edge,faam_altitude,faam_temp,cip=15,title="Histogram",filename="./hist.png"):
	print "plot_hist_ice_water_edge"
	fig, ax = plt.subplots(sharey=False,figsize=(5,3))
	if(cip==15):width = 10
	elif(cip==100):width = 60
	values_ice = values_ice.astype(int) 
	values_water = values_water.astype(int)
	p1 = plt.bar(bins, values_water , width, color='y')#, yerr=menStd)
	p2 = plt.bar(bins, values_ice, width, bottom=values_water)#, color='b')# , yerr=womenStd)
	p3 = plt.bar(bins, values_edge, width, bottom=values_water+values_ice, color='grey')#, color='b')# , yerr=womenStd)

	values = values_water + values_ice + values_edge
	print "values"
	print values
	#print('Mean:', np.nansum(values*bins)/np.nansum(values))
	#print('Standard Deviation:', np.std(values))
	#ax.text(0.8, .3, r'$\mu={},\ \sigma={},\ kurtosis={}$'.format(np.mean(values),np.std(values),kurtosis(values)),{'fontsize': 10})
	#ax.text(2.5, 2.5,'at {:.2f} m'.format(faam_altitude),{'fontsize': 8})
	#plt.bar(bins,values,width)
	#hist_mean = np.nansum(values*bins)/np.nansum(values)
	#hist_std = math.sqrt(1/np.nansum(values) * np.sum(values*np.square(bins-hist_mean)))
	zip_values = np.array(zip(bins,values))
	#print zip_values
	stats=Stats(zip_values)
	#print vars(stats)
	hist_mean = vars(stats)["mean"]
	hist_std = vars(stats)["stdv"]
	hist_skew = vars(stats)["skew"]
	hist_kurtosis = vars(stats)["kurtosis"]
	
	#plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp-273.15)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f},\ Kurtosis={:.2f}$'.format(hist_mean ,hist_std,hist_skew,hist_kurtosis),{'fontsize': 9})
	#plt.title('Mean particles size histogram as observed at {:.0f} m AMSL and T = {:.0f} '.format(faam_altitude,faam_temp)+u'\xb0' + 'C\n' +r'Statistics: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f}$'.format(hist_mean ,hist_std,hist_skew),{'fontsize': 9})
	
	plt.title('Mean particles size histogram: $\mu={:.0f},\ \sigma={:.0f},\ Skewness={:.2f}$'.format(hist_mean ,hist_std,hist_skew),{'fontsize': 9})
	#if(cip==15):
	ax.set_yscale('log')
	##ax.set_ylim(0.,np.max(values))
	plt.ylabel('Mean particles number ($\#\ m^{-3}$)',{'fontsize': 9})
	plt.xlabel('Particles size ($\mu$m)',{'fontsize': 9})
	#plt.xticks(bins)#, ('Bill', 'Fred', 'Mary', 'Sue'))
	plt.legend((p1[0], p2[0], p3[0]), ('Liquid', 'Ice','Edge'),fontsize='small')
	plt.tight_layout()
	plt.savefig(filename , format='png',dpi=300)
	print filename
	#print faam_altitude
    
	return 	
#--------------

	
def plot_clusters_in_data(root, timeticks_cases, alts_cases, alt_shape, time_shape, colors_list=[], file_name="clusters_in_data.png",output_dir="",hs=None,hvis=None,timeserie_cases=[],verbose=False):
	"""
	input:
		root - 				cluster id
		timeticks_cases - 	list of timeticks per case to use in the plot
		alt_shape -			list of altitude shapes to use in reshaping the data
		time_shape -		list of time shapes to use in reshaping the data
		indexes_cases -		lest of indexes for the right positioning of the labels in the plot
	
	"""
	if(verbose):
		print "!!!!!! we are in plot_clusters_in_data !!!!!!!"
		
	(all_childrens_list,all_data,all_labels,all_medoids,all_indxs) = get_active_clusteres_data(root,childrens_list=[])
	print "all_data[0] = {}".format(all_data[0,:])
	clusters_number = len(all_childrens_list)
	
	cluster_labels = get_all_children_full_name(root,childrens_list=[])
	F = open('{}/labels_{}.txt'.format(output_dir,clusters_number),'w') 
	F.write('\nOld clusters.\n')
	F.write(', '.join(cluster_labels)) 
	F.write('\nNew clusters.\n')
	
	cluster_labels = ([s.replace('cl', '') for s in cluster_labels])
	cluster_labels = ([s.replace('root', '0') for s in cluster_labels])
	new_labels = (['f_cl{}'.format(s+1) for s in np.arange(all_medoids.shape[0])])
	F.write(', '.join(new_labels))
	F.close()
	if(verbose):
		print ', '.join(cluster_labels)
		print ', '.join(new_labels)
	
	two_levels_list = []
	for a_label in cluster_labels:
		a_label = a_label.replace("cl", "")
		if(len(a_label.split('.',2))==3): 
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append('.'.join(a_label.split('.',2)[:-1]))
		else:
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append(a_label)
	first_part_labels = ['.'.join(x.split('.',2)[:-1]) for x in cluster_labels]
	
	
	plots_number = np.unique(all_data[:,-2])
	cluster_colors = []
	
	series = []
	
	for p in plots_number.astype(int):
		if(verbose): print 'p is {} and new_labels are {}'.format(p,new_labels)
		this_labels = all_labels[all_data[:,-2]==p]
		# get the original indxs for this case. The list of indxs is appended list of all indexes for all cases
		if p>1: prev_max_index = sum([(x*y) for (x,y) in zip(alt_shape[:(p-1)],time_shape[:(p-1)])])
		else:prev_max_index = 0
		this_indxs = all_indxs[all_data[:,-2]==p] - prev_max_index
		#print "timeticks_cases[p-1][0] = {}".format(timeticks_cases[p-1][0])
		event = timeticks_cases[p-1][0].strftime('%Y%m%d')
		#print "timeticks_cases[p-1][0] = {}; event is {}, directory is {}{}/plot".format(timeticks_cases[p-1][0],event,output_dir,event)
		
		
		# create directory if it doesn't exist and create a filename with the full path
		if not os.path.exists("{}{}".format(output_dir,event)):os.makedirs("{}{}".format(output_dir,event))
		savefile = "{}{}/case{}_{}_{}".format(output_dir,event,p,clusters_number,file_name)
		
		if(verbose):print "the plot will be saved in the file {}".format(savefile)
		
		#fig = plt.figure(figsize=(10,5))#True)
		
		#ax = fig.add_subplot(111)
		
		altitudes = alts_cases[p-1]
		if(colors_list == []):
			cluster_colors = sns.color_palette("hls", clusters_number)
		else:
			white = Color("white")
			main_clusters = np.unique(two_levels_list)
			
			if (len(colors_list) <> len(two_levels_list)):
				for c,a_color in enumerate(colors_list):
					number_of_subclusters = len(re.findall(main_clusters[c], ':'.join(two_levels_list)))
					if(number_of_subclusters > 1):cluster_colors.extend([x.hsl for x in list(white.range_to(Color(hsl=a_color),number_of_subclusters))][:])
					else: cluster_colors.append(a_color)
			else:cluster_colors = colors_list
				
		cluster_map,cluster_norm = qvpp.leveled_discreet_cmap(np.arange(clusters_number),cluster_colors)
		
		#z_colors = sns.color_palette("hls", 12)
		
		# print "!!! cluster_colors length"
		# print cluster_colors
		
		# print "np.arange(clusters_number)"
		# print np.arange(clusters_number)
		
		# print "cluster_map"
		# print cluster_map
		
		# print "cluster_norm"
		# print cluster_norm
		# print "clusters_number"
		
		# print clusters_number
		
		cm = colors.LinearSegmentedColormap.from_list('my_list', cluster_colors, N=clusters_number)
			
		timeticks = timeticks_cases[p-1]
		
		x_lims = mdates.date2num([timeticks[0],timeticks[-1]])
		#[nc.num2date((this_data[:,-1][0]-this_data[:,-1][0]%3600)+x) for x in range(0,int(this_data[:,-1][-1]-this_data[:,-1][0]),3600)]
		#[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
		#print "timeticks"
		#print timeticks
		
		Z = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=int)*np.nan#np.empty_like(this_data[:,0], dtype=int)*np.nan
		#T = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=float)*np.nan
		
		#print T.shape
		#print this_labels.shape
		
		Z[this_indxs] = this_labels
		
		# Extract temperature information for this case (p)
		#T[this_indxs] = all_data[all_data[:,-2]==p,-4]
		#print T.shape
		#print all_data[all_data[:,-2]==p,-4].shape    
		#T = all_data[all_data[:,-2]==p,-4]
		
		
		#print  "!!!        Temperature !!!!!!! "
		
		#print  T
		#print "in Fahrenheit min(T) = {}, max(T) = {}".format(np.nanmin(T), np.nanmax(T))		
		
		#T = (T - 273.15)
		#print T.shape
		Z = Z.reshape(alt_shape[p-1],time_shape[p-1])
		#T = T.reshape(alt_shape[p-1],time_shape[p-1])
		
		
		
		#fig, ax = plt.subplots(figsize=(10,5))
	
		#im = qvpp.plot_cluster_field(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),alts=altitudes,ax=ax)
		im = plot_cluster_field_with_faam(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),new_labels,timeticks,savefile,alts=altitudes,hs=hs[p-1],hvis=hvis[p-1],timeserie=timeserie_cases[p-1])
		
		#if(hs is not None):hs.plot(ax=ax,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
		
		#qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax, im=im)
		
		#print "hs !!!!"
		#print hs
		#if(hvis is not None):hvis.plot(ax=ax,color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    
		#print "qvpp.t_contourplt(T, {}, {}, ax=None)".format(x_lims,np.unique(all_data[:,-3]))
		
		#
		#qvpp.t_contourplt(T, x_lims, altitudes[0], altitudes[-1], ax=ax)
		
		#plt.colorbar(ticks=np.unique(this_labels), label='clusters')
		#plt.tight_layout()
		
		#savefile2 = output_dir + "{}{}_{}".format(output_dir,clusters_number,file_name)
		#plt.savefig(savefile , format='png',dpi=300)
		
		#print ""
		#plt.close(savefile)
		#if(verbose):print "figure %s is ok" %savefile
		#print "function is done"
		
		
	return cluster_colors
	
	

def read_long_series(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=None,verbose=False,debug=False,trap_interval=None):
				   #(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height, pd_X,elev,plot=True,output_dir=output_dir,verbose=verbose,debug=debug,trap_interval=None)		
	
	""" 
	Read all data for all cases from input file into one pandas dataset
	"""
	#print "in read_long_series"
	with open(cases_file, 'rb') as csvfile:
		spamreader = csv.DictReader(csvfile, delimiter=',')#reader(csvfile, delimiter=' ', quotechar='|'
		"""
		CSV file into dictionary row per row operation
		"""
		case_id = 1
		timeticks_cases = []
		
		alts_cases = []
		alt_shape = []
		hs_cases = []
		hvis_cases = []
		time_shape = []
		indexes_cases = []
		timeserie_cases = []
		
		sys_phi_list = []
		h_list = []
		
		new_serie = True
		#i_new_serie = True
		timeticks_series = []
		
		alts_series = []
		altseries_shape = []
		hs_series = []
		hvis_series = []
		timeseries_shape = []
		#indexes_series = []
		timeserie_series = []
		indexserie_series = []
		
		serie_id = 0
		#i_serie_id = 0
		
		startserie_time = []
		endserie_time = []
		
		timeticks_series = []
		
		events_in_serie = []
		
		
		for row in spamreader:
			# Create the name of the qvp - file
			event = row['case']
			zoom_interval = [row['start'],row['end']]
			(zoom_start,zoom_end) = event_time(zoom_interval,verbose)
			
			# 
			if(trap_interval is not None):(trap_start,trap_end) = trap_time(event,trap_interval,verbose)
			
			qvp_file = read_file(qvp_directory_name,event,elev,verbose=verbose,debug=debug)	
			
			
			# time in the QVP_file	
			times = qvp_file.variables['Time']
			int_times = np.asarray(times[:]).astype(int)
			
			
			
			faam_file_name = row['faam']
			
			hvis = None
			hs = None
			jd = None
			
			# Altitudes in the QVP file are
			alts = qvp_file.variables['Height']
			
			## QVP xlims
			start_time = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
			end_time = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			if(verbose):
				print "start_time is {}".format(start_time)
				print "end_time is {}".format(end_time)
				
			x_lims = mdates.date2num([start_time,end_time])
			
			#print "!!----timeseries_shape beginning"
			#print timeseries_shape
			#print "len(qvp_file['Time'][:])"
			#print len(qvp_file['Time'][:])
			
			if(new_serie):
				if(serie_id > 0):
					endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
					timeseries_shape.append(len(qvp_file['Time'][:]))
				else:
					startserie_time = [nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)]
					endserie_time = [nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)]
					timeseries_shape = [len(qvp_file['Time'][:])]
				new_serie = False
				serie_id +=1
			else:
				diff = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]
				#print "nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]"
				#print "{} - {}".format(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar),startserie_time[serie_id-1])
				
				#print "diff.days"
				#print diff.days
				if abs(diff.days)>1: 
					serie_id +=1
					startserie_time.append(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar))
					#endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
					endserie_time.append(nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar))
					new_serie = False
					timeseries_shape.append(len(qvp_file['Time'][:]))
				else:
					#startserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)1
					timeseries_shape[serie_id-1] += len(qvp_file['Time'][:])
					endserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			
			
			
			datetime_time = [nc.num2date(x,qvp_file['Time'].units) for x in qvp_file['Time']]
			timeticks = ccplib.get_timeticks(datetime_time)#[begin,first_quater,center,last_quater,end]
			#[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
			
			#print "timeseries_shape end"
			#print timeseries_shape
			#print "[nc.num2date(({}-{}%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int({}-{}),3*7200)]".format(qvp_file['Time'][0],qvp_file['Time'][0],qvp_file['Time'][-1],qvp_file['Time'][0])
			
			
			timeticks_cases.append(timeticks)
			alts = qvp_file['Height']
			if(verbose):
				print "Height alts[:].shape is {}".format(alts[:].shape)
				print "timeticks {}".format(timeticks)
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units).type
				#print timeticks[0].type 
				#print timeticks[1].type
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units)
			# get the height closest to the max_height (we will not look above)
			
			#tempr_qvp = qvp_file['temperature_2']
			
			
			
			upper_indx = np.where(np.abs(alts[:] - max_height) <= 25)[0][0]
			if(verbose):
				print "upper_indx is {}".format(upper_indx)
				print "alts[0:upper_indx].shape"
				print alts[0:upper_indx].shape
			alt_shape.append(len(alts[0:upper_indx]))
			alts_cases.append(alts[0:upper_indx])
			time_shape.append(len(qvp_file['Time'][:]))
			
			#if(timeseries_shape == []):timeseries_shape = len(qvp_file['Time'][:])
			#else:
		#		if
			
			
			hs_cases.append(hs)
			hvis_cases.append(hvis)
			
			timeserie_cases.append(nc.num2date(qvp_file['Time'][:],qvp_file['Time'].units,qvp_file['Time'].calendar))
			
			
			
			
			# X is main array with all dual-pol data, temperature, time and the case number  
			# +2 extra  for height and case number
			X = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),len(field_list)))
			Y = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1))
			if(verbose):
				print "X.shape before read_input_data called for the {} time".format(case_id)
				print X.shape 
				#print "!!!Working with the general matrix to be used in PCA !!!! {}".format(X.shape)
			
			(pd_X,indexes,sys_phi_list,night_mask) = read_input_data_night(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case_id,elev,sys_phi_list=sys_phi_list,verbose=verbose)
										   #read_input_data_night(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case,elev,verbose=False,debug=True,sys_phi_list=[]):
			indexes_cases.append(indexes)
			
			
			if(plot):qvpp.plot_first_4(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=case_id,hmax = max_height,hmin=0,x0=0,x1=-1,count_thr=count_thr,hs=hs,hvis=hvis,sys_phi_list=sys_phi_list,night_mask=night_mask)
			
				#qvpp.plot_standard_4(qvp_file,timeticks,output_dir,case_id=case_id,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,sys_phi_list=sys_phi_list)
			if(verbose):
				print "pd_X.shape after read_input_data called for the {} time".format(case_id)
				print pd_X.shape
			case_id = case_id + 1
			#pd.DataFrame(data = X[~np.isnan(X).any(axis=1)], columns = field_lables_list)
			
	        """
	        Add data from this file to the total dataset
	        """
		
		xserie_lims = mdates.date2num(zip(startserie_time,endserie_time))
		
		n=0
		#print "len(timeserie_cases)"
		#print len(timeserie_cases)
		for i in range(0,len(xserie_lims)):
			diff = endserie_time[i] - startserie_time[i]
			
			timeticks_ser = []
			timeserie_ser = []
			indexserie_ser = []
			
			for j in range(0,diff.days+1):
				
				if (n == 0):
					timeticks_ser = timeticks_cases[n]
					timeserie_ser = timeserie_cases[n]
					#indexserie_ser = indexes_cases[n]
					indexserie_ser = [n+1]
				else:
					
					timeticks_ser = timeticks_ser[0:-1]
					timeticks_ser.extend(timeticks_cases[n])
					timeserie_ser = np.append(timeserie_ser, timeserie_cases[n])#timeserie_ser.extend(timeserie_cases[n])
					#timeserie_ser.extend(timeserie_cases[n])#timeserie_ser.extend(timeserie_cases[n])
					
					#indexserie_ser = np.append(indexserie_ser,indexes_cases[n])
					indexserie_ser.extend([n+1])
				
				print "len(timeserie_cases[{}])".format(n)
				print len(timeserie_cases[n])
				n +=1
				print "n= {} len(timeserie_ser) for j = {}".format(n,j)
				print len(timeserie_ser)
				
				print timeserie_ser[0]
				print timeserie_ser[-1]
					
			if(n == 0):
				timeticks_series = timeticks_ser
				timeserie_series = timeserie_ser 
				indexserie_series = indeserie_series
			else:
				timeticks_series.append(timeticks_ser)# = #list(set(timeticks_series + timeticks_ser))# timeticks_series.append(timeticks_ser)
				#timeserie_series = np.append(timeserie_series,timeserie_ser)
				timeserie_series.append(timeserie_ser)
				#
				#indexserie_series = np.append(indexserie_series,indexserie_ser)
				indexserie_series.append(indexserie_ser)
				print "len(timeserie_ser) for i = {}".format(i)
				print len(timeserie_ser)
				print timeserie_ser[0]
				print timeserie_ser[-1]
				
			#print "!!!!!! ------------ timeserie_series for i = {}".format(i)
			#print indexserie_series
			#print "len(timeserie_series)"
			#print len(timeserie_series)
			#print "len(timeserie_series[i])"
			#print len(timeserie_series[i])
			#print timeserie_series[i][0]
			#print timeserie_series[i][-1]
			
			
		
		start_end_series = zip(startserie_time,endserie_time)
		
		
	return(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,timeseries_shape,indexes_cases,indexserie_series,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeserie_series,h_list,start_end_series)


def read_long_series_withnight(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=None,verbose=False,debug=False,trap_interval=None):
	#(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height, pd_X,elev,plot=True,output_dir=output_dir,verbose=verbose,debug=debug,trap_interval=None)		
	
	""" 
	Read all data for all cases from input file into one pandas dataset
	"""
	print "in read_long_series_withnight"
	#print "cases_file"
	#print cases_file
	#print field_lables_list
	print field_list
	print field_lables_list
	
	
	with open(cases_file, 'rb') as csvfile:
		spamreader = csv.DictReader(csvfile, delimiter=',')#reader(csvfile, delimiter=' ', quotechar='|'
		"""
		CSV file into dictionary row per row operation
		"""
		case_id = 1
		timeticks_cases = []
		
		alts_cases = []
		alt_shape = []
		hs_cases = []
		hvis_cases = []
		time_shape = []
		indexes_cases = []
		timeserie_cases = []
		
		sys_phi_list = []
		h_list = []
		
		new_serie = True
		#i_new_serie = True
		timeticks_series = []
		
		alts_series = []
		altseries_shape = []
		hs_series = []
		hvis_series = []
		timeseries_shape = []
		#indexes_series = []
		timeserie_series = []
		indexserie_series = []
		
		serie_id = 0
		#i_serie_id = 0
		
		startserie_time = []
		endserie_time = []
		
		timeticks_series = []
		
		events_in_serie = []
		
		
		for row in spamreader:
			# Create the name of the qvp - file
			event = row['case']
			print "event"
			print event
			zoom_interval = [row['start'],row['end']]
			(zoom_start,zoom_end) = event_time(zoom_interval,verbose)
			
			# 
			if(trap_interval is not None):(trap_start,trap_end) = trap_time(event,trap_interval,verbose)
			
			qvp_file = read_file(qvp_directory_name,event,elev,verbose=verbose,debug=debug)	
			
			
			# time in the QVP_file	
			times = qvp_file.variables['Time']
			int_times = np.asarray(times[:]).astype(int)
			
			
			
			faam_file_name = row['faam']
			
			hvis = None
			hs = None
			jd = None
			
			# Altitudes in the QVP file are
			alts = qvp_file.variables['Height']
			
			## QVP xlims
			start_time = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
			end_time = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			if(verbose):
				print "start_time is {}".format(start_time)
				print "end_time is {}".format(end_time)
				
			x_lims = mdates.date2num([start_time,end_time])
			
			#print "!!----timeseries_shape beginning"
			#print timeseries_shape
			#print "len(qvp_file['Time'][:])"
			#print len(qvp_file['Time'][:])
			
			if(new_serie):
				if(serie_id > 0):
					endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
					timeseries_shape.append(len(qvp_file['Time'][:]))
				else:
					startserie_time = [nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)]
					endserie_time = [nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)]
					timeseries_shape = [len(qvp_file['Time'][:])]
				new_serie = False
				serie_id +=1
			else:
				diff = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]
				#print "nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]"
				#print "{} - {}".format(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar),startserie_time[serie_id-1])
				
				#print "diff.days"
				#print diff.days
				if abs(diff.days)>1: 
					serie_id +=1
					startserie_time.append(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar))
					#endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
					endserie_time.append(nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar))
					new_serie = False
					timeseries_shape.append(len(qvp_file['Time'][:]))
				else:
					#startserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)1
					timeseries_shape[serie_id-1] += len(qvp_file['Time'][:])
					endserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			
			
			
			datetime_time = [nc.num2date(x,qvp_file['Time'].units) for x in qvp_file['Time']]
			timeticks = ccplib.get_timeticks(datetime_time)#[begin,first_quater,center,last_quater,end]
			#[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
			
			#print "timeseries_shape end"
			#print timeseries_shape
			#print "[nc.num2date(({}-{}%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int({}-{}),3*7200)]".format(qvp_file['Time'][0],qvp_file['Time'][0],qvp_file['Time'][-1],qvp_file['Time'][0])
			
			
			timeticks_cases.append(timeticks)
			alts = qvp_file['Height']
			if(verbose):
				print "Height alts[:].shape is {}".format(alts[:].shape)
				print "timeticks {}".format(timeticks)
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units).type
				#print timeticks[0].type 
				#print timeticks[1].type
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units)
			# get the height closest to the max_height (we will not look above)
			
			#tempr_qvp = qvp_file['temperature_2']
			
			
			
			upper_indx = np.where(np.abs(alts[:] - max_height) <= 25)[0][0]
			if(verbose):
				print "upper_indx is {}".format(upper_indx)
				print "alts[0:upper_indx].shape"
				print alts[0:upper_indx].shape
			alt_shape.append(len(alts[0:upper_indx]))
			alts_cases.append(alts[0:upper_indx])
			time_shape.append(len(qvp_file['Time'][:]))
			
			#if(timeseries_shape == []):timeseries_shape = len(qvp_file['Time'][:])
			#else:
		#		if
			
			
			hs_cases.append(hs)
			hvis_cases.append(hvis)
			
			timeserie_cases.append(nc.num2date(qvp_file['Time'][:],qvp_file['Time'].units,qvp_file['Time'].calendar))
			
			#print "field_list"
			#print field_list
			
			# X is main array with all dual-pol data, temperature, time and the case number  
			# +2 extra  for height and case number
			X = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),len(field_list)))
			#print "X.shape"
			#print X.shape
			Y = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1))
			if(verbose):
				print "X.shape before read_input_data called for the {} time".format(case_id)
				print X.shape 
				#print "!!!Working with the general matrix to be used in PCA !!!! {}".format(X.shape)
			
			#(pd_X,indexes,sys_phi_list) = read_input_data(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case_id,serie=serie_id,sys_phi_list=sys_phi_list,verbose=verbose)
			(pd_X,indexes,sys_phi_list,night_mask) = read_input_data_night(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case_id,elev,sys_phi_list=sys_phi_list,verbose=verbose)
													#read_input_data_night(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case,elev,verbose=False,debug=True,sys_phi_list=[]):
	
			
			indexes_cases.append(indexes)
			
			
			if(plot):qvpp.plot_first_4(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=case_id,hmax = max_height,hmin=0,x0=0,x1=-1,count_thr=count_thr,hs=hs,hvis=hvis,sys_phi_list=sys_phi_list,night_mask=night_mask)
			
				#qvpp.plot_standard_4(qvp_file,timeticks,output_dir,case_id=case_id,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,sys_phi_list=sys_phi_list)
			if(verbose):
				print "pd_X.shape after read_input_data called for the {} time".format(case_id)
				print pd_X.shape
			case_id = case_id + 1
			#pd.DataFrame(data = X[~np.isnan(X).any(axis=1)], columns = field_lables_list)
			
	        """
	        Add data from this file to the total dataset
	        """
		
		xserie_lims = mdates.date2num(zip(startserie_time,endserie_time))
		
		n=0
		#print "len(timeserie_cases)"
		#print len(timeserie_cases)
		for i in range(0,len(xserie_lims)):
			diff = endserie_time[i] - startserie_time[i]
			
			timeticks_ser = []
			timeserie_ser = []
			indexserie_ser = []
			
			for j in range(0,diff.days+1):
				
				if (n == 0):
					timeticks_ser = timeticks_cases[n]
					timeserie_ser = timeserie_cases[n]
					#indexserie_ser = indexes_cases[n]
					indexserie_ser = [n+1]
				else:
					
					timeticks_ser = timeticks_ser[0:-1]
					timeticks_ser.extend(timeticks_cases[n])
					timeserie_ser = np.append(timeserie_ser, timeserie_cases[n])#timeserie_ser.extend(timeserie_cases[n])
					#timeserie_ser.extend(timeserie_cases[n])#timeserie_ser.extend(timeserie_cases[n])
					
					#indexserie_ser = np.append(indexserie_ser,indexes_cases[n])
					indexserie_ser.extend([n+1])
				
				#print "len(timeserie_cases[{}])".format(n)
				#print len(timeserie_cases[n])
				n +=1
				#print "n= {} len(timeserie_ser) for j = {}".format(n,j)
				#print len(timeserie_ser)
				
				#print timeserie_ser[0]
				#print timeserie_ser[-1]
					
			if(n == 0):
				timeticks_series = timeticks_ser
				timeserie_series = timeserie_ser 
				indexserie_series = indeserie_series
			else:
				timeticks_series.append(timeticks_ser)# = #list(set(timeticks_series + timeticks_ser))# timeticks_series.append(timeticks_ser)
				#timeserie_series = np.append(timeserie_series,timeserie_ser)
				timeserie_series.append(timeserie_ser)
				#
				#indexserie_series = np.append(indexserie_series,indexserie_ser)
				indexserie_series.append(indexserie_ser)
				#print "len(timeserie_ser) for i = {}".format(i)
				#print len(timeserie_ser)
				#print timeserie_ser[0]
				#print timeserie_ser[-1]
				
			#print "!!!!!! ------------ timeserie_series for i = {}".format(i)
			#print indexserie_series
			#print "len(timeserie_series)"
			#print len(timeserie_series)
			#print "len(timeserie_series[i])"
			#print len(timeserie_series[i])
			#print timeserie_series[i][0]
			#print timeserie_series[i][-1]
			
			
		
		start_end_series = zip(startserie_time,endserie_time)
		
		
	return(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,timeseries_shape,indexes_cases,indexserie_series,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeserie_series,h_list,start_end_series,qvp_file['Time'])


def fill_in_time_gap(time_diff,count,this_timeserie,the_count,Zc,alt_shape):
	
		# we start to take an action if there is a gap in the time series.
	if (int(max(time_diff)/count.most_common(1)[0][0])>1):
		# number of missing elements
		#clmn_nmb = int(max(time_diff)/the_count)-1
		#timeserie_len = int((this_timeserie[-1] - this_timeserie[0]).seconds/the_count)+1#int(24*60*60/the_count)
		
		# index of the last element before the time gap
		ind = time_diff.index(max(time_diff))
		
		# start to fill in the timeserie by all stqmps from the beginning till the beginning of the gap
		ts = (this_timeserie[0:(ind+1)].tolist())
		
		# fill in the first part of the mask with ones (should be zeros???)
		ts_mask = np.ones(len(ts)).tolist()
		
		# generate list of new timestamps to fill in the gap
		list_toext = pd.date_range(start=this_timeserie[ind]+dt.timedelta(seconds=the_count), end=this_timeserie[time_diff.index(max(time_diff))+1]-dt.timedelta(seconds=the_count/2), freq='{}S'.format(the_count))
		list_toext =list_toext.to_pydatetime() #[dt.datetime.fromtimestamp(timestamp) for timestamp in list_toext]
		
		# addt the list of the new dates to the ts list
		ts.extend(list_toext)
		
		# form Zc array, filling in the part of the data, that were missing as np.nan
		Z_concat = np.concatenate((Zc[:,0:(ind+1)], np.zeros((alt_shape,len(list_toext)))*(np.nan)), axis=1) 
		Z_concat = np.concatenate((Z_concat,Zc[:,(ind+1):]), axis=1) 
		Zc = Z_concat
		
		# add the same number of zeros to the mask (shouljd be ones???)
		ts_mask.extend(np.zeros(len(ts)-len(ts_mask)).tolist())
		
		# fill in the rest of the timeserie
		ts.extend(this_timeserie[(ind+1):].tolist())
		
		# fill the rest of the mask with ones (should be zeros???)
		ts_mask.extend(np.ones(len(this_timeserie[(ind+1):])).tolist())
	else:
		# if there were no gap the ts and the mask stay the same as this_timeserie.
		ts = this_timeserie
		# (should be zeros???)
		ts_mask = np.ones(len(this_timeserie))
	return (Zc,ts,ts_mask)
	
def get_night_hours(c,timeserie_cases,alt_shape):
	
	a = Astral()
	a.solar_depression = "civil"
	
	#set city to Lnodon and change coordinates to Porton Down
	city = a["London"]
	# providing coordinates of Porton Down,  UK
	city.latitude = float(51.1167)
	city.longitude = float(1.7000)
	utc=pytz.UTC
	# get the date of the case/day to select sunrise and the sunset
	#print "timeserie_cases[c-1][0]"
	#print timeserie_cases[c-1][0]
	case_day = timeserie_cases[c-1][0]
	#print "case_day"
	#print case_day
	
	next_day = timeserie_cases[c-1][0]+dt.timedelta(days=1)
	#print "next_day"
	#print next_day
	
	civil_twilight_this_day = city.sun(date=case_day, local=False)
	#print "civil_twilight_this_day"
	#print civil_twilight_this_day
	civil_twilight_next_day = city.sun(date=next_day, local=False)

	#night_hours_mask = (timeserie_cases[c-1] > civil_twilight_this_day['dusk'].replace(tzinfo=None)) & (timeserie_cases[c-1] < civil_twilight_next_day['dawn'].replace(tzinfo=None))
	#morning_hours_mask = (timeserie_cases[c-1] < civil_twilight_this_day['dawn'].replace(tzinfo=None))
	night_hours_mask = (timeserie_cases[c-1] > civil_twilight_this_day['sunset'].replace(tzinfo=None)) & (timeserie_cases[c-1] < civil_twilight_next_day['sunrise'].replace(tzinfo=None))
	morning_hours_mask = (timeserie_cases[c-1] < civil_twilight_this_day['sunrise'].replace(tzinfo=None))
	
	night_hours_mask = np.tile(night_hours_mask, (alt_shape[c-1],1))
	morning_hours_mask  = np.tile(morning_hours_mask, (alt_shape[c-1],1))
	
	hours_mask = np.logical_or(night_hours_mask, morning_hours_mask)
	
	return (hours_mask,night_hours_mask,civil_twilight_this_day['sunset'].replace(tzinfo=None),civil_twilight_next_day['sunrise'].replace(tzinfo=None))


def get_night_hours_foroneday(timeserie):
	
	a = Astral()
	a.solar_depression = "civil"
	
	#set city to Lnodon and change coordinates to Porton Down
	city = a["London"]
	# providing coordinates of Porton Down,  UK
	city.latitude = float(51.1167)
	city.longitude = float(1.7000)
	utc=pytz.UTC
	# get the date of the case/day to select sunrise and the sunset
	print "timeserie"
	print timeserie
	print "timeserie[0]"
	print timeserie[0]
	case_day = timeserie[0]
	next_day = timeserie[0]+dt.timedelta(days=1)
	
	civil_twilight_this_day = city.sun(date=case_day, local=False)
	civil_twilight_next_day = city.sun(date=next_day, local=False)

	night_hours_mask = (timeserie > civil_twilight_this_day['dusk'].replace(tzinfo=None)) & (timeserie < civil_twilight_next_day['dawn'].replace(tzinfo=None))
	morning_hours_mask = (timeserie< civil_twilight_this_day['dawn'].replace(tzinfo=None))
	
	
	#night_hours_mask = np.tile(night_hours_mask, (alt_shape,1))
	#morning_hours_mask  = np.tile(morning_hours_mask, (alt_shape,1))
	
	hours_mask = np.logical_or(night_hours_mask, morning_hours_mask)
	print "np.where(hours_mask == True)"
	print "np.where(hours_mask == True)"
	print np.where(hours_mask == True)
	
	return (hours_mask,night_hours_mask,civil_twilight_this_day['dusk'].replace(tzinfo=None),civil_twilight_next_day['dawn'].replace(tzinfo=None))
	

def get_cluster_colors(clusters_number,two_levels_list,c,cluster_colors,colors_list=[]):
	#print "we generate color list"
	if(colors_list == []):cluster_colors = sns.color_palette("hls", clusters_number)
	else:
		main_clusters = np.unique(two_levels_list)
		if (len(colors_list) <> len(two_levels_list)):
			for c,a_color in enumerate(colors_list):
				number_of_subclusters = len(re.findall(main_clusters[c], ':'.join(two_levels_list)))
				if(number_of_subclusters > 1):cluster_colors.extend([x.hsl for x in list(white.range_to(Color(hsl=a_color),number_of_subclusters))][:])
				else: cluster_colors.append(a_color)
		else:cluster_colors = colors_list
			
	cluster_map,cluster_norm = qvpp.leveled_discreet_cmap(np.arange(clusters_number),cluster_colors)
	
	cm = colors.LinearSegmentedColormap.from_list('my_list', cluster_colors, N=clusters_number)
	return (cluster_map,cluster_norm,cm,cluster_colors)

def get_night_cluster_colors(Z_night,clusters_number,cluster_colors):
	clusters_number_night = len(np.unique(Z_night[~np.isnan(Z_night)]))
	
	if(clusters_number_night !=  clusters_number):
		cluster_colors_night =  [] 
		for j in np.unique(Z_night[~np.isnan(Z_night)]).astype(int):
			if (j==0):cluster_colors_night = [cluster_colors[j]]
			else:cluster_colors_night.append(cluster_colors[j])
	else:cluster_colors_night = cluster_colors
		
	cluster_map_night,cluster_norm_night = qvpp.leveled_discreet_cmap(np.arange(clusters_number_night),cluster_colors_night)
	cm_night = colors.LinearSegmentedColormap.from_list('my_list', cluster_colors_night, N=len(np.unique(Z_night[~np.isnan(Z_night)])))
	
	return (cluster_map_night,cluster_norm_night,cm_night,cluster_colors_night)


def plot_clusters_in_data_long_series(root, timeticks_series, alts_cases, alt_shape, time_shape,timeseries_shape, colors_list=[], file_name="clusters_in_data.png",output_dir="",hs=None,hvis=None,timeserie_cases=[],timeseries_series=[],casesserie =[],verbose=False,not_final = None):
		
	"""
	input:
		root - 				cluster id
		timeticks_cases - 	list of timeticks per case to use in the plot
		alt_shape -			list of altitude shapes to use in reshaping the data
		time_shape -		list of time shapes to use in reshaping the data
		indexes_cases -		lest of indexes for the right positioning of the labels in the plot
	
	"""
	levels = np.arange(200,1200,100)
	
	if(verbose):print "---------------------!!!!!! plot_clusters_in_data_long_series !!!!!!!"
		
	(all_childrens_list,all_data,all_labels,all_medoids,all_indxs) = get_active_clusteres_data(root,childrens_list=[])
	clusters_number = len(all_childrens_list)
	
	cluster_labels = get_all_children_full_name(root,childrens_list=[])
	F = open('{}/labels_{}.txt'.format(output_dir,clusters_number),'w') 
	F.write('\nOld clusters.\n')
	F.write(', '.join(cluster_labels)) 
	F.write('\nNew clusters.\n')
	
	cluster_labels = ([s.replace('cl', '') for s in cluster_labels])
	cluster_labels = ([s.replace('root', '0') for s in cluster_labels])
	if(not_final == None):new_labels = (['f_cl{}'.format(s+1) for s in np.arange(all_medoids.shape[0])])
	else: new_labels = cluster_labels
	F.write(', '.join(new_labels))
	F.close()
	if(verbose):
		print ', '.join(cluster_labels)
		print ', '.join(new_labels)
	
	two_levels_list = []
	for a_label in cluster_labels:
		a_label = a_label.replace("cl", "")
		if(len(a_label.split('.',2))==3): 
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append('.'.join(a_label.split('.',2)[:-1]))
		else:
			#print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			two_levels_list.append(a_label)
	first_part_labels = ['.'.join(x.split('.',2)[:-1]) for x in cluster_labels]
	
	
	plots_number = len(timeticks_series)#np.unique(all_data[:,-2])
	
	cluster_colors = []
	
	series = []
	
	# a = Astral()
	# a.solar_depression = "civil"
	
	# #set city to Lnodon and change coordinates to Porton Down
	# city = a["London"]
	# # providing coordinates of Porton Down,  UK
	# city.latitude = float(51.1167)
	# city.longitude = float(1.7000)
	# utc=pytz.UTC
	
	date_format = mdates.DateFormatter('%b %d\n%H')
	
	# nighttime list for the function output
	list_of_Z_night = []
	
	# Excel file for all data from the night intervals
	xlsxfile = "{}/nigt_serie_data_{}_{}.xlsx".format(output_dir,clusters_number,file_name[:-4])
	
	#with pd.ExcelWriter(xlsxfile) as writer:

		# go through all timeseries (each will be plotted in a separate plot)
	for p in range(0,plots_number):
		xlsxfile = "{}/nigt_serie{}_data_{}_{}.xlsx".format(output_dir,p,clusters_number,file_name[:-4])
	
		with pd.ExcelWriter(xlsxfile) as writer:
				
			if(verbose): print 'p is {} and new_labels are {}'.format(p,new_labels)
			
			# for each plot select cases (days) that will be combined together to form the plot
			cases_list = casesserie[p]
			
			print "cases_list the days participating in the plot"
			print cases_list
			
			# variable to hold all the data belonging to this time series
			current_timeseres_serie = []
			current_timeseres_night_mask = []

			list_dusk = []
			list_dawn = []
			
			# for each case/day in the time serie p
			for c in cases_list:
				alts = alts_cases[0]
				#print alts
				level_ind = np.zeros(len(levels), dtype=int)
				
				# levels indexes in the altitude data
				for l,lev in enumerate(levels):level_ind[int(l)] = np.argmin(np.absolute(alts-lev))
				
				# select the timeline for this case/day
				this_timeserie = timeserie_cases[c-1]
				
				# get the difference between two consequent steps (timestep) in the timeline
				time_diff = [int((t - s).seconds) for s, t in zip(this_timeserie, this_timeserie[1:])]
				
				# count how often each timestep happends in this day
				count = Counter(time_diff)
				
				# varieable to hold the timeseris
				ts=[]
				
				# most common timestep (scanning strategy)
				the_count = count.most_common(1)[0][0] 
				
				# gewt all klabels for this case/day
				this_labels = all_labels[all_data[:,-2]==c]
				
				# get the first index of this case/day in the total array of data
				if c>1:prev_max_index = sum([(x*y) for (x,y) in zip(alt_shape[:(c-1)],time_shape[:(c-1)])])
				else:prev_max_index = 0
				
				# get all indexes of the opints in the total array, that belong to this case/day
				this_indxs = all_indxs[all_data[:,-2]==c] - prev_max_index
				
				# get all altitudes for this day
				altitudes = alts_cases[c-1]
				
				# create an empty flat array to fill in with the labels data
				Zc = np.empty([time_shape[c-1]*alt_shape[c-1]], dtype=int)*np.nan#np.empty_like(this_data[:,0], dtype=int)*np.nan
				
				# fill the data from the total dataset
				Zc[this_indxs] = this_labels
				
				# reshape the falt array and create 2D field with cluster labels with x-axes time and y-axes altitude 
				Zc = Zc.reshape(alt_shape[c-1],time_shape[c-1])
				# Zc is a blok of data for one day "c" from the case list
				
				# check if the data have a time gap and fill it in
				(Zc,ts,ts_mask) = fill_in_time_gap(time_diff,count,this_timeserie,the_count,Zc,alt_shape[c-1])
	
				# apdate timeserie_cases[c-1] with the current ts
				timeserie_cases[c-1] = np.array(ts)
				
				# convert array to integer
				ts_mask = np.array(ts_mask).astype(int)
				
				(hours_mask,night_hours_mask,dusk,dawn) = get_night_hours(c,timeserie_cases,alt_shape)
				if c==0:
					list_dusk = dusk
					list_dawn = dawn
				else:
					list_dusk.append(dusk)
					list_dawn.append(dawn)
				
				if(c == cases_list[0]):
					Z = Zc
					#mZ = mZc
					series = timeserie_cases[c-1]
					current_timeseres_night_mask = night_hours_mask
				else:
					Z = np.append(Z,Zc, axis=1)
					#mZ = np.append(mZ,mZc,axis=1)
					series = np.append(series,timeserie_cases[c-1])
					current_timeseres_night_mask = np.append(current_timeseres_night_mask,hours_mask, axis=1)
				
				# timeline with filled gaps	
				current_timeseres_serie = series
			
			
			#f = nc4.Dataset('./sampleZ.nc','w', format='NETCDF4')
			#tempgrp = f.createGroup('Clusters')
			nc_file='{}Clusters{}_for_serie{}.h5'.format(output_dir,clusters_number,p+1)
			nc = Dataset(nc_file, 'w', format='NETCDF4')#h5py.File('{}Clusters_for_serie{}.h5'.format(output_dir,p+1), 'w')
			nc.createDimension('Height', len(alts))
			nc.createDimension('Time', len(series))
			times = nc.createVariable('Time', np.float64,('Time',))
			heights = nc.createVariable('Height', np.float32, ('Height',))
			time_units = 'seconds since %i-01-01T00:00:00Z' %series[0].year
			number_times = [date2num(stime,time_units,calendar='gregorian') for stime in series]
			heights[:] = alts
			times[:] = number_times
			times.units = time_units
			times.calendar = "gregorian"
			heights.units = "metres above sea level"
			#print "Z.shape"
			#print Z.shape
			Z_nc = nc.createVariable('Clusters', 'u8', ('Height','Time'))
			Z_nigh_mask = nc.createVariable('Mask_night', 'u8', ('Height','Time'))
			
			#intZ = Z
			#intZ[np.isnan(Z)]=-9999
			
			Z_nc[:] = Z.astype(int)
			Z_nigh_mask[:] = current_timeseres_night_mask
			Z_nc.long_name = 'Cluster indexes in the QVP time serie {}'.format(p+1)
			Z_nigh_mask.long_name = 'Nighttime mask for the QVP time serie {}'.format(p+1)
			nc.close()
	    
	
	
			
			#h5f.create_dataset('Clusters', data=Z)
			#h5f.create_dataset('Mask_night', data=current_timeseres_night_mask)
			#print len(series)
			#arr_serie =  np.array(series).astype('<i8')
			#print arr_serie.shape
			#print arr_serie.view('<i8')
			#print np.array(series).view('<i8')
			#dset = h5f.create_dataset('Time_line', arr_serie.shape, '<i8')
			#dset[:,:] = arr_serie.view('<i8')
			#create_dataset('Time_line', data=series)
			#h5f.create_dataset('Altitudes', data=alts)
			#h5f.close()
			
			# update timeseries_series[p] if there were any gaps filled in the data
			if(len(current_timeseres_serie) != len(timeseries_series[p])):
				timeseries_series[p] = current_timeseres_serie
			
			# create directory if it doesn't exist and create a filename with the full path
			if not os.path.exists("{}".format(output_dir)):os.makedirs("{}".format(output_dir))
			
			savefile = "{}/serie{}_{}_{}".format(output_dir,p+1,clusters_number,file_name)
			savefile2 = "{}/nigt_serie{}_{}_{}".format(output_dir,p+1,clusters_number,file_name)
	
	
			if(verbose):print "the plot will be saved in the file {}".format(savefile)
			
			(cluster_map,cluster_norm,cm,cluster_colors) = get_cluster_colors(clusters_number,two_levels_list,c,cluster_colors,colors_list=colors_list)
		
			# need to get timeticks at this level from timeserie_cases[c-1]	
			timeticks = timeticks_series[p]
			
			x_lims = mdates.date2num([timeticks[0],timeticks[-1]])
			
			#im = qvpp.plot_cluster_field(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),alts=altitudes,ax=ax)new_labels
			im = plot_cluster_field_with_faam(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),cluster_labels,timeticks,savefile,alts=altitudes,hs=None,hvis=None,timeserie=timeseries_series[p])
			
			
			Z_night = np.zeros(Z.shape) * np.nan
			Z_night[current_timeseres_night_mask] = Z[current_timeseres_night_mask]
			
			print "____________!!!! Working with the night data to plot characteristics"
			print "Z_night.shape {} inside cycle for plotnumber {}".format(Z_night.shape,p)
			
			
			
			im1 = plot_cluster_field_with_faam(Z_night,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),cluster_labels,timeticks,savefile2,alts=altitudes,hs=None,hvis=None,timeserie=timeseries_series[p])
			
			df_cl1 = []
			df_cl1_night = []
			index = pd.DatetimeIndex(timeseries_series[p])
			
			(cluster_map_night,cluster_norm_night,cm_night,cluster_colors_night) = get_night_cluster_colors(Z_night,clusters_number,cluster_colors)
			
			# figure for level per level plot 
			fig_levels = plt.figure(figsize=(200,100))
			fig_levels,axa_levels = plt.subplots(5, 2)
			
			
			# figure for all levels together
			fig = plt.figure(figsize=(10,5))#True)
			axa = fig.add_subplot(111)
			margin_bottom = []#np.zeros(len(index))
			#clc = 0
			margin_bottom_level = []
			
			# create numpy array with the number of clusters x number of levels x number of cases in the plot dimentions and filled with counts 
			array_of_night_counts = np.zeros((clusters_number,len(level_ind),len(cases_list)+1))
			
			# try to plot cluster per cluster
			for i in range(0,clusters_number):
				print "!!!______________________________________cluster is {}".format(i)
				
				# create empty counts sets of the same size as original Z set
				Counts = np.zeros(Z.shape, dtype=int)
				Counts_night = np.zeros(Z_night.shape, dtype=int)
				Counts_level = np.zeros((len(level_ind),Z_night.shape[0],Z_night.shape[1]), dtype=int)#level_ind
				
				# where Z has this cluster set counts to 1	
				Counts[Z==i] = 1
				Counts_night[Z_night==i] = 1
				low = 0
				for l, lev in enumerate(level_ind):
					Counts_level[l,low:lev,:] = Counts_night[low:lev]
					low = lev
				
				# ??
				if(np.nonzero(Counts_night)[0].size != 0): 
					tmp_Z = Z[Counts>0]
					tmp_Z =tmp_Z[~np.isnan(tmp_Z)]
					tmp_Z_night = Z_night[Counts_night>0]
					tmp_Z_night = tmp_Z_night[(~np.isnan(tmp_Z_night))]
					
				# create df_ck1_night_level and set it based on Counts_level may be per level
				
				
				h=None
				k=None
				leveled_data = []
				
				# if the cluster is the first one
				if(i==0):
					print "-----i==0------"
					df_cl1 = pd.DataFrame(data=np.sum(Counts, axis=0),index=index,columns=[i])#[cluster_labels[i]])#[new_labels[i]])
					df_cl1_night = pd.DataFrame(data=np.sum(Counts_night, axis=0),index=index,columns=[i])#[cluster_labels[i]])#[new_labels[i]])	
					
					print "l, lev in enumerate(levels)"
					print  list(enumerate(levels))
					
					prev_lev = 0
					# for each level
					for l, lev in enumerate(levels):
						df_night_levels = pd.DataFrame(data=np.sum(Counts_level[l], axis=0),index=index,columns=[lev])#[cluster_labels[i]])#[new_labels[i]])
						df_night_levels_sum = df_night_levels.resample('24H', base=12).sum()#, label='right'
						
						array_of_night_counts[i,l,:] = df_night_levels_sum.values.reshape(len(cases_list)+1)
						
						print "!!!!!!!! df_night_levels_sum.index.date"
						print  df_night_levels_sum.index.date
						
						df_night_levels_sum.index = df_night_levels_sum.index.date
						print "df_night_levels_sum.index"
						print df_night_levels_sum.index
						
						lev_values = np.array(df_night_levels_sum.values)
						
						if(l==0):margin_bottom_level = [np.zeros(len(lev_values))]
						else:margin_bottom_level.append(np.zeros(len(lev_values)))
						
						if (l<5):
							axa_levels[l,0].set_title("{} - {} m".format(prev_lev,lev))
							
							prev_lev = lev
							
							df_night_levels_sum[1:].plot.bar(ax=axa_levels[l,0], stacked=True, bottom = margin_bottom_level[l][1:], color=cluster_colors[i], label="f_cl{}".format(i), align='center')#cluster_colors_night[clc]
							if(l<4):
								for ticklabel in axa_levels[l,0].xaxis.get_ticklabels():
									ticklabel.set_visible(False)
							else:
								h, k = axa_levels[l,0].get_legend_handles_labels()
								print "-----------------h = {} and k = {}".format(h,k)
								#axa_levels[l,0].xaxis.set_major_formatter(date_format)
								date_form = mdates.DateFormatter("%m-%d")
								axa_levels[l,0].xaxis.set_major_formatter(date_form)

								axa_levels[l,0].xaxis.set_tick_params(rotation=30, labelsize=8)
								
							axa_levels[l,0].get_legend().remove()

						else:
							axa_levels[l-5,1].set_title("{} - {} m".format(prev_lev,lev))
							
							prev_lev = lev
							
							if((l-5)<4):
								for ticklabel in axa_levels[l-5,1].xaxis.get_ticklabels():
									ticklabel.set_visible(False)
							else:
								date_form = mdates.DateFormatter("%m-%d")
								axa_levels[l-5,1].xaxis.set_major_formatter(date_form)
								
							df_night_levels_sum[1:].plot.bar(ax=axa_levels[l-5,1], stacked=True, bottom = margin_bottom_level[l][1:], color=cluster_colors[i], label="f_cl{}".format(i), align='center')#cluster_colors_night[clc]
							
							axa_levels[l-5,1].get_legend().remove()
							
						margin_bottom_level[l] = lev_values.flatten()
						
						print "margin_bottom_level after plotting cluster {} in level {}".format(i,l)
						print margin_bottom_level
						leveled_data = df_night_levels_sum.values#margin_bottom_level
						print "=1=1=1=1=1= leveled_data"
						print leveled_data
						
				else:
					print "-------else------"
					
					df_cl2 = pd.DataFrame(data=np.sum(Counts, axis=0),index=index,columns=[i])#[cluster_labels[i]])#[new_labels[i]])
					
					df_cl2_night = pd.DataFrame(data=np.sum(Counts_night, axis=0),index=index,columns=[i])#[cluster_labels[i]])#[new_labels[i]])
					
					print "df_cl2_night"
					print df_cl2_night.describe
					df_cl1 = pd.concat([df_cl1,df_cl2], axis=1, sort=False).reindex(df_cl1.index)#.merge(df_cl2)#
					df_cl1_night = pd.concat([df_cl1_night,df_cl2_night], axis=1, sort=False).reindex(df_cl1.index)#.merge(df_cl2)#
					
					
					for l, lev in enumerate(levels):
						#if (l==0) :df_night_levels = pd.DataFrame(data=np.sum(Counts_level[l], axis=0),index=index,columns=[lev])#[cluster_labels[i]])#[new_labels[i]])
						#else:
						df_night_levels = pd.DataFrame(data=np.sum(Counts_level[l], axis=0),index=index,columns=[lev])#[cluster_labels[i]])#[new_labels[i]])
							#df_night_levels = pd.concat([df_night_levels,another_column], axis=1, sort=False).reindex(df_night_levels.index)#.merge(df_cl2)#
						df_night_levels_sum = df_night_levels.resample('24H', base=12).sum()#, label='right'
						lev_values = np.array(df_night_levels_sum.values)
						
						array_of_night_counts[i,l,:] = df_night_levels_sum.values.reshape(len(cases_list)+1)
						
						margin_bottom_level = np.array(margin_bottom_level)
						
						if (l<5):
							df_night_levels_sum[1:].plot.bar(ax=axa_levels[l,0], stacked=True, bottom = margin_bottom_level[l][1:], color=cluster_colors[i], label="f_cl{}".format(i), align='center')#cluster_colors_night[clc]
							axa_levels[l,0].get_legend().remove()
							#axa_levels[l,0].xaxis.set_major_formatter(date_format)
						else:
							df_night_levels_sum[1:].plot.bar(ax=axa_levels[l-5,1], stacked=True, bottom = margin_bottom_level[l][1:], color=cluster_colors[i], label="f_cl{}".format(i), align='center')#cluster_colors_night[clc]
							axa_levels[l-5,1].get_legend().remove()
							#axa_levels[l-5,1].xaxis.set_major_formatter(date_format)
	
						
						margin_bottom_level[l] = margin_bottom_level[l] + lev_values.flatten()
						print "margin_bottom_level[{}] after adding lev_values".format(l)
						print margin_bottom_level[l]
						
					# add another 2D array to 
					leveled_data.append(df_night_levels_sum.values)#margin_bottom_level)
					print "=2=2=2=2=2= leveled_data"
					print leveled_data
						
				del Counts_night
				del Counts
				df_to_print2 = df_cl1_night.resample('24H', base=12).sum()#, label='right'
				
				#df_to_print2.to_excel("./df_to_print2.xlsx")  
				excel_df = df_to_print2[1:].copy()
				print "list_dusk"
				print list_dusk
				print "excel_df.index"
				print excel_df.index
				
				excel_df['dusk'] = list_dusk
				excel_df['dawn'] = list_dawn
				excel_df.to_excel(writer, sheet_name='Time_Serie_{}'.format(p+1))
				
				
				if(i in tmp_Z_night):
					Values = np.array(df_to_print2[i].values)
					print Values
					
					if(margin_bottom == []):margin_bottom = np.zeros(len(Values))
					#print ""
					df_to_print2[i][1:].plot.bar(ax=axa, stacked=True, bottom = margin_bottom[1:], color=cluster_colors[i], label="{}".format(i), align='center')#cluster_colors_night[clc]
					print "margin_bottom"
					print margin_bottom
				
					margin_bottom += Values
					
				
			print "///-------------array_of_night_counts-----------///"
			print array_of_night_counts
			
			# Excel file for all data from the night intervals
			#xlsx_level_file = "{}/nigt_serie_{}_data_levels_{}_{}.xlsx".format(output_dir,p,clusters_number,file_name[:-4])
	
			#with pd.ExcelWriter(xlsx_level_file) as writer_level:
				#df_level = df_cl1_night.resample('24H', base=12).sum()#, label='right'
			low=0
			for l, lev in enumerate(levels):
				#df_to_print2.to_excel("./df_to_print2.xlsx")  
				
				excel_level_df = pd.DataFrame(data=np.transpose(array_of_night_counts[:,l,:]),index=df_to_print2.index)
				excel_level_df = excel_level_df[1:]
				excel_level_df['dusk'] = list_dusk
				excel_level_df['dawn'] = list_dawn
				excel_level_df.to_excel(writer, sheet_name='Time_Serie_{}_level_{}-{}'.format(p+1,low,lev))
				low = lev
				
				#writer_level

			
			#df_to_print = df_cl1_night['2017-06-18 22:00:00':'2017-06-19 06:00:00']
			#print  "df_to_print"
			#print df_to_print
			
			#df_to_print2 = df_cl1_night.resample('24H', base=12).sum()#, label='right'
			#print "df_to_print2"
			#print df_to_print2
			#print df_to_print2.describe
			#print "df_to_print2.index"
			#print df_to_print2.index
			#print "df_to_print2.columns"
			#print df_to_print2.columns
			#print "df_to_print2.values"
			#print df_to_print2.values
			#print "df_to_print2.head()"
			#print df_to_print2.head()
			#print "df_to_print2.columns.values"
			#print df_to_print2.columns.values
			#print "df.index.name"
			#print df_to_print2.index.name
			#df_to_print2.index.name = "Index"
			#print df_to_print2.index.name
			#print  df_to_print2.describe
			#print df_to_print2.columns.values
			#print "cluster_colors"
			#print cluster_colors
			#print "cluster_colors[df_to_print2.columns.values]"
			#print cluster_colors[df_to_print2.columns.values]
			
			#pivot_df = df_to_print2.pivot(columns=df_to_print2.columns.values, values=df_to_print2.values)
			
			#!!!df_to_print2.loc[:,df_to_print2.columns.values].plot.bar(stacked=True, color=cluster_colors_night)
			
			#date_format = mdates.DateFormatter('%b %d\n%H') #mdates.DateFormatter('%Y-%m')#mdates.DateFormatter('%b %d')
			
			#axa.xaxis.set_major_formatter(date_format),colormap=cm_night
			
			#axa = df_to_print2.plot.bar(stacked=True,ax=axa,colormap=cm)#, align='center')#,vmin=np.nanmin(Z), vmax=np.nanmax(Z))
			#axa.xaxis.set_major_formatter(date_format)
			#print "plt.legend(h,k, loc=(1.05, 0.8))"
			#print h,k
			#plt.legend(h,k, loc=(1.05, 0.8))
			fig_levels.savefig("{}/serie{}_barplots_levels_{}_night_image.png".format(output_dir,p+1,clusters_number) , format='png',dpi=100)
			print "savefig({}/serie{}_barplots_levels_{}_night_image.png)".format(output_dir,p+1,clusters_number)
			plt.close(fig_levels)
			
			
			
			
			axa.legend()
			plt.xticks(rotation=10)
			
			plt.draw()
			
			plt.savefig("{}/serie{}_barplot{}_night_image.png".format(output_dir,p+1,clusters_number) , format='png',dpi=100)
			print "savefig({}/serie{}_barplot{}_night_image.png)".format(output_dir,p+1,clusters_number)
			#"{}/nigt_serie{}_{}_{}".format(output_dir,p+1,clusters_number,file_name)
	
			plt.close("{}/serie{}_barplot{}_night_image.png".format(output_dir,p+1,clusters_number))
			
			fig2 = plt.figure(figsize=(10,5))#True)
			
			axa2 = fig2.add_subplot(111)
			df_to_print1 = df_cl1.resample('24H', base=12).sum()#, label='right'
			
			df_to_print1[1:].plot.bar(stacked=True,ax=axa2,colormap=cm, align='center')
			plt.xticks(rotation=10)
			plt.draw()
			
			plt.savefig("{}/serie{}_barplot{}_day_image.png".format(output_dir,p+1,clusters_number) , format='png',dpi=100)
			plt.close("{}/serie{}_barplot{}_day_image.png".format(output_dir,p+1,clusters_number))
			
			
			
					
			#df2.plot.bar(stacked=True)
			
			df_Z = pd.DataFrame(data = Z)#, columns = [altitudes,timeseries_series])#{}
			df_Z_night = pd.DataFrame(data = Z_night)
			if (p==0): list_of_Z_night = [Z_night]
			else:list_of_Z_night.append(Z_night)
			
			
	return cluster_colors,df_Z,list_of_Z_night
	
	

def plot_clusters_in_data_v1(root, timeticks_cases, alts_cases, alt_shape, time_shape, colors_list=[], file_name="clusters_in_data.png",output_dir="",verbose=True):
	"""
	input:
		root - 				cluster id
		timeticks_cases - 	list of timeticks per case to use in the plot
		alt_shape -			list of altitude shapes to use in reshaping the data
		time_shape -		list of time shapes to use in reshaping the data
		indexes_cases -		lest of indexes for the right positioning of the labels in the plot
	
	"""
	if(verbose):print "we are in plot_clusters_in_data"
	(all_childrens_list,all_data,all_labels,all_medoids,all_indxs) = get_active_clusteres_data(root,childrens_list=[])
	
	print all_data.shape
	print all_labels.shape
	print all_medoids.shape
	print all_indxs.shape
	
	cluster_labels = get_all_children_full_name(root,childrens_list=[])
	print cluster_labels
	
	cluster_labels = ([s.replace('cl', '') for s in cluster_labels])
	cluster_labels = ([s.replace('root', '0') for s in cluster_labels])
	new_labels = (['f_cl{}'.format(s+1) for s in np.arange(all_medoids.shape[0])])
	
	
	# two_levels_list = []
	# for a_label in cluster_labels:
		# #a_label = a_label.replace("cl", "")
		# if(len(a_label.split('.',2))==3): 
			# #print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			# two_levels_list.append('.'.join(a_label.split('.',2)[:-1]))
		# else:
			# #print "label name is {} and it should be replaced by {}".format(a_label,a_label.replace("cl", ""))
			# two_levels_list.append(a_label)
	# first_part_labels = ['.'.join(x.split('.',2)[:-1]) for x in cluster_labels]
	
	clusters_number = len(all_childrens_list)
	plots_number = np.unique(all_data[:,-2])
	
	print "!!!plots_number {}".format(plots_number)
	
	cluster_colors = []
	for p in plots_number.astype(int):
		
		this_labels = all_labels[all_data[:,-2]==p]
		# get the original indxs for this case. The list of indxs is appended list of all indexes for all cases
		if p>1: prev_max_index = sum([(x*y) for (x,y) in zip(alt_shape[:(p-1)],time_shape[:(p-1)])])
		else:prev_max_index = 0
		
		print "prev_max_index is {}".format(prev_max_index)
		
		this_indxs = all_indxs[all_data[:,-2]==p] - prev_max_index
		print "this_indxs is {}".format(this_indxs)
		
		#print "timeticks_cases[p-1][0] = {}".format(timeticks_cases[p-1][0])
		event = timeticks_cases[p-1][0].strftime('%Y%m%d')
		print "timeticks_cases[p-1][0] = {}; event is {}, directory is {}{}/plot".format(timeticks_cases[p-1][0],event,output_dir,event)
		if not os.path.exists("{}{}".format(output_dir,event)):os.makedirs("{}{}".format(output_dir,event))
	
		savefile = "{}{}/case{}_{}_{}".format(output_dir,event,p,clusters_number,file_name)
		
		if(verbose):print "the plot will be saved in the file {}".format(savefile)
		
		fig = plt.figure(figsize=(10,5))#True)
		ax = fig.add_subplot(111)
		hmax = max_height
		altitudes = alts_cases[p-1]
		if(colors_list == []):cluster_colors = sns.color_palette("hls", clusters_number)
		# else:
			# white = Color("white")
			# main_clusters = np.unique(two_levels_list)
			
			# if (len(colors_list) <> len(two_levels_list)):
				# for c,a_color in enumerate(colors_list):
					# number_of_subclusters = len(re.findall(main_clusters[c], ':'.join(two_levels_list)))
					# if(number_of_subclusters > 1):cluster_colors.extend([x.hsl for x in list(white.range_to(Color(hsl=a_color),number_of_subclusters))][:])
					# else: cluster_colors.append(a_color)
			# else:cluster_colors = colors_list
				
		cluster_map,cluster_norm =qvpp.leveled_discreet_cmap(np.arange(clusters_number),cluster_colors)
		
		#z_colors = sns.color_palette("hls", 12)
		
		cm = colors.LinearSegmentedColormap.from_list('my_list', cluster_colors, N=clusters_number)
			
		timeticks = timeticks_cases[p-1]
		print "timeticks {}".format(timeticks)
		
		x_lims = mdates.date2num([timeticks[0],timeticks[-1]])
		
		#[nc.num2date((this_data[:,-1][0]-this_data[:,-1][0]%3600)+x) for x in range(0,int(this_data[:,-1][-1]-this_data[:,-1][0]),3600)]
		#[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
		#print "timeticks"
		#print timeticks
		
		Z = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=int)*np.nan#np.empty_like(this_data[:,0], dtype=int)*np.nan
		T = np.empty([time_shape[p-1]*alt_shape[p-1]], dtype=float)*np.nan
		
		Z[this_indxs] = this_labels
		# Extract temperature information for this case (p)
		T[this_indxs] = all_data[all_data[:,-2]==p,-4]
		#print "in Fahrenheit min(T) = {}, max(T) = {}".format(np.nanmin(T), np.nanmax(T))	
		
		print 	"alt_shape[p-1]={},time_shape[p-1]={} and Z.shape = {}".format(alt_shape[p-1],time_shape[p-1],Z.shape)
		
		#T[this_indxs] = (T[this_indxs] - 273.15)
		Z = Z.reshape(alt_shape[p-1],time_shape[p-1])
		T = T.reshape(alt_shape[p-1],time_shape[p-1])
		
		date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
		
		print "Z.shape"
		print Z.shape
		print "T.shape"
		print T.shape
	
		im = qvpp.plot_cluster_field(Z,cm,cluster_norm,x_lims,np.nanmin(Z),np.nanmax(Z),alts=altitudes)
		qvpp.annotate_cluster_plot('clusters',new_labels,hmax,date_format,timeticks,ax=ax, im=im)
		#print "qvpp.t_contourplt(T, {}, {}, ax=None)".format(x_lims,np.unique(all_data[:,-3]))
		
		#qvpp.t_contourplt(T, x_lims, altitudes[0], altitudes[-1], ax=ax)
		
		#plt.colorbar(ticks=np.unique(this_labels), label='clusters')
		
		#savefile2 = output_dir + "{}{}_{}".format(output_dir,clusters_number,file_name)
		plt.savefig(savefile , format='png',dpi=200)
		plt.close()
		if(verbose):print "figure %s is ok" %savefile
		#print "function is done"
		
		
	return cluster_colors
	
	
	
# def read_long_series(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=None,verbose=False,debug=False,trap_interval=None):
	# """ 
	# Read all data for all cases from input file into one pandas dataset
	# """
	# with open(cases_file, 'rb') as csvfile:
		# spamreader = csv.DictReader(csvfile, delimiter=',')#reader(csvfile, delimiter=' ', quotechar='|'
		# """
		# CSV file into dictionary row per row operation
		# """
		# case_id = 1
		# timeticks_cases = []
		
		# alts_cases = []
		# alt_shape = []
		# hs_cases = []
		# hvis_cases = []
		# time_shape = []
		# indexes_cases = []
		# timeserie_cases = []
		
		# sys_phi_list = []
		# h_list = []
		
		# new_serie = True
		# timeticks_series = []
		
		# alts_series = []
		# altseries_shape = []
		# hs_series = []
		# hvis_series = []
		# timeseries_shape = []
		# indexes_series = []
		# timeserie_series = []
		
		# serie_id = 0
		
		# startserie_time = []
		# endserie_time = []
		
		# timeticks_series = []
		
		
		# for row in spamreader:
			# # Create the name of the qvp - file
			# event = row['case']
			# zoom_interval = [row['start'],row['end']]
			# (zoom_start,zoom_end) = event_time(zoom_interval,verbose)
			
			# # 
			# if(trap_interval is not None):(trap_start,trap_end) = trap_time(event,trap_interval,verbose)
			
			# qvp_file = read_file(qvp_directory_name,event,elev,verbose=verbose,debug=debug)	
			
			
			# # time in the QVP_file	
			# times = qvp_file.variables['Time']
			# int_times = np.asarray(times[:]).astype(int)
			
			# faam_file_name = row['faam']
			
			# hvis = None
			# hs = None
			# jd = None
			
			# # Altitudes in the QVP file are
			# alts = qvp_file.variables['Height']
			
			# ## QVP xlims
			# start_time = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
			# end_time = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			# if(verbose):
				# print "start_time is {}".format(start_time)
				# print "end_time is {}".format(end_time)
				# print "qvp_file['Time'].units"
				# print qvp_file['Time'].units
				# print "qvp_file['Time'].calendar"
				# print qvp_file['Time'].calendar
			# x_lims = mdates.date2num([start_time,end_time])
			# #date_format = mdates.DateFormatter('%d %b %Y\n%H:%M:%S')
			
			# print "serie_id"
			# print serie_id
			# print startserie_time
			# print endserie_time
			# if(new_serie):
				# if(serie_id > 0):
					# endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
				# else:
					# startserie_time = [nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)]
					# endserie_time = [nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)]
				# new_serie = False
				# serie_id +=1
			# else:
				# diff = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]
				# #print "nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar) - endserie_time[serie_id-1]"
				# #print "{} - {}".format(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar),startserie_time[serie_id-1])
				
				# #print "diff.days"
				# #print diff.days
				# if abs(diff.days)>1: 
					# serie_id +=1
					# startserie_time.append(nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar))
					# #endserie_time[serie_id - 1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
					# endserie_time.append(nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar))
					# new_serie = False
				# else:
					# #startserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)1
					
					# endserie_time[serie_id-1] = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			
			# #print "startserie_time[{}] is {}".format(serie_id,startserie_time[serie_id])
			# #print "endserie_time[{}] is {}".format(serie_id,endserie_time[serie_id])
			# #print startserie_time
			# #print endserie_time
			
			
			
			
	
			# timeticks = [nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			# print "[nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]"
			# print "[nc.num2date(({}-{}%3600)+x,qvp_file['Time'].units) for x in range(0,3600+({}-{}),3*7200)]".format(qvp_file['Time'][0],qvp_file['Time'][0],qvp_file['Time'][-1],qvp_file['Time'][0])
			
			# timeticks_cases.append(timeticks)
			# alts = qvp_file['Height']
			# if(verbose):
				# print "Height alts[:].shape is {}".format(alts[:].shape)
				# print "timeticks {}".format(timeticks)
				# #print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units).type
				# #print timeticks[0].type 
				# #print timeticks[1].type
				# #print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units)
			# # get the height closest to the max_height (we will not look above)
			
			# #tempr_qvp = qvp_file['temperature_2']
			
			
			
			# upper_indx = np.where(np.abs(alts[:] - max_height) <= 25)[0][0]
			# if(verbose):
				# print "upper_indx is {}".format(upper_indx)
				# print "alts[0:upper_indx].shape"
				# print alts[0:upper_indx].shape
			# alt_shape.append(len(alts[0:upper_indx]))
			# alts_cases.append(alts[0:upper_indx])
			# time_shape.append(len(qvp_file['Time'][:]))
			# hs_cases.append(hs)
			# hvis_cases.append(hvis)
			
			# timeserie_cases.append(nc.num2date(qvp_file['Time'][:],qvp_file['Time'].units,qvp_file['Time'].calendar))
			
			
			
			
			# # X is main array with all dual-pol data, temperature, time and the case number  
			# # +2 extra  for height and case number
			# X = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),len(field_list)))
			# Y = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1))
			# if(verbose):
				# print "X.shape before read_input_data called for the {} time".format(case_id)
				# print X.shape 
				# #print "!!!Working with the general matrix to be used in PCA !!!! {}".format(X.shape)
			
			# (pd_X,indexes,sys_phi_list) = read_input_data(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case_id,sys_phi_list=sys_phi_list,verbose=verbose)
			# indexes_cases.append(indexes)
			# if(plot):qvpp.plot_first_4(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=case_id,hmax = 1500,hmin=0,x0=0,x1=-1,count_thr=count_thr,hs=hs,hvis=hvis,sys_phi_list=sys_phi_list)
			
				# #qvpp.plot_standard_4(qvp_file,timeticks,output_dir,case_id=case_id,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,sys_phi_list=sys_phi_list)
			# if(verbose):
				# print "pd_X.shape after read_input_data called for the {} time".format(case_id)
				# print pd_X.shape
			# case_id = case_id + 1
			# #pd.DataFrame(data = X[~np.isnan(X).any(axis=1)], columns = field_lables_list)
			
	        # """
	        # Add data from this file to the total dataset
	        # """
		# print "zip(startserie_time,endserie_time)"
		# print zip(startserie_time,endserie_time)
		# print "mdates.date2num(zip(startserie_time,endserie_time))"
		# print mdates.date2num(zip(startserie_time,endserie_time))
		# xserie_lims = mdates.date2num(zip(startserie_time,endserie_time))
		# #timeticks = [nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			
		
		
		# n=0
		# #timeticks_ser = []
		# for i in range(0,len(xserie_lims)):
			# diff = endserie_time[i] - startserie_time[i]
			
			# timeticks_ser = []
			# timeserie_ser = []
			
			# for j in range(0,diff.days+1):
				
				# if (n == 0):
					# timeticks_ser = timeticks_cases[n]
					# timeserie_ser = timeserie_cases[n]
				# else:
					
					# timeticks_ser = timeticks_ser[0:-1]
					# timeticks_ser.extend(timeticks_cases[n])
					# timeserie_ser = np.append(timeserie_ser, timeserie_cases[n])#timeserie_ser.extend(timeserie_cases[n])
				
				# n +=1
					
			# if(n == 0):
				# timeticks_series = timeticks_ser
				# timeserie_series = timeserie_ser 
			# else:
				# timeticks_series.append(timeticks_ser)# = #list(set(timeticks_series + timeticks_ser))# timeticks_series.append(timeticks_ser)
				# timeserie_series = np.append(timeserie_series,timeserie_ser)
			
			
			
			
			
		# # print "timeticks_series"
		# # print timeticks_series
		# # print "timeserie_series"
		# # print timeserie_series
		
		
		
		# #print nc.num2date((xserie_lims[0][0]-xserie_lims[0][0]%3600),qvp_file['Time'].units)
		# #print [nc.num2date((xserie_lims[0][0]-xserie_lims[0][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(xserie_lims[0][1]-xserie_lims[0][0]),3*7200)]
			
		
		# #timeticks_series = [nc.num2date((xserie_lims[0][0]-xserie_lims[0][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(xserie_lims[0][1]-xserie_lims[0][0]),3*7200)]
			

	# return(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,h_list,zip(startserie_time,endserie_time))

		
	
def read_cases(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=None,verbose=False,debug=False):
	""" 
	Read all data for all cases from input file into one pandas dataset
	"""
	
	with open(cases_file, 'rb') as csvfile:
		spamreader = csv.DictReader(csvfile, delimiter=',')#reader(csvfile, delimiter=' ', quotechar='|'
		"""
		CSV file into dictionary row per row operation
		"""
		case_id = 1
		timeticks_cases = []
		alts_cases = []
		alt_shape = []
		hs_cases = []
		hvis_cases = []
		time_shape = []
		indexes_cases = []
		timeserie_cases = []
		
		sys_phi_list = []
		h_list = []
		
		new_serie = True
		
		
		for row in spamreader:
			# Create the name of the qvp - file
			event = row['case']
			zoom_interval = [row['start'],row['end']]
			(zoom_start,zoom_end) = event_time(zoom_interval,verbose)
			qvp_file = read_file(qvp_directory_name,event,elev,verbose=verbose,debug=debug)	
			
			
			# time in the QVP_file	
			times = qvp_file.variables['Time']
			int_times = np.asarray(times[:]).astype(int)
			#print "int_times"
			#print int_times
			
			faam_file_name = row['faam']
			#print faam_file_name
			hvis = None
			hs = None
			jd = None
			
			
			
			# Temperature fromQVP_file saved for the plots
			#T = np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	        #            (qvp_file['temperature_2/Means'][:] - 273.15),
	        #            np.nan)
			#T_cases.append(T)
			
			# Altitudes in the QVP file are
			alts = qvp_file.variables['Height']
			
			#print "T"
			#print T
			#print len(T)
			#print T.shape
			#print qvp_file['temperature_2/Means'][:].shape
			#print "mean T"
			#print np.nanmean(T)
			
			
			#print "count_thr"
			
			#print  count_thr
			

			## QVP xlims
			start_time = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
			end_time = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
			if(verbose):
				print "start_time is {}".format(start_time)
				print "end_time is {}".format(end_time)
				print "qvp_file['Time'].units"
				print qvp_file['Time'].units
				print "qvp_file['Time'].calendar"
				print qvp_file['Time'].calendar
			x_lims = mdates.date2num([start_time,end_time])
			#date_format = mdates.DateFormatter('%d %b %Y\n%H:%M:%S')
	
			timeticks = [nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
			timeticks_cases.append(timeticks)
			alts = qvp_file['Height']
			if(verbose):
				print "Height alts[:].shape is {}".format(alts[:].shape)
				print "timeticks {}".format(timeticks)
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units).type
				#print timeticks[0].type 
				#print timeticks[1].type
				#print nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+01,qvp_file['Time'].units)
			# get the height closest to the max_height (we will not look above)
			
			#tempr_qvp = qvp_file['temperature_2']
			
			
			
			upper_indx = np.where(np.abs(alts[:] - max_height) <= 25)[0][0]
			if(verbose):
				print "upper_indx is {}".format(upper_indx)
				print "alts[0:upper_indx].shape"
				print alts[0:upper_indx].shape
			alt_shape.append(len(alts[0:upper_indx]))
			alts_cases.append(alts[0:upper_indx])
			time_shape.append(len(qvp_file['Time'][:]))
			hs_cases.append(hs)
			hvis_cases.append(hvis)
			
			
			timeserie_cases.append(nc.num2date(qvp_file['Time'][:],qvp_file['Time'].units,qvp_file['Time'].calendar))
			
			
			
			
			# X is main array with all dual-pol data, temperature, time and the case number  
			# +2 extra  for height and case number
			X = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),len(field_list)))
			Y = np.zeros((len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1))
			if(verbose):
				print "X.shape before read_input_data called for the {} time".format(case_id)
				print X.shape 
				#print "!!!Working with the general matrix to be used in PCA !!!! {}".format(X.shape)
			
			(pd_X,indexes,sys_phi_list) = read_input_data(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case_id,sys_phi_list=sys_phi_list,verbose=verbose)
			indexes_cases.append(indexes)
			if(plot):qvpp.plot_first_4(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=case_id,hmax = max_height,hmin=0,x0=0,x1=-1,count_thr=count_thr,hs=hs,hvis=hvis,sys_phi_list=sys_phi_list)
			
				#qvpp.plot_standard_4(qvp_file,timeticks,output_dir,case_id=case_id,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,sys_phi_list=sys_phi_list)
			if(verbose):
				print "pd_X.shape after read_input_data called for the {} time".format(case_id)
				print pd_X.shape
			case_id = case_id + 1
			#pd.DataFrame(data = X[~np.isnan(X).any(axis=1)], columns = field_lables_list)
			
	        """
	        Add data from this file to the total dataset
	        """

	return(pd_X,timeticks_cases,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,h_list)

def read_input_data(qvp_file,field_list,upper_indx,count_thr,alts,X,Y,pd_X,field_lables_list,case,serie=None,verbose=False,debug=False,sys_phi_list=[]):
	"""
	Reads input data from the qvp file and generates the panda dataset
	
	Parameters:
    -----------------------------------------
    qvp_file - file-object daataset from the read_file() function
    field_list - list of fields from the input file to be analysed
    upper_indx,count_thr,alts
    
    Returns:
    -----------------------------------------
    file_obj a dataset formed from the opend file
	"""
	# get uPhiDP
	#f = qvp_file['uPhiDP/Means'][0:upper_indx]
	uPhiDP = np.reshape(np.where(qvp_file['uPhiDP/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['uPhiDP/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	#PhiDP = np.reshape(np.where(qvp_file['PhiDP/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['PhiDP/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	# get RhoHV for the non-meteo mask
	RhoHV = np.reshape(np.where(qvp_file['RhoHV/Counts'[:]][0:upper_indx]>count_thr,(qvp_file['RhoHV/Means'][0:upper_indx]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
    
    #  Get height for range calculation to remove nearest ranges up to 400m
	Height = qvp_file['Height'][0:upper_indx]
	
    # Convert height to range assumint 20 deg elevation
	Range = np.cos(np.radians(20)) * Height
	
	sys_phi_list.append(np.nanmedian(uPhiDP))		
	
	for i,field in enumerate(field_list[0:-5]):
		if(debug):print "adding %d field (%s) to the matrix" %(i+1,field)
		field_count = field+'/Counts'
		field_mean = field+'/Means'
		
		f = qvp_file[field_mean][0:upper_indx]
		
		# remove bins near the radar up to 400 m range (due to influence by the side lobes) 
		#f[Range<=400,:] = np.nan
		
		if(debug):print "original shape is "
		if(debug):print f.shape
		X[:,i] = np.reshape(np.where(qvp_file[field_count[:]][0:upper_indx]>count_thr,f,np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
		#X[RhoHV < 0.8,i] = np.nan                
		#if (field is "uPhiDP"):
			#if(debug):
				#print "we are modifying the field"
				#print "np.nanmean(X[:,i]) is {}, i is {}, field is {}".format(np.nanmean(X[:,i]),i,field)
			#X[:,i] = X[:,i] - np.nanmedian(uPhiDP)
			#if(debug):
				#print "np.nanmean(PhiDP)"
				#print np.nanmedian(PhiDP)
				#print "np.nanmean(X[:,i])"
				#print np.nanmedian(X[:,i])
		if (field is "temperature_2"):
			X[:,i] = X[:,i] - 273.15
		if(debug):
			print "X[:,i].shape"
			print X[:,i].shape
			
		#if(debug):plot_hist(X,i,output_dir,field)
		
	
		if(i==0):
			# update Y time-vector
			X[:,-1] = np.reshape(np.where(qvp_file[field_count[:]][0:upper_indx]>count_thr,(qvp_file['Time'][:]),np.nan),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	
		
	
	# add Height field
	#X[:,-3] = np.reshape(np.tile(alts[0:upper_indx], (len(qvp_file['Time'][:]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	X[:,-5] = np.reshape(np.tile(alts[0:upper_indx], (len(qvp_file['Time'][:]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	
	
	
	
	# add Case field
	#X[:,-2] = np.reshape(np.tile(case, (len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	X[:,-4] = np.reshape(np.tile(case, (len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	
	if(debug):print "case index is np.unique(X[:,-3]) {}".format(np.unique(X[:,-3]))
	
	# add Serie field
	X[:,-3] = np.reshape(np.tile(serie, (len(qvp_file['Time'][:])*len(alts[0:upper_indx]),1)).transpose(),(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	if(debug):print "serie index is np.unique(X[:,-2]) {}".format(np.unique(X[:,-2]))
	
	
	
	# get night mask
	#print "nc.num2date(qvp_file['Time'],qvp_file['Time'].units,qvp_file['Time'].calendar)"
	#print [nc.num2date(qvp_file['Time'],qvp_file['Time'].units,qvp_file['Time'].calendar)]
	#print "len([nc.num2date(x,qvp_file['Time'].units,qvp_file['Time'].calendar) for x in qvp_file['Time']])"
	#print len([nc.num2date(x,qvp_file['Time'].units,qvp_file['Time'].calendar) for x in qvp_file['Time']])
	(hours_mask,night_hours_mask,dusk,dawn) = get_night_hours(1,np.asarray([[nc.num2date(x,qvp_file['Time'].units,qvp_file['Time'].calendar) for x in qvp_file['Time']]]),alts[0:upper_indx].shape)
	#print hours_mask,night_hours_mask,dusk,dawn
	#print "night_hours_mask.shape"
	#print night_hours_mask.shape
	#print "np.tile(night_hours_mask, (len(alts[0:upper_indx]),1))"
	#print np.tile(night_hours_mask, (len(alts[0:upper_indx]),1))
	
	
	#print np.tile(night_hours_mask, (len(alts[0:upper_indx]),1)).shape
	
	#print "np.tile(night_hours_mask, (len(alts[0:upper_indx]),1)).transpose()"
	#print np.tile(night_hours_mask, (len(alts[0:upper_indx]),1)).transpose()


	#qvp_file['Time'][:]
	
	print "qvp_file['Time'].units"
	print qvp_file['Time'].units
	# add Night field, which is the 
	X[:,-2] = np.reshape(night_hours_mask,(len(qvp_file['Time'][:])*len(alts[0:upper_indx])))
	if(debug):print "nightmask index is np.unique(X[:,-2]) {}".format(np.unique(X[:,-2]))
	
	
	indexes = ~np.isnan(X).any(axis=1)
	
	# print indexes.shape
	# print X.shape
	# print X
	# print field_lables_list
	# update pd_X with the X data
	if(pd_X is None):
		pd_X = pd.DataFrame(data = X[indexes], columns = field_lables_list)
		
	else: 
		df_X = pd.DataFrame(data = X[indexes], columns = field_lables_list)
		pd_X = pd_X.append(df_X)

	#if(serie is not None):
			#pd_X['case'] = case
	#		pd_X['serie'] = serie
			
	if(debug):
		print "describe the data"
		print pd_X.describe()
		print pd_X.keys().tolist()
	#return
	return (pd_X,indexes,sys_phi_list)
	
def plot_box_plots(df,columns,labels,titles,all_medoids,cluster_labels,output_dir='./',palette=None,not_final=None):
	
	for variable,label,title in zip(columns,labels,titles):
		if (all_medoids.shape[0]>=2):
			if not_final == None: 
				cluster_labels = ([s.replace('cl', '') for s in cluster_labels])
				cluster_labels = ([s.replace('root', '0') for s in cluster_labels])
				new_labels = (['f_cl{}'.format(s+1) for s in np.arange(all_medoids.shape[0])]) 
			else:new_labels = cluster_labels
			
		filename = qvpp.box_plot(df,variable, label, title, new_labels, output_dir,palette=palette)
		
		
	return filename

def get_nightmask_for_all_data(df,qvp_file_Time,cases):
	print "inside get_nightmask_for_all_data"
	#print "df {}".format(df)
	#print "unique days in df['time']"
	#print df['time']
	
	#print "qvp_file_Time.units,qvp_file_Time.calendar"
	#print qvp_file_Time.units,qvp_file_Time.calendar
	time = np.array([nc.num2date(x,qvp_file_Time.units,qvp_file_Time.calendar) for x in df['time']], dtype='datetime64')
	#print "time"
	#print time
	time_dates = np.array([nc.num2date(x,qvp_file_Time.units,qvp_file_Time.calendar) for x in df['time']], dtype='datetime64[D]')
	
	#print "time_dates"
	#print time_dates
	
	df['time_dates'] = time_dates
	df['night'] = np.tile(0,len(time_dates))
	
	# df_time = pd.arrays.DatetimeArray(values=time)
	# df_time_dates = pd.arrays.DatetimeArray(values=time_dates)
	# print "df_time_dates"
	# print df_time_dates
	print "np.unique(time_dates)"
	print np.unique(time_dates)
	for case in np.unique(time_dates):
		this_day_df = df.loc[df['time_dates'] == case]
		#print "this_day_df"
		#print this_day_df
		#print case 
		(hours_mask,night_hours_mask,dusk,dawn) = get_night_hours_foroneday(np.asarray([nc.num2date(x,qvp_file_Time.units,qvp_file_Time.calendar) for x in this_day_df['time']]))
		#print len(case)
		#print "hours_mask"
		#print hours_mask
		#print len(hours_mask)
		#print np.where(hours_mask == 1)
		#print "df.loc[df['time_dates'] == case,'night']"
		#print df.loc[df['time_dates'] == case,'night']
		#print len(df.loc[df['time_dates'] == case,'night'])
		df.loc[df['time_dates'] == case,'night'] = hours_mask # should it be hours or night_hours_mask?
		#print "df.loc[df['night']==1] "
		#print df.loc[df['night']==1] 
		#print len(df.loc[df['night']==1].index) 
		#rose_mask = df.index == 'rose'
	
	
	#print "df.loc[df['night']==1] "
	#print df.loc[df['night']==1] 
	#(hours_mask,night_hours_mask,dusk,dawn) = get_night_hours(1,np.asarray([[nc.num2date(x,qvp_file_Time.units,qvp_file_Time.calendar) for x in df['time']]]),alts_shape)
	
	
	
	return df


def main():
	usage = """%prog [options] /input/directory/QVP_file.nc 
    Use the option -h or --help to get all possible options
    """
    #determine basename of the application
	appbasename = os.path.basename(os.path.splitext(sys.argv[0])[0])

    #determine the basedir  of the application
	basedir = os.path.dirname(sys.argv[0]) + os.sep
	if(basedir == "/"):basedir = "./"

    #the default logfile is located in the same directory as the program and
    #has the program name with the extension '.log'
	dfltlogfile = basedir + appbasename + '.log'

	warnings.filterwarnings("ignore",category=DeprecationWarning)

	# default velues, that will be changed by the parser
	#event = '20170517'
	#zoom_interval = ['20170517T000000','20170517T235959']
	#field_list = ['dBuZ','dBZ','KDP_UKMO','KDP','ZDR','ZDRu','RhoHV','RhoHVu','uPhiDP','PhiDP','SQI','SNR','V','Vu','W','DOP','DOPu']
	

	# parsing options and arguments
	parser = OptionParser(usage=usage)
	parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
	parser.add_option("-e", "--elevation", dest="elev", default = 3, help="elevation angle")
	parser.add_option("-z", "--zoom time-interval", dest="zoom_interval", nargs=2,  help="Time interval tozoom, Formated as HHMMSS, HHMMSS")
	parser.add_option("-n", "--final", action="store_true", dest="not_final", help="The clustering is not final")
	
	(options, files) = parser.parse_args()
	debug = options.debug
	verbose = options.verbose
	elev = options.elev
	#event = options.event
	zoom_interval = None#options.zoom_interval
	not_final = options.not_final
	print "not_final is {}".format(not_final)
	
	#Make aware datetime objects from the start and stop dates (with time zone info included in the object)
	utc = pytz.timezone("UTC")
	zoom_start = None
	zoom_end = None
	if((zoom_interval is not None)):
		if (len(zoom_interval)==2):
			#Make aware datetime objects from the zoom_interval dates (with time zone info included in the object)
			zoom_start = zoom_interval[0]
			#yyyymmdd'T'HHMMSS 
					#dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
						#			int(zoom_interval[0][6:8]), int(zoom_interval[0][9:11]),
						#			int(zoom_interval[0][11:13]),int(zoom_interval[0][13:15]), tzinfo=utc)
			zoom_end = zoom_interval[1]
			#dt.datetime(int(zoom_interval[1][0:4]),int(zoom_interval[1][4:6]),
						#			int(zoom_interval[1][6:8]), int(zoom_interval[1][9:11]),
						#			int(zoom_interval[1][11:13]),int(zoom_interval[1][13:15]), tzinfo=utc)
		else:
			if verbose: print "Zoom interval should be a list with two dates. The event's start and end are used instead."
			
		
	# Option to save the png files
	savefile = True
	
	if(len(sys.argv)<2):print "The number of input values should be 3: clusters hdf-file, original list of cases and output directory."
	# './data' for the file ./data/20170517_QVP/20170517_QVP_90deg.nc
	cluster_file_name = sys.argv[1]
	if verbose: print ("First argument is %s") %cluster_file_name	
	
	#./data/20170517_QVP/plots/
	cases_file = sys.argv[2]
	if verbose: print ("Second argument is %s") %cases_file
	
	#./data/20170517_QVP/plots/
	output_dir = sys.argv[3]
	if verbose: print ("Third argument is %s") %output_dir
	
	qvp_directory_name = sys.argv[4]
	if verbose: print ("Forth argument is %s") %qvp_directory_name
	
	
	# Read hdf file with clusters
	print "reaidng the file {}".format(cluster_file_name)
	root = read_clusters(cluster_file_name)
	
	
	#./data/20170517_QVP/FAAM/core_faam_20170517_v004_r0_c013_1hz.nc
	
	
	(all_childrens_list,all_data,all_labels,all_medoids,all_indxs) = get_active_clusteres_data(root,childrens_list=[])
	cluster_labels = get_all_children_full_name(root,childrens_list=[])
	columns=root.variables.tolist()
	print columns
	print "get_active_clusteres_data"
	print "{}".format(len(all_childrens_list))
	print all_data.shape
	print all_data[0,:]
	#print all_labels.shape
	#print all_medoids
	#print all_medoids.shape
	#print all_indxs
	
	if root is not None:
		plotting_tree(root,parent=None,root_node=None, filename="{}tree_read{}.png".format(output_dir,len(cluster_labels)))
	 	if verbose: print "the png-file {} is saved the loop is finished".format("{}tree_read{}.png".format(output_dir,len(cluster_labels)))
	
	
	columns.append('labels')
	
	data = np.append(np.asarray(all_data),np.reshape(np.asarray(all_labels),(-1,1)), axis=1)
	print "data.shape"
	print data.shape
	
	#columns.append('series')
	#data = np.append(np.asarray(all_data),all_data[:,-1], axis=1)
	
	active_cl_pd = pd.DataFrame(data = data, columns=columns) 
	#print "active_cl_pd.describe()"
	#print active_cl_pd.describe()
	
	df = pd.DataFrame(data = np.asarray(all_data), columns=root.variables.tolist()) 
	df['cluster']=pd.Series((all_labels+1).tolist())
	
	# labels = [r'$\mathrm{Z_{H}}$ (dBZ)',r'$\mathrm{Z_{V}}$ (dBZ)',
	# r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
	# r'$\mathrm{K_{DP}}$ ($^\circ km^{-1}$)',"SNR","Altitude (m)","Cases"]
	# titles = [r'Horizontal Reflectivity Factor',r'Vertical Reflectivity Factor',
	# r'Differential Reflectivity',r'Copolar Correlation Coefficient', 
	# r'Specific Differential Phase',"Signal to noise ratio","Altitude","Cases"]
	
	
	
	labels = [r'$\mathrm{Z_{H}}$ (dBZ)',r'$\mathrm{Z_{V}}$ (dBZ)',
	r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
	r'$\mathrm{K_{DP}}$ ($^\circ km^{-1}$)',"Altitude (m)","Cases"]
	titles = [r'Horizontal Reflectivity Factor',r'Vertical Reflectivity Factor',
	r'Differential Reflectivity',r'Copolar Correlation Coefficient', 
	r'Specific Differential Phase',"Altitude","Cases"]
	# labels = [r'$\mathrm{Z_{V}}$ (dBZ)',
	# r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
	# r'$\mathrm{K_{DP}}$ ($^\circ km^{-1}$)',"Altitude (m)","Cases"]
	# titles = [r'Vertical Reflectivity Factor',
	# r'Differential Reflectivity',r'Copolar Correlation Coefficient', 
	# r'Specific Differential Phase',"Altitude","Cases"]
	
	print "columns"
	print columns
	
	palette = sns.color_palette("hls", all_medoids.shape[0])
	
	
	#filename = qvpp.cat_plot(df,variable, label, title, new_labels, output_dir,palette=palette)
	#print all_medoids.shape
	#print "!!!!!!!!!_________!!!!!!!!!__________!!!!!!!!!!!"
	#print "df {}".format(df)
	
	filename = plot_box_plots(df,columns,labels,titles,all_medoids,cluster_labels,output_dir=output_dir,palette=palette,not_final=not_final)
	
		

	
	print "____________!!!! Working with the night data to plot characteristics"
	#print "df columns {} inside main np.asarray(all_data).shape {} and root.variables.tolist() are {}".format(columns,np.asarray(all_data).shape,root.variables.tolist())
	#print "df {} and df.columns {}".format(df,df.columns)		
	#print "root.variables.tolist()"
	#print root.variables.tolist()
	#print "df.describe()"
	#print df.describe()
	
	
	field_list = ['dBZ','dBZv','ZDR','RhoHV','KDP_UKMO'] #,'PhiDP','KDP_UKMO'
	field_lables_list = ['dBZ','dBZv','ZDR','RhoHV','KDP']#,'PhiDP','KDP'
	#field_list = ['dBZv','ZDR','RhoHV','KDP_UKMO'] #,'PhiDP','KDP_UKMO'
	#field_lables_list = ['dBZv','ZDR','RhoHV','KDP']#,'PhiDP','KDP'
	if('KDP_UKMO' in field_list):field_lables_list = [field.replace('KDP_UKMO', 'KDP') for field in field_list]
	elif('uPfhiDP' in field_list):field_lables_list = [field.replace('uPfhiDP', 'PfhiDP') for field in field_list]
	# Option to save the png files
	savefile = True
	faam = True
	sys_phi = 109
	count_thr = 100#270

	pd_X = None
	#field_list.append("height")
	field_list.append("case")
	field_list.append("serie")
	field_list.append("night")
	field_list.append("time")
	#field_lables_list.append("height")
	field_lables_list.append("case")
	field_lables_list.append("serie")
	field_lables_list.append("night")
	field_lables_list.append("time")
	
	#field_lables_list
	#field_list
	# read all cases from the QVP files
	#(pd_X,timeticks_cases,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs,hvis,timeserie_cases,h_list) = read_cases(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height, pd_X,elev,verbose=verbose,debug=debug,plot=True,output_dir=output_dir)		
	
	#colors_list = plot_clusters_in_data(root, timeticks_cases, alts_cases, alt_shape, time_shape, colors_list=[], file_name="clusters_in_data_loop.png",output_dir=output_dir,hs=hs,hvis=hvis,timeserie_cases=timeserie_cases)
	
	
	#(pd_X,timeticks_series,alts_series,alt_shape,time_shape,indexes_series,sys_phi_list,hs,hvis,timeserie_series,h_list) 
	#(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeserie_series,h_list,start_end_series) 
	#(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeserie_series,h_list,(startserie,endserie))
	#(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,indexes_cases,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeserie_series,h_list,startserie,endserie,onemore)= read_long_series(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=output_dir,verbose=verbose,debug=debug)		
	
	
	
	print "calling read_long_series_withnight"
	#print "field_lables_list"
	#print field_lables_list
	#print "field_list"
	#print field_list
	
	(pd_X,timeticks_cases,timeticks_series,alts_cases,alt_shape,time_shape,timeseries_shape,indexes_cases,indexserie_series,sys_phi_list,hs_cases,hvis_cases,timeserie_cases,timeseries_series,h_list,start_end_series,qvp_file_Time)= read_long_series_withnight(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=output_dir,verbose=verbose,debug=debug,trap_interval=None)
	#read_long_series(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=output_dir,verbose=verbose,debug=debug,trap_interval=None)		
	
	
	
	
	df_night = get_nightmask_for_all_data(df,qvp_file_Time,timeserie_cases)
	print "df_night {}".format(df_night)
	print "df_night.loc[df_night['night'] == 1]"
	print df_night.loc[df_night['night'] == 1]
	
	
	filename_night = plot_box_plots(df_night.loc[df_night['night'] == 1],columns,labels,titles,all_medoids,cluster_labels,output_dir=output_dir+"/night/",palette=palette,not_final=not_final)
	print "filename_night"
	print filename_night
	
	filename_day = plot_box_plots(df_night.loc[df_night['night'] == 0],columns,labels,titles,all_medoids,cluster_labels,output_dir=output_dir+"/day/",palette=palette,not_final=not_final)
	print "filename_day"
	print filename_day
	
	
	
	print "!!!!!!!!!!___________!!!!!!!!!!_____________!!!!!!!!!!!" 																																																#(cases_file,qvp_directory_name,field_lables_list,field_list,count_thr,max_height,pd_X,elev,plot=True,output_dir=None,verbose=False,debug=False)
	print "pd_X {}".format(pd_X)
	
	
	colors_list,df_Z,list_of_Z_night = plot_clusters_in_data_long_series(root, timeticks_series, alts_cases, alt_shape,time_shape,timeseries_shape, colors_list=[], file_name="clusters_in_data_long_series.png",output_dir=output_dir,hs=hs_cases,hvis=hvis_cases,timeserie_cases=timeserie_cases,timeseries_series=timeseries_series,casesserie = indexserie_series,not_final = not_final)#timeserie_cases,timeserie_series=timeserie_series)
	colors_list = plot_clusters_in_data(root, timeticks_cases, alts_cases, alt_shape, time_shape, colors_list=[], file_name="clusters_in_data_loop.png",output_dir=output_dir,hs=hs_cases,hvis=hvis_cases,timeserie_cases=timeserie_cases)
	
	list_Z_night_flatten = []
	for p in range(0,len(list_of_Z_night)):
		print "for different time series p = {}".format(p)
		print list_of_Z_night[p].shape
		print list_of_Z_night[p].flatten()
		if p == 0: list_Z_night_flatten = [list_of_Z_night[p].flatten()]
		else:list_Z_night_flatten.append(list_of_Z_night[p].flatten())
	
	
	#print "df_Z.describe()"
	#print df_Z.describe()
	
	Z = df_Z
	print "Z.shape"
	print Z.shape
	print "len(timeseries_series)"
	print len(timeseries_series)
	unique, counts = np.unique(Z, return_counts=True, axis=1)
	print "unique"
	print unique.shape
	print "counts"
	print counts.shape
	
	#ts_counts = pd.Series(counts[0], timeseries_series)

	
	# plot of centroids found in data (height and temperature to add)
	#colors_list = plot_clusters_centroids_in_data(root, timeticks_cases, alts_cases, alt_shape, time_shape,h_list,None, colors_list=[], file_name="clusters_in_data_loop.png",output_dir=output_dir,hs=hs_cases,hvis=hvis_cases,timeserie_cases=timeserie_cases)
	
	
	
	
	#print "df"
	#print df
	
	
	# for k in np.unique(df['case']):
		# print "case is {}".format(k)
		# a_case = pd.DataFrame(df.query("case=='%s'" %k))
		# for i in np.unique(a_case['cluster']):
			# print "unique cluster is {}".format(i)
			# a_cluster = pd.DataFrame(a_case.query("cluster=='%s'" %i))
			# print "a_cluster.describe()"
			# print a_cluster.describe()
			# counts = []
			# temp = []
			# dBZ = []
			# ZDR = []
			# RhoHV  = []
			# KDP = []
			# height = []
			# temp_std = []
			# dBZ_std = []
			# ZDR_std = []
			# RhoHV_std  = []
			# KDP_std = []
			
			# for level in np.unique(a_cluster['height']):
				# a_level = a_cluster.query("height=='%s'" %level)
				# #print "a_level.describe()"
				# #print a_level.describe()
				# #print a_level.describe()['dBZ']['count']
				# counts.append(a_level.describe()['dBZ_ac']['count'])
				# temp.append(a_level.describe()['temperature_2']['mean'])
				# dBZ.append(a_level.describe()['dBZ_ac']['mean'])
				# ZDR.append(a_level.describe()['ZDR_ac']['mean'])
				# RhoHV.append(a_level.describe()['RhoHV']['mean'])
				# KDP.append(a_level.describe()['KDP']['mean'])
				# height.append(level)
				# temp_std.append(a_level.describe()['temperature_2']['std'])
				# dBZ_std.append(a_level.describe()['dBZ_ac']['std'])
				# ZDR_std.append(a_level.describe()['ZDR_ac']['std'])
				# RhoHV_std.append(a_level.describe()['RhoHV']['std'])
				# KDP_std.append(a_level.describe()['KDP']['std'])
				
				
			# print 'counts {}'.format(counts)
			# #(df,y,labels, title, cluster_labels, outputdir, palette=None)
			# #a_cluster['counts'] = pd.Series(counts)
			# a_cluster_data = pd.DataFrame({"counts" : counts, "temperature_2" : temp, 
			# "dBZ" : dBZ, "ZDR" : ZDR, "RhoHV" : RhoHV, "KDP" : KDP,'height':height,
			# "temperature_2_std" : temp_std, "dBZ_std" : dBZ_std, "ZDR_std" : ZDR_std, "RhoHV_std" : RhoHV_std, "KDP_std" : KDP_std})
			# labels = ['Counts',"$\mathrm{T}$ ($^\circ$C)",r'$\mathrm{Z_{H}}$ (dBZ)',
			# r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{CC}$', 
			# r' $\mathrm{K_{DP}}$ ($deg/km$)']
			# case_label = "case_{}".format(k)
		
			# #file_name = qvpp.line_plot(a_cluster_data,'height',['counts','temperature_2','dBZ','ZDR','RhoHV','KDP'],labels, 'Cluster Description', new_labels[i-1], output_dir, palette=palette[i-1])
			# limits = [[-40,20],[-20,40],[0,1],[0.9,1],[-0.3,0.6]]
			# file_name = qvpp.line_and_error_plot_per_case(a_cluster_data,'height',['counts','temperature_2','dBZ','ZDR','RhoHV','KDP'],labels, limits,'mean per altitude', new_labels[i-1],case_label,output_dir, palette=palette[i-1])
			# print "file_name qvp.line_plot"
			# print file_name
			
	print "np.unique(df['cluster'])"
	print np.unique(df['cluster'])
	
	
	
	

	# lists for the centroids
		
	c_dBZ = []
	c_dBZv = []
	c_ZDR = []
	c_RhoHV = []
	c_KDP = []
	c_height = []
	
	days_c_dBZ = []
	days_c_dBZv = []
	days_c_ZDR = []
	days_c_RhoHV = []
	days_c_KDP = []
	days_c_height = []
	days = []
	days_c_clusters = []
	
	print "np.unique(df['case'])"
	print np.unique(df['case'])
	print "np.unique(df['cluster'])"
	print np.unique(df['cluster'])
	number_of_points = np.zeros((len(np.unique(df['cluster'])),len(np.unique(df['case']))))
	print number_of_points.shape
		
	for i in np.unique(df['cluster']):
		print ("cluster=='%s'" %i)
		a_cluster = pd.DataFrame(df.query("cluster=='%s'" %i))
		print "a_cluster.describe()"
		print a_cluster.describe()
		
		#print "np.unique(a_cluster['height'])"
		#print np.unique(a_cluster['height'])
		
		counts = []
		height = []
		dBZ = []
		dBZv = []
		ZDR = []
		RhoHV  = []
		KDP = []
		height = []
		#height_std = []
		dBZ_std = []
		dBZv_std = []
		ZDR_std = []
		RhoHV_std  = []
		KDP_std = []
		
		
		# append centroids lists
		#print a_cluster.describe()['dBZ_ac']['mean']
		
		
		#c_dBZ.append(a_cluster.describe()['dBZ']['mean'])
		c_dBZv.append(a_cluster.describe()['dBZv']['mean'])
		c_ZDR.append(a_cluster.describe()['ZDR']['mean'])
		c_RhoHV.append(a_cluster.describe()['RhoHV']['mean'])
		c_KDP.append(a_cluster.describe()['KDP']['mean'])
		#c_height.append(a_cluster.describe()['height']['mean'])
		
		days.extend(np.unique(a_cluster['case'])[:])
		print "days"
		print days
		print "np.unique(a_cluster['case'])"
		print np.unique(a_cluster['case'])
		days_c_clusters.extend(np.repeat(i,len(np.unique(a_cluster['case'])[:])))
		
		print "days_c_clusters"
		print days_c_clusters
		
		for a_day in np.unique(a_cluster['case']):
			a_day_cluster = pd.DataFrame(a_cluster.query("case=='%s'" %a_day))
			#days_c_dBZ.append(a_day_cluster.describe()['dBZ']['mean'])
			days_c_dBZv.append(a_day_cluster.describe()['dBZv']['mean'])
			days_c_ZDR.append(a_day_cluster.describe()['ZDR']['mean'])
			days_c_RhoHV.append(a_day_cluster.describe()['RhoHV']['mean'])
			days_c_KDP.append(a_day_cluster.describe()['KDP']['mean'])
			#days_c_height.append(a_day_cluster.describe()['height']['mean'])
			#print number_of_points.shape
			#print "number_of_points[{},{}] = {}".format(int(i)-1,int(a_day)-1,len(a_day_cluster.index))
			number_of_points[int(i)-1,int(a_day)-1] = len(a_day_cluster.index)
			
			
		
		#print "number_of_points"
		#print number_of_points
		
		print "days_c_dBZ"
		print len(days_c_dBZ)
		
		print "days_c_ZDR"
		print len(days_c_ZDR)
		
		print "len(days)"
		print len(days)
		
		
		
		
		for level in np.unique(a_cluster['height']):
			a_level = a_cluster.query("height=='%s'" %level)
			#print "a_level.describe()"
			#print a_level.describe()
			#print a_level.describe()['dBZ']['count']
			counts.append(a_level.describe()['dBZv']['count'])
			#height.append(a_level.describe()['height']['mean'])
			#dBZ.append(a_level.describe()['dBZ']['mean'])
			dBZv.append(a_level.describe()['dBZv']['mean'])
			ZDR.append(a_level.describe()['ZDR']['mean'])
			RhoHV.append(a_level.describe()['RhoHV']['mean'])
			KDP.append(a_level.describe()['KDP']['mean'])
			height.append(level)
			#height_std.append(a_level.describe()['height_std']['std'])
			#dBZ_std.append(a_level.describe()['dBZ']['std'])
			#print a_level.describe()['dBZv']['std']
			#print len(dBZv_std)
			
			dBZv_std.append(a_level.describe()['dBZv']['std'])
			#print len(dBZv_std)
			ZDR_std.append(a_level.describe()['ZDR']['std'])
			RhoHV_std.append(a_level.describe()['RhoHV']['std'])
			KDP_std.append(a_level.describe()['KDP']['std'])
			
			
		#print 'counts {}'.format(counts)
		#(df,y,labels, title, cluster_labels, outputdir, palette=None)
		#a_cluster['counts'] = pd.Series(counts)
		# a_cluster_data = pd.DataFrame({"counts" : counts, 
		# "dBZ" : dBZ,"dBZv" : dBZv, "ZDR" : ZDR, "RhoHV" : RhoHV, "KDP" : KDP,
		# "dBZ_std" : dBZ_std,"dBZv_std" : dBZv_std, "ZDR_std" : ZDR_std, "RhoHV_std" : RhoHV_std, "KDP_std" : KDP_std})
		# labels = ['Counts',"$height",r'$\mathrm{Z_{H}}$ (dBZ)',r'$\mathrm{Z_{V}}$ (dBZ)',
		# r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
		# r' $\mathrm{K_{DP}}$ ($^\circ km^{-1}$)']
		# a_cluster_data = pd.DataFrame({"counts" : counts, 
		# "dBZv" : dBZv, "ZDR" : ZDR, "RhoHV" : RhoHV, "KDP" : KDP,
		# "dBZ_std" : dBZ_std,"dBZv_std" : dBZv_std, "ZDR_std" : ZDR_std, "RhoHV_std" : RhoHV_std, "KDP_std" : KDP_std})
		# labels = ['Counts',"$height",r'$\mathrm{Z_{V}}$ (dBZ)',
		# r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
		# r' $\mathrm{K_{DP}}$ ($^\circ km^{-1}$)']
		
		print "dBZv len {}, counts len {}, ZDR len {},  RhoHV len {}, KDP len {}, dBZv_std len({}),ZDR_std len{}, RhoHV_std len {},KDP_std {}".format(len(dBZv),len(counts),len(ZDR),len(RhoHV),len(KDP),len(dBZv_std),len(ZDR_std),len(RhoHV_std),len(KDP_std))
		
		
		
		a_cluster_data = pd.DataFrame({"counts" : counts, "dBZv" : dBZv, "ZDR" : ZDR, "RhoHV" : RhoHV, "KDP" : KDP,"dBZv_std" : dBZv_std, "ZDR_std" : ZDR_std, "RhoHV_std" : RhoHV_std, "KDP_std" : KDP_std})
		labels = ['Counts',"$height",r'$\mathrm{Z_{V}}$ (dBZ)',
		r'$\mathrm{Z_{DR}}$ (dB)',r'$\mathrm{\rho_{HV}}$', 
		r' $\mathrm{K_{DP}}$ ($^\circ km^{-1}$)']
		
		
		
		
		#file_name = qvpp.line_plot(a_cluster_data,'height',['counts','temperature_2','dBZ','ZDR','RhoHV','KDP'],labels, 'Cluster Description', new_labels[i-1], output_dir, palette=palette[i-1])
		limits = [[-40,20],[-20,40],[-0.5,1],[0.95,1]]#,[-0.3,0.6]]
		#limits = [[-40,20],[-20,40],[-0.5,1],[0.95,1],[-0.3,0.6]]
		#file_name = qvpp.line_and_error_plot(a_cluster_data,'height',['counts','dBZ','ZDR','RhoHV','KDP','height'],labels, limits,'mean per altitude', new_labels[i-1], output_dir, palette=palette[i-1])
		#print "file_name qvp.line_plot"
		#print file_name
	#print cluster_labels
	#print new_labels
	
	#c_clusters = pd.DataFrame({'labels': cluster_labels,'dBZ': c_dBZ,'dBZv': c_dBZv,'ZDR': c_ZDR,'RhoHV': c_RhoHV,'KDP': c_KDP})
	c_clusters = pd.DataFrame({'labels': cluster_labels,'dBZv': c_dBZv,'ZDR': c_ZDR,'RhoHV': c_RhoHV,'KDP': c_KDP})
	
	print c_clusters
	print c_clusters.describe()
	
	#print "len(np.repeat(new_labels,8))"
	#print len(np.repeat(new_labels,8))
	
	print len(days_c_dBZ)
	print len(days_c_ZDR)
	print len(days)
	
	#days_clusters = pd.DataFrame({'clusters': days_c_clusters,'dBZ': days_c_dBZ,'dBZv': days_c_dBZv,'ZDR': days_c_ZDR,'RhoHV': days_c_RhoHV,'KDP': days_c_KDP,'days':days})
	days_clusters = pd.DataFrame({'clusters': days_c_clusters,'dBZv': days_c_dBZv,'ZDR': days_c_ZDR,'RhoHV': days_c_RhoHV,'KDP': days_c_KDP,'days':days})
	
	print "days_clusters"
	print days_clusters.describe()
	
	
	print "!!! --------------number_of_points!!!-------------"
	
	#print number_of_points
	#number_of_points /= number_of_points.sum(keepdims=True)
	#print number_of_points
	#number_of_points_df = pd.DataFrame({'2017-02-01': number_of_points[0,:],'2017-02-03': number_of_points[1,:],'2017-03-03': number_of_points[2,:],'2017-03-22': number_of_points[3,:],'2017-05-17': number_of_points[4,:],'2018-01-24': number_of_points[5,:],'2018-02-13':number_of_points[6,:],'2018-02-14':number_of_points[7,:]})
	#number_of_points_df = pd.DataFrame( number_of_points,columns=['2017-02-01', '2017-02-03', '2017-03-03', '2017-03-22','2017-05-17','2018-01-24','2018-02-13','2018-02-14'])
	#{'2017-02-01': number_of_points[0,:],'2017-02-03': number_of_points[1,:],'2017-03-03': number_of_points[2,:],'2017-03-22': number_of_points[3,:],'2017-05-17': number_of_points[4,:],'2018-01-24': number_of_points[5,:],'2018-02-13':number_of_points[6,:],'2018-02-14':number_of_points[7,:]})
	
	#print "days_clusters"
	#print days_clusters.describe()
	
	
	
	
	def make_spider( row, title, color,with_days=False,handles=[]):
		
		print "make_spider"
		# small_size = 26
		# normal_size = 36
		# big_size = 52
		small_size = 28
		normal_size = 34
		big_size = 36
		linewidth=1.5
		
		# dBZ, ZDR, CC, KDP, T
		#(-20,40),(-1.5,2.0),(0.7,1.0),(-0.3,0.6),(-20,10)
	
		#'KDP', 'RhoHV', 'DOP', 'ZDR', 'dBZ'[0.5,1],
		ranges = np.array([[-7.5,10.0],[0.5,0.9],[-2.0,9.0],[-15,15]])
		
		#labels = [r'$\mathbf{K_{DP}}$',r'$\mathbf{\rho_{HV}}$',"$\mathbf{T}$" ,r'$\mathbf{Z_{DR}}$',r'$\mathbf{Z_{H}}$']
		#labels = [r'$\mathrm{K_{DP}}$',r'$\mathrm{\rho_{HV}}$',r'$\mathrm{Z_{DR}}$',r'$\mathrm{Z_{H}}$']
		labels = [r'$\mathrm{K_{DP}}$',r'$\mathrm{\rho_{HV}}$',r'$\mathrm{Z_{DR}}$',r'$\mathrm{Z_{V}}$']
		
		vas = ["top","top","bottom","bottom","top"]
		has = ["center","left","left","right","right"]
	
		# print "c_clusters"
		# print c_clusters
		# print "ranges"
		# print ranges
	
		
		values = np.array(c_clusters.loc[row].drop('labels').values.flatten().tolist())
		# print "values"
		# print values
		#days_values = []
		#df1.loc[:, df1.loc['a'] > 0]
		
		
		
		scaled_values = (values-ranges[:,0])/(ranges[:,1] - ranges[:,0])
		
		# print "scaled_values"
		# print scaled_values
		
		
		
		values = scaled_values.tolist() #_scale_data(c_clusters.loc[row].drop('labels').values.flatten().tolist(),ranges)
		# print "values"
		# print values 
		values.extend(values[:1])
		#values += values[:1]
		# print values
		# print len(values)
		
		#print "_scale_data"
		#print scaled_values
		#_scale_data(c_clusters.loc[row].drop('labels').values.flatten().tolist(),ranges)

		# number of variable
		categories=list(c_clusters)[:-1]
		N = len(categories)
		#print categories
		
		# We are going to plot the first line of the data frame.
		# But we need to repeat the first value to close the circular graph:
		#values= c_clusters.loc[0].drop('labels').values.flatten().tolist()
		#values += values[:1]
		#print "values"
		#print values
		#print "values"
		#print len(values)
		
		# What will be the angle of each axis in the plot? (we divide the plot / number of variable)
		angles = [n / float(N) * 2 * pi for n in range(N)]
		deg_angels = [n / float(N) * 360 for n in range(N+1)]
		# print "angles"
		# print angles
		angles += angles[:1]
		# print "angles"
		# print angles
		# print "len(angles)"
		# print len(angles)
		
		# Initialise the spider plot
		ax = plt.subplot(5,3,row+1, polar=True )
		 
		# If you want the first axis to be on top:
		ax.set_theta_offset(3*pi / 2)
		ax.set_theta_direction(-1)

		 
		# Draw one axe per variable + add labels labels yet
		print "angles[:-1]"
		print angles[:-1]
		
		print "deg_angels[:-1]"
		print deg_angels[:-1]
		
		print "labels"
		print labels
		#locs, labels = plt.xticks(angles[:-1], labels, color='grey', size=normal_size, frac=1.3)#,va="baseline",ha=has )
		ax.set_thetagrids(deg_angels[:-1], labels=labels, color='grey', size=normal_size, frac=2.5, fontname="Times New Roman")
		
		#print "locs"
		#print locs
		# ang = angles[:-1]
		# for l,lab in enumerate(labels): 
			# loc,label = plt.xticks([ang[l]], [lab], color='grey', size=normal_size,va=vas[l],ha=has[l] )
			# print "locs, label"
			# print loc,label 
		#ax.set_xticks(spacing=5.)
		 
		# Draw ylabels
		#
		ax.set_rlabel_position(0)
		#ax.set_rgrids([0.2,0.4,0.6,0.8,1.0], angle=angles[0])
		
		# print "locs"
		# print locs
		
		# print "labels"
		# print labels
		#plt.yticks([])   dBZ  -20 40 => [-8,4,16,28,40]
		
		#KDP [-0.3,0.6] => -0.12, -.06, 0.24, 0.42, 0.60 
		#np.array([[-7.5,10.0],[0.5,0.9],[0.5,1],[-2.0,9.0],[-15,15]])
		plt.yticks([-8.0,-4.5,1.0,4.5,10.0], [-0.12, -.06, 0.24, 0.42, 0.60], color="grey", size=small_size, fontname="Times New Roman")
		
		ax.set_ylim([0, 1])
		
		#ax2 = polar_twin(ax,ticks=[0.2,0.4,0.6,0.8,1.0],labels=[-0.12, -.06, 0.24, 0.42, 0.60])
		#ax2.set_ylim([0, 1])
		
		# dBZ
		ax.text(angles[4]-radians(5),0.15,'-15',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[4]-radians(5),0.35,'-7.5',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[4]-radians(5),0.55,'0',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[4]-radians(5),0.75,'7.5',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[4]-radians(5),0.95,'15',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.set_ylim([0, 1])
		
		# print "angles[1]"
		# print angles[1]
		#RhoHV
		ax.text(angles[1]-radians(5),0.15,'0.5',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[1]-radians(5),0.35,'0.6',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[1]-radians(5),0.55,'0.7',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[1]-radians(5),0.75,'0.8',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[1]-radians(5),0.95,'0.9',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.set_ylim([0, 1])
		#Altitude
		#ax.text(angles[2]-radians(5),0.15,'0',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		#ax.text(angles[2]-radians(5),0.35,'0.35',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		#ax.text(angles[2]-radians(5),0.55,'0.7',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		#ax.text(angles[2]-radians(5),0.75,'1.05',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		#ax.text(angles[2]-radians(5),0.95,'1.4',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		#ax.set_ylim([0, 1])
		
		# #ZDR
		ax.text(angles[3]-radians(5),0.15,'-14',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[3]-radians(5),0.35,'-8',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[3]-radians(5),0.55,'-2',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[3]-radians(5),0.75,'4',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.text(angles[3]-radians(5),0.95,'10',ha='center',va='center', color="grey", size=small_size, fontname="Times New Roman")
		ax.set_ylim([0, 1])
		
		
		
		h, = ax.plot(angles, values, color=color, linewidth=linewidth, linestyle='solid', label='total dataset centroid')
		ax.fill(angles, values, color=color, alpha=0.2)
		handles.append(copy.copy(h))
		with_days=False

		if(with_days):
			linestyles = [
			     (0, (1, 10)),
			     (0, (1, 1)),
			
			     
			     (0, (5, 5)),
			     (0, (5, 1)),
			
			     (0, (3, 5, 1, 5)),
			     (0, (3, 1, 1, 1)),
			
			     (0, (3, 5, 1, 5, 1, 5)),
			    
			     (0, (3, 1, 1, 1, 1, 1))]
			
			
			labels = ["centroid for 2017-02-01","centroid for 2017-02-03","centroid for 2017-03-03","centroid for 2017-03-22","centroid for 2017-05-17","centroid for 2018-01-24","centroid for 2018-02-13","centroid for 2018-02-14"]
			# 20170201,        20170203 20170303     20170322     20170527       20180124            20180213        20180214    
			#"loosely dotted",'dotted','dashed','densely dashed',"dashdotted","densely dashdotted","dashdotdotted","densely dashdotdotted"]
			a_cluster = days_clusters[days_clusters.clusters==row+1]
			days = np.unique(days_clusters[days_clusters.clusters==row+1]['days'])
			print "days"
			print days
			print "days_clusters[(days_clusters.clusters==row)]"
			print days_clusters[(days_clusters.clusters==row)]
			
			for d in days:
				print "d = {}, row = {}".format(d,row)
				days_values = days_clusters[(days_clusters.clusters==row+1) & (days_clusters.days==d)]
				print "days_values"
				print days_values
				d = int(d)
				days_values = days_values.values.flatten().tolist()
				
				print days_values
				days_values_5 = days_values[0:4]
				print days_values_5
				days_values_5.append(days_values[5])
				scaled_days_values = (days_values_5-ranges[:,0])/(ranges[:,1] - ranges[:,0])
				scaled_days_values = scaled_days_values.tolist() 
				scaled_days_values.extend(scaled_days_values[:1])
				h, = ax.plot(angles, scaled_days_values, color=color, linewidth=linewidth, linestyle=linestyles[d-1], label=labels[d-1])
				handles.append(copy.copy(h))
			
			 
		# Add a title
		plt.title(title, size=big_size, color=color, y=1.03, fontname="Times New Roman Bold")
		return handles
		
 


	# initialize the figure
	# my_dpi=98
	# #plt.figure(figsize=(100, 150))
	# plt.figure(figsize=(3000/my_dpi, 4500/my_dpi), dpi=300)
	# cluster_colors = sns.color_palette("hls", 13)
	my_dpi=300
	#plt.figure(figsize=(100, 150))
	plt.figure(figsize=(9000/my_dpi, 14000/my_dpi), dpi=300)
	cluster_colors = sns.color_palette("hls", 13)
	
	my_palette = plt.cm.get_cmap("hsv", len(c_clusters.index)+1)
	
	# dBZ, ZDR, CC, KDP, T
	#(-20,40),(-1.5,2.0),(0.7,1.0),(-0.3,0.6),(-20,10)
	
	#'KDP', 'RhoHV', 'T', 'ZDR', 'dBZ'
	


 

	# Loop to plot
	print "Loop to plot"
	print range(0, len(c_clusters.index))
	
	for row in range(0, len(c_clusters.index)):
		print "c_clusters['labels'][{}]".format(row)
		print c_clusters['labels'][row]
		
		handles = make_spider( row=row, title='Cluster '+c_clusters['labels'][row], color=my_palette(row),with_days=True,handles=[])
		
	
	ax = plt.subplot(5,3,14, frameon=False, xticks=[], yticks=[])
	for h in handles:
		h.set_color("grey")
		
	#normal_size = 36
	#small_size = 10
	normal_size = 34
	#	big_size = 16
	#	linewidth=0.5
	ax.legend(handles=handles, loc='lower left', fontsize=normal_size,labelspacing=0.6);
	#plt.gca().legend(handles=handles, loc='lower right', fontsize=10, bbox_to_anchor=(1, 0.8))
	
	plt.tight_layout()
	plt.savefig(output_dir + "all_clusters_with_days_spider" , format='png', dpi=my_dpi)
	print output_dir + "all_clusters_with_days_spider"
	
	my_dpi=98
	plt.figure(figsize=(3000/my_dpi, 4500/my_dpi), dpi=my_dpi)
	#plt.figure(figsize=(100, 150))
	
	cluster_colors = sns.color_palette("hls", 13)
	
	my_palette = plt.cm.get_cmap("hsv", len(c_clusters.index)+1)
	
	# dBZ, ZDR, CC, KDP, T
	#(-20,40),(-1.5,2.0),(0.7,1.0),(-0.3,0.6),(-20,10)
	
	#'KDP', 'RhoHV', 'T', 'ZDR', 'dBZ'
	print "c_clusters.index"
	print c_clusters.index

	# Loop to plot
	for row in range(0, len(c_clusters.index)):
		print "c_clusters['labels'][{}]".format(row)
		print c_clusters['labels'][row]
		
		#make_spider( row=row, title='Cluster '+c_clusters['labels'][row], color=my_palette(row))
		make_spider( row=row, title=c_clusters['labels'][row], color=my_palette(row))
		
	
	plt.tight_layout()
	plt.savefig(output_dir + "all_clusters_spider" , format='png', dpi=my_dpi)
	print output_dir + "all_clusters_spider"
	
	
	
	
	
	
	
	print  "!!!!!!! CIP data !!!!!!!!"
	
	CIP_filename = "core-cloud-phy_faam_20180214_v005_r0_c082_cip15.nc"
	
	
	
	
	print conc_cip100_ice
	print conc_cip15_ice
	
	

	#(cip_conc, cip_faam_times) = read_CIP(qvp_directory_name, CIP_filename, event = '20180214', time_start=faam_times_list[-1][-2], time_end=faam_times_list[-1][-1], CIP=15 )
	
	# my_dpi=98
	# plt.figure(figsize=(500/my_dpi, 300/my_dpi), dpi=my_dpi)
	# cluster_colors = sns.color_palette("hls", 13)
	
	# my_palette = plt.cm.get_cmap("hsv", len(c_clusters.index)+1)
	# print "number_of_points_df.columns"
	# print number_of_points_df.columns
	# print number_of_points_df.values
	# print 
	# parallel_coordinates(number_of_points_df.values, number_of_points_df.columns, colormap=my_palette)
	# plt.title("Percents of total number of points per event, per cluster", size=14, fontname="Times New Roman Bold")
		
	

	# plt.tight_layout()
	# plt.savefig("percents_per_cluster_per_case" , format='png', dpi=my_dpi)
	# print "percents_per_cluster_per_case"
	

		
		
		
		
		
		
		
		
		
		
		
	return 0
	
	 
		
	times = qvp_file.variables['Time']
	int_times = np.asarray(times[:]).astype(int)
	#[int(i) for i in times[:]])
	print "experiment with int_times_list"
	print int_times
	faam_dates_in_times = np.asarray(nc.date2num(jd,times.units,calendar=times.calendar)).astype(int)
	print faam_dates_in_times
	int_times_list = np.in1d(faam_dates_in_times,int_times).nonzero()
	print "int_times_list"
	print int_times_list
	print int_times_list[0].shape
	print faam_dates_in_times[int_times_list[0]]
	print h[int_times_list[0]]
	h_in_obs = np.nan * np.ones(shape=int_times.shape)
	h_in_obs[np.in1d(int_times,faam_dates_in_times[int_times_list[0]]).nonzero()[0]] = h[int_times_list[0]]
	h_in_obs[h_in_obs<0] = np.nan
	print h_in_obs
	
	#np.asarray(qvp_file['dBZ/Means']).transpose()
	
	td =  np.asarray(nc.num2date(list(int_times),times.units,times.calendar), dtype=object)
	
	# make panda.series
	hs = pd.Series(h_in_obs,index=td)
	
	
	# TEst
	
	# fig = plt.figure(figsize=(12,4))
	# ax = fig.add_subplot(111)
	# hs.plot(ax=ax,title='%s' % (h.long_name))
	# ax.set_ylabel(h.units)
	
	# savefile0 = "./test.png"
	# plt.savefig(savefile0 , format='png',dpi=200)
	# plt.close()
	
	plt.plot(qvp_file.variables['Time'][1:]-qvp_file.variables['Time'][:-1])
	#plt.plot(ver_file['Time'][1:]-ver_file['Time'][:-1])
	plt.ylim(150,500)
	# savefile1 = "./qvp_file_1minlast.png"
	# if savefile:
	   # plt.savefig(savefile1 , format='png',dpi=200)
	# plt.close()
	# In[12]:
	
	
	#zdrmap,zdrnorm = qvpp.leveled_cmap([-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2,2.5,3],plt.cm.gist_ncar)
	
	#colors = ("#664F2B","#86745D","#AEA295","#F1F1F1","#90A9A0","#537F70","#095E49")
	# colors = ("#001889","#72008D","#AB1488","#D24E71","#E8853A","#ECC000","#DAFF47")
	# zdrmap,zdrnorm = qvpp.leveled_cmap([-1.5,-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],plt.cm.gist_ncar)
	# rhomap,rhonorm = qvpp.leveled_cmap([0.7,0.8,0.9,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1],plt.cm.gist_ncar)
	# kdpmap,kdpnorm = qvpp.leveled_cmap([-0.3,-0.1,0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5],plt.cm.gist_ncar)
	# phimap,phinorm = qvpp.leveled_cmap([-2,-1,0,1,2,3,4,5,6,8,10,12],plt.cm.gist_ncar)
	# Zmap,Znorm = qvpp.leveled_cmap([-20,-10,0,4,8,12,16,20,24,28,32,36,40],plt.cm.gist_ncar)
	# Vmap,Vnorm = qvpp.leveled_cmap([-2,0,0.2,0.4,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.2,4.0,6.0],plt.cm.gist_ncar)
	#____________________________________________________________________________________________________________________________
	zdr_colors = sns.color_palette("hls", 15)
	zdrmap,zdrnorm = qvpp.leveled_discreet_cmap([-1.5,-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],zdr_colors)
	
	
	rhohv_colors = sns.color_palette("hls", 11)
	rhomap,rhonorm = qvpp.leveled_discreet_cmap([0.7,0.8,0.9,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1],rhohv_colors)
	
	kdp_colors = sns.color_palette("hls", 15)
	kdpmap,kdpnorm = qvpp.leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
	
	phi_colors = sns.color_palette("hls", 11)
	phimap,phinorm =qvpp.leveled_discreet_cmap([-2,-1,0,1,2,3,4,5,6,8,10,12],phi_colors)
	
	z_colors = sns.color_palette("hls", 12)
	Zmap,Znorm = qvpp.leveled_discreet_cmap([-20,-10,0,4,8,12,16,20,24,28,32,36,40],z_colors)
	
	v_colors = sns.color_palette("hls", 16)
	Vmap,Vnorm = qvpp.leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
	
	print "leveled_cmap are ready"
	
	# In[13]:
	
	
	## QVP xlims
	start_time = nc.num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
	end_time = nc.num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
	x_lims = mdates.date2num([start_time,end_time])
	
	## FAAM xlims
	faam_start_time = nc.num2date(faam_file['Time'][0],faam_file['Time'].units,faam_file['Time'].calendar)
	faam_end_time = nc.num2date(faam_file['Time'][-1],faam_file['Time'].units,faam_file['Time'].calendar)
	faam_x_lims = mdates.date2num([faam_start_time,faam_end_time])
	
	## Vertical xlims
	vstart_time = nc.num2date(ver_file['Time'][0],ver_file['Time'].units,qvp_file['Time'].calendar)
	vend_time = nc.num2date(ver_file['Time'][-1],ver_file['Time'].units,qvp_file['Time'].calendar)
	vx_lims = mdates.date2num([vstart_time,vend_time])
	
	date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
	
	print " date_format is ready"
	# In[14]:
	
	
	sys_phi = 109
	zdr_cal = +2.1
	timeticks = [nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,
	                      qvp_file['Time'].units) for x in range(0,3600+int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3*7200)]
	faam_timeticks = [nc.num2date((faam_file['Time'][0]-faam_file['Time'][0]%3600)+x,
	                      faam_file['Time'].units) for x in range(0,3600+int(faam_file['Time'][-1]-faam_file['Time'][0]),3*7200)]
	
	print "faam_timeticks are ready"
	print faam_timeticks
	
	# In[15]:
	
	
	# In[16]:
	
	# In[17]:
	
	
	fig, ax_arr = plt.subplots(5, 1, sharey=False,figsize=(17,22))#True)
	print ax_arr
	#plt.figure(figsize=(17,22))
	hmax = max_height
	count_thr = 134
	
	## Plot Reflectivity - Subplot 1 ----------------------------------------------
	#plt.subplot(511)
	#qvpp.plot_faam_data(faam_file,'ALT_GIN')
	hs.plot(ax=ax_arr[0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight')
	im = qvpp.plot_field(qvp_file,'dBuZ',Znorm,Zmap,x_lims,-20,40,count_threshold=count_thr,ax=ax_arr[0])
	qvpp.annotate_field_plot('Reflectivity (dBZ)',Znorm,hmax,date_format,timeticks,ax=ax_arr[0], im=im)
	#qvpp.z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[0])
	
	print "subplot(511)"
	
	# ## Plot ZDR - Subplot 2 -------------------------------------------------------
	# plt.subplot(512)
	# qvpp.plot_field(qvp_file,'ZDRu',zdrnorm,zdrmap,x_lims,-1,3,offset=zdr_cal,count_threshold=count_thr)
	# qvpp.annotate_field_plot('ZDRu (dB)',zdrnorm,hmax,date_format,timeticks)
	# qvpp.z_contourplt(qvp_file, count_thr, x_lims)
	
	## Plot ZDR - Subplot 2 -------------------------------------------------------
	#plt.subplot(512)
	hs.plot(ax=ax_arr[1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight')
	im = qvpp.plot_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-1.5,2,offset=0,count_threshold=count_thr, ax=ax_arr[1])
	qvpp.annotate_field_plot('ZDR (dB)',zdrnorm,hmax,date_format,timeticks, ax=ax_arr[1], im=im)
	#qvpp.z_contourplt(qvp_file, count_thr, x_lims, ax=ax_arr[1])
	
	
	print "subplot(512)"
	
	## Plot RhoHV - Subplot 3 ----------------------------------------------------
	#plt.subplot(513)
	hs.plot(ax=ax_arr[2],color = '0.25', alpha=0.5, lw=3,label='FAAM flight')
	im = qvpp.plot_field(qvp_file,'RhoHVu',rhonorm,rhomap,x_lims,0.7,1,count_threshold=count_thr, ax=ax_arr[2])
	qvpp.annotate_field_plot('CC',rhonorm,hmax,date_format,timeticks, ax=ax_arr[2], im=im)
	#qvpp.z_contourplt(qvp_file, count_thr,x_lims, ax=ax_arr[2])
	
	
	print "subplot(513)"
	
	## Plot KDP - Subplot 4 ----------------------------------------------------
	#plt.subplot(514)
	hs.plot(ax=ax_arr[3],color = '0.25', alpha=0.5, lw=3,label='FAAM flight')
	im = qvpp.plot_field(qvp_file,'KDP_UKMO',kdpnorm,kdpmap,x_lims,-0.3,0.5,count_threshold=count_thr, ax=ax_arr[3])
	qvpp.annotate_field_plot(r'$\mathrm{K_{DP}} (^\circ km^{-1})$',kdpnorm,hmax,date_format,timeticks, ax=ax_arr[3], im=im)
	#qvpp.z_contourplt(qvp_file, count_thr, x_lims, ax=ax_arr[3])
	
	
	print "subplot(514)"
	
	## Plot uPhiDP - Subplot 5 ----------------------------------------------------
	#plt.subplot(515)
	#print ver_file
	hs.plot(ax=ax_arr[4],color = '0.45', alpha=0.5, lw=3,label='FAAM flight')
	im = qvpp.plot_field(ver_file,'Vu',Vnorm,Vmap,vx_lims,-0.5,8,scaler=-1,count_threshold=count_thr, ax=ax_arr[4])
	qvpp.annotate_field_plot('Fall speed (m/s)',Vnorm,hmax,date_format,timeticks, ax=ax_arr[4], im=im)
	#qvpp.z_contourplt(qvp_file, count_thr, x_lims, ax=ax_arr[4])
	
	
	print "subplot(515)"
	
	plt.tight_layout()
	
	if savefile:
		savefile2 = output_dir + "qvp_5subplots.png"
		plt.savefig(savefile2 , format='png',dpi=200)
	#plt.show()
	plt.close()
	
	print "%sqvp_5subplots.png is saved" %output_dir
	# In[18]:
	
	# In[19]:
	
	
	timeticks = [nc.num2date((qvp_file['Time'][0]-qvp_file['Time'][0]%3600)+x,
	                      qvp_file['Time'].units) for x in range(0,int(qvp_file['Time'][-1]-qvp_file['Time'][0]),3600)]
	
	print "timeticks are ready"
	
	# In[20]:
	if(zoom_start is not None):X0 = zoom_start.hour
	else: X0 = faam_start_time.hour
	if(zoom_end is not None):X1 = zoom_end.hour
	else:X1 = faam_end_time.hour
	print '%d, %d' %(X0,X1)
	qvpp.Plot_standard_6(qvp_file,ver_file,timeticks,output_dir,faam_file,savefile=savefile,x0=X0,x1=X1,count_thr=100,hmin=0,hmax = max_height,hs=hs)
	print "Plot_standard_6 is ready"
	
	
	
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
        main()
	
