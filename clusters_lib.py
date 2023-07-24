#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  clusters_lib.py
#  
#  Copyright 2018 Maryna Lukach <earmlu@SEE-L11936.local>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import numpy as np
import h5py

from scipy.spatial import distance
from scipy import ndimage
import sklearn.preprocessing as prp
import sklearn.utils as utl
import functools
from sklearn.metrics import pairwise_distances_chunked
from sklearn.preprocessing.data import QuantileTransformer

from sklearn.decomposition import PCA
#, FastICA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from sklearn import cluster

from pandas.plotting import scatter_matrix

from sklearn.metrics import silhouette_samples, silhouette_score, pairwise_distances_argmin_min,calinski_harabaz_score,davies_bouldin_score
from anytree import Node, RenderTree
from anytree.render import ContStyle
#from anytree.exporter import DotExporter
from anytree.dotexport import RenderTreeGraph
import matplotlib.pyplot as plt

import QVP_plot_functions as qvpp
plot_kwds = {'alpha' : 0.5, 's' : 10, 'linewidths':0}

# added for prints and tests 

import pandas as pd

#import skfuzzy as fuzz
import seaborn as sns

from sklearn.feature_extraction import image
from sklearn.cluster import spectral_clustering
import warnings
from sklearn.neighbors.nearest_centroid import NearestCentroid

from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh

from numpy.linalg import norm
import datetime as dt
import pytz
from astral import Astral

#from rpy2.robjects.packages import importr
#import rpy2.robjects as ro
#import pandas.rpy.common as com

class clusternode(object):
	
	def __init__(self,name="root",parent=None,n_points=0,data=None,pca=None,parent_pca=None,variables=[],medoid=[],center=None,indxs=None,n_subclusters=0,subclusters=None):
		#,root_pca=None
		self.name = name
		self.parent = parent
		self.n_points = n_points
		self.variables = variables
		self.data = data
		self.pca = pca
		self.parent_pca = parent_pca
		self.medoid = medoid
		self.center = center
		
		#if((indxs is None)&(name is "root")):self.indxs = numpy.ones((self.n_points), dtype=bool)#range(self.n_points)
		#else:
		self.indxs = indxs
		self.n_subclusters = n_subclusters
		self.subclusters = []
		if subclusters is not None:
			for subcluster in subclusters:
				self.add_subnode(subcluster)
	def __repr__(self):
		return self
	def add_subnode(self, clusternode):
		assert isinstance(name,clusternode)
		self.subclusters.append(name)
	def get_children(self):
		return self.subclusters
	def get_parent():
		return self.parent

class spectr_model(object):
	def __init__(self,labels,cntr):
		self.labels_=labels
		self.cluster_centers_=cntr
		self.colors=[]
		self.labels_list=[]
	def fit(self,X):
		return self.X
	def add_list(self,colors,labels_list):
		self.colors=colors
		self.labels_list=labels_list

class fpc_model(object):
	def __init__(self,cntr, u, u0, d, jm, p):
		self.X=np.transpose(u)
		self.labels_=np.argmax(u,axis=0)
		self.cluster_centers_=cntr
		self.colors=[]
		self.labels_list=[]
	def fit(self,X):
		return self.X
	def add_list(self,colors,labels_list):
		self.colors=colors
		self.labels_list=labels_list
		
# Input    data : 2d array, size (S, N) Data to be clustered. N is the number of data sets; S is the number of features within each sample vector.
# Output			
#cntr : 2d array, size (S, c)# Cluster centers. Data for each center along each feature provided for every cluster (of the c requested clusters).

# u : 2d array, (S, N) # Final fuzzy c-partitioned matrix.

# u0 : 2d array, (S, N) # Initial guess at fuzzy c-partitioned matrix (either provided init or random guess used if init was not provided).

# d : 2d array, (S, N) # Final Euclidian distance matrix.

# jm : 1d array, length P # Objective function history.

# p : int # Number of iterations run.

	# for model in algorithms.values():
		# if debug:
			# print "algorithms.values()"
			# print model
		# #model.fit(X_new[:, 0:2]) -  for ERAD plots it was like that
		# model.fit(X_new)
		# results.append(list(model.labels_))
		# algorithms['kmeans'].fit(X_new)#[:, 0:2])
		# if debug:
			# print  model.labels_
			# print model.inertia_
			# print "model.cluster_centers_"
			# print model.cluster_centers_

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


			
def get_night_hours_for_one_day(timeserie):
	
	a = Astral()
	a.solar_depression = "civil"
	
	#set city to Lnodon and change coordinates to Porton Down
	city = a["London"]
	# providing coordinates of Porton Down,  UK
	city.latitude = float(51.1167)
	city.longitude = float(1.7000)
	utc=pytz.UTC
	# get the date of the case/day to select sunrise and the sunset
	print "In get_night_hours_foroneday"
	#print "timeserie"
	#print timeserie
	#print "timeserie[0]"
	#print timeserie[0]
	case_day = timeserie[0]
	next_day = timeserie[0]+dt.timedelta(days=1)
	
	civil_twilight_this_day = city.sun(date=case_day, local=False)
	civil_twilight_next_day = city.sun(date=next_day, local=False)
	#print timeserie > civil_twilight_this_day['sunset'].replace(tzinfo=None)

	#night_hours_mask = (timeserie > civil_twilight_this_day['dusk'].replace(tzinfo=None)) & (timeserie < civil_twilight_next_day['dawn'].replace(tzinfo=None))
	#morning_hours_mask = (timeserie< civil_twilight_this_day['dawn'].replace(tzinfo=None))
	night_hours_mask = (timeserie > civil_twilight_this_day['sunset'].replace(tzinfo=None)) & (timeserie < civil_twilight_next_day['sunrise'].replace(tzinfo=None))
	morning_hours_mask = (timeserie < civil_twilight_this_day['sunrise'].replace(tzinfo=None))
	
	
	#night_hours_mask = np.tile(night_hours_mask, (alt_shape,1))
	#morning_hours_mask  = np.tile(morning_hours_mask, (alt_shape,1))
	
	hours_mask = np.logical_or(night_hours_mask, morning_hours_mask)
	#print "np.where(hours_mask == True)"
	#print "np.where(hours_mask == True)"
	#print np.where(hours_mask == True)
	
	#return (hours_mask,night_hours_mask,civil_twilight_this_day['dusk'].replace(tzinfo=None),civil_twilight_next_day['dawn'].replace(tzinfo=None))
	return (hours_mask,night_hours_mask,civil_twilight_this_day['sunset'].replace(tzinfo=None),civil_twilight_next_day['sunrise'].replace(tzinfo=None))
	


	
def plotting_tree(clust,parent=None,root_node=None, filename="test.png"):	
	"""
	plotting the tree starting from the provided cluster - clust
	prints the rendertree  and exports it to the graph
	
	the function runs in parallel recursive calls
	"""

	if(parent==None):cl_fullname = clust.name
	else:
		cl_fullname_list = get_cluster_parent_names(clust)
		cl_fullname = ".".join(cl_fullname_list)
		
	node = Node(cl_fullname, parent=parent)
	if(root_node==None):root_node = node
	if(len(clust.subclusters)>0):
		for j in range(clust.n_subclusters):
			plotting_tree(clust.subclusters[j], parent=node,root_node=root_node,filename=filename)
			
	if(parent==None):
		print(RenderTree(root_node, style=ContStyle()))
		print "tree filename is {}".format(filename)
		
		try:
			RenderTreeGraph(root_node).to_picture(filename)
		except :
			print '\n*** INFO: exception occured during printing the tree'
			print ' - **EXITING**'


def getAffinityMatrix(coordinates, k = 7):
    """
    Calculate affinity matrix based on input coordinates matrix and the numeber
    of nearest neighbours.
    
    Apply local scaling based on the k nearest neighbour
        References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    """
    # calculate euclidian distance matrix
    dists = distance.squareform(distance.pdist(coordinates)) 
    
    # for each row, sort the distances ascendingly and take the index of the 
    #k-th position (nearest neighbour)
    knn_distances = np.sort(dists, axis=0)[k]
    knn_distances = knn_distances[np.newaxis].T
    
    # calculate sigma_i * sigma_j
    local_scale = knn_distances.dot(knn_distances.T)

    affinity_matrix = dists * dists
    affinity_matrix = -affinity_matrix / local_scale
    # divide square distance matrix by local scale
    affinity_matrix[np.where(np.isnan(affinity_matrix))] = 0.0
    # apply exponential
    affinity_matrix = np.exp(affinity_matrix)
    np.fill_diagonal(affinity_matrix, 0)
    return affinity_matrix

def get_number_of_points_per_cluster(number_of_subclusters,labels):
	n_points = []
	for i in range(number_of_subclusters):
		n_points.append(len(labels[labels==i]))
	return n_points


def write_cluster_group(cluster, parent_group):
	cluster_group = parent_group.create_group(cluster.name)
	cluster_group.attrs["n_points"] = cluster.n_points
	cluster_group.attrs["n_subclusters"] = cluster.n_subclusters
	data_id = cluster_group.create_dataset("data", (cluster.data.shape), data=cluster.data)
	if(cluster.pca is not None):pca_id = cluster_group.create_dataset("pca", data=cluster.pca)#(cluster.pca.shape), 
	if(cluster.medoid is not None):medoid_id = cluster_group.create_dataset("medoid", (np.asarray(cluster.medoid).shape), data=np.asarray(cluster.medoid))
	if(cluster.center is not None):center_id = cluster_group.create_dataset("center", (cluster.center.shape), data=cluster.center)
	if(cluster.parent_pca is not None):parent_pca_id = cluster_group.create_dataset("parent_pca", (cluster.parent_pca.shape), data=cluster.parent_pca)
	if(len(cluster.variables) is not 0):
		#print "(len(cluster.variables) is not 0) and is {} for the cluster {}".format(len(cluster.variables),cluster.name)
		variables_id = cluster_group.create_dataset("variables", (np.asarray(cluster.variables).shape), data=np.asarray(cluster.variables))
		#print variables_id
	#else:print "len(cluster.variables) is  {} nad {} ".format(len(cluster.variables),cluster.variables)
	indxs_id = cluster_group.create_dataset("indxs", (cluster.indxs.shape), data=cluster.indxs)
	
	return cluster_group
			
def write_all_subclusters(cluster, parent_group):
	if(cluster.n_subclusters>0):
		for subcluster in cluster.subclusters:
			group_id = write_cluster_group(subcluster, parent_group) 
			write_all_subclusters(subcluster, group_id)
	
def write_clusters(root, output_dir,input_filename=None,verbose=True):
	"""
	TODO: will save the tree to the HDF5-file
	
	"""
	childrens_list = get_active_clusteres_list(root,childrens_list=[])
	#print childrens_list
	if input_filename is not None:output_filename = "{}{}_{}_data_{}_active_clusters_v1.hdf".format(output_dir,input_filename,root.name,len(childrens_list))
	else:output_filename = "{}{}_data_{}_active_clusters_v1.hdf".format(output_dir,root.name,len(childrens_list))
	f = h5py.File(output_filename, "w")
	parent_group = write_cluster_group(root, f)
	write_all_subclusters(root, parent_group)
	if verbose:print "writing the {} file is finished".format(output_filename)
	f.close()

def read_group_into_cluster(group,parent=None):
	if(parent == None):
		cluster = clusternode(name=group.name.split("/",-1)[-1],n_points=group.attrs['n_points'],
							data=group['data'],indxs=group['indxs'],
							variables=group['variables'][:],
							n_subclusters = group.attrs['n_subclusters'])
	else:
		cluster = clusternode(name=group.name.split("/",-1)[-1],n_points=group.attrs['n_points'],
							data=group['data'],indxs=group['indxs'],
							variables=group['variables'][:],
							n_subclusters = group.attrs['n_subclusters'],
							pca=group['pca'],medoid=group['medoid'],parent=parent)
	
	for cl in [x for x in group.keys() if x.startswith('cl')]:	
		cl_id = read_group_into_cluster(group[cl],parent=cluster)
		cluster.subclusters.append(cl_id)
	return cluster
	
	

def read_clusters(file_in):
	"""
	TODO: reads tree from the HDF5 file
	
	
	"""
	print "in read_clusters, reading file {}".format(file_in)
	f = h5py.File(file_in, "r")
	print f
	Datasetnames=f.keys()
	print Datasetnames
	root = None
	for group in Datasetnames:
		root = read_group_into_cluster(f[group],parent=None)
	return root
    
	
def get_active_clusteres_data(clust,childrens_list=[],data=[],labels=[],medoids=[],indxs=[],debug=False):
	"""
	recursive function
	gets cluster id
	returns list of active clusters id's, full data, cluster labels and cluster medoids if they were calculated
	"""
	if debug:
		print "in get_active_clusteres_data"
		print "function was called for clust {}, which has {} subclusters".format(clust.name,len(clust.subclusters))
		print "labels {}; len(childrens_list) {}".format(labels,len(childrens_list))
	
	if(len(clust.subclusters)==0):
		childrens_list.append(clust)
		if debug:
			print "this is an active cluster {}".format(clust.name)
			print "len(childrens_list) is {}".format(len(childrens_list))
		
		
		if(len(data)>0):data = np.append(np.asarray(data),np.asarray(clust.data),axis=0)
		else:data = np.asarray(clust.data)
		if(len(labels)>0):labels = np.append(np.asarray(labels),np.full((clust.data.shape[0]),(len(childrens_list)-1)),axis=0)
		else:labels = np.full((clust.data.shape[0]),(len(childrens_list)-1))
		
		if(len(medoids) != 0):
			
			medoids = np.append(np.asarray(medoids),np.asarray([clust.medoid]),axis=0)
			
		else:
			medoids = np.asarray([clust.medoid])
			
		if(len(indxs)>0):indxs = np.append(np.asarray(indxs),np.asarray(clust.indxs),axis=0)
		else:indxs = np.asarray(clust.indxs)
		
	else:
		if debug:print "clust {} has {} subclusters".format(clust.name,len(clust.subclusters))
		
		for subcluster in clust.subclusters:
			if debug:
				print "working with subclaster {}".format(subcluster.name)
				print "labels {}; len(childrens_list) {}".format(labels,len(childrens_list))
				for child in childrens_list:
					print child.name
				print "before adding np.unique(labels) = {} and len(labels) = {};len(indxs) = {}".format(np.unique(labels),len(labels),len(indxs))
			(childrens_list,data,labels,medoids,indxs) = get_active_clusteres_data(subcluster,childrens_list=childrens_list,data=data,labels=labels,medoids=medoids,indxs=indxs)
			if debug:
				print "after adding np.unique(labels) = {} and len(labels) = {};len(indxs) = {}".format(np.unique(labels),len(labels),len(indxs))
				print "len(childrens_list) {}".format(len(childrens_list))
	#print "get_active_clusteres_data"
	return (childrens_list,data,labels,medoids,indxs)
	 
def sort_active_clusters_size_based(X_new,labels, debug=True):
	"""
	gets cluster data and cluster labels
	returns list of cluster labels, and sorted list of cluster nodes
	sort is based on the size of the cluster (number of points in the cluster)
	"""
	cluster_size = []
	cluster_indexes = []
	if debug:
		print "in sort_active_clusters_size_based"
		print np.unique(labels)
	
	for k in np.unique(labels):
		cluster_indexes.append(np.where(labels == k)[0])
		cluster_size.append(len(labels[labels == k]))
	sorted_list_of_clusteres = sorted(range(len(cluster_size)), key=lambda k: cluster_size[k], reverse=True)
	#if debug:
	print "cluster_size is {} and sorted_list_of_clusteres {}".format(cluster_size, sorted_list_of_clusteres)
	return (cluster_indexes,sorted_list_of_clusteres)
		
def sort_active_clusters_decending(X_new,labels):
	"""
	gets cluster data and cluster labels
	returns list of cluster labels, and sorted list of cluster nodes
	"""
	inter_cluster_max_dist = []
	cluster_indexes = []
	#print "in sort_active_clusters_decending"
	#print np.unique(labels)
	for k in np.unique(labels):
		#print k
		cluster_indexes.append(np.where(labels == k)[0])
		#cluster_dist = distance.cdist(X_new[cluster_indexes[k],:], X_new[cluster_indexes[k],:])
		cluster_dist = distance.cdist(X_new[cluster_indexes[k],:], X_new[cluster_indexes[k],:])**2
		inter_cluster_max_dist.append(np.amax(cluster_dist))
	#print inter_cluster_max_dist
	
	sorted_list_of_clusteres = sorted(range(len(inter_cluster_max_dist)), key=lambda k: inter_cluster_max_dist[k], reverse=True)
	#print sorted_list_of_clusteres
	return (cluster_indexes,sorted_list_of_clusteres)
	
def sort_clusters_decending(model,X_new,selected_pca,verbose=False):
	"""
	gets K-mean model object, PCA-transformed data and the number of selected PCA components
	returns list of cluster labels, and sorted list of cluster labels
	"""
	inter_cluster_max_dist = []
	cluster_indexes = []
	number_of_clusters = len(np.unique(model.labels_))
	for k in range(number_of_clusters):
		cluster_indexes.append(np.where(model.labels_ == k)[0])
		cluster_dist = distance.cdist(X_new[cluster_indexes[k],0:selected_pca], X_new[cluster_indexes[k],0:selected_pca])
		inter_cluster_max_dist.append(np.amax(cluster_dist))
	
	sorted_list_of_clusteres = sorted(range(len(inter_cluster_max_dist)), key=lambda k: inter_cluster_max_dist[k], reverse=True)
	if verbose:
		print "inter_cluster_max_dist"
		print inter_cluster_max_dist
		print sorted_list_of_clusteres
	return (cluster_indexes,sorted_list_of_clusteres)	

def get_subclusters(data,columns=['dBZ','ZDR','RhoHV','PhiDP','temperature_2'],verbose=False):
	""" gets an optimal number of clusters in the data
		start with 2 clusters, increases the number of clusters while the optimal number is reached
	    input:		data
	    output:		number_of_clusters (int), PCA-transformed dataset, number of selected PCA, K-mean model object (sklearn)
	"""
	if verbose:print "in get_subclusters, the input data.shape is {}".format(data.shape)
	number_of_clusters = 2
	sil_diff = 0
	sil = []
	model_list = []
	# Get transformed data
	X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
	#print X_std.shape
	#print columns
	dataframe_X_std = pd.DataFrame(data=X_std, columns=columns)
	
	if verbose:
		print "describe X_std"
		print dataframe_X_std.describe()
	
	# Get PCA for X_std
	(selected_pca, X_new) = get_PCA(X_std)
	new_field_list = ["X_{}".format(i) for i in range(len(columns))]#["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	pd_X_new = pd.DataFrame(data = X_new, columns = new_field_list[0:len(columns)])
	if verbose:
		print ("selected_pca")
		print (selected_pca)
		print ("describe X_new")
		print (pd_X_new.describe())
		print ("X_new.shape")
		print (X_new.shape)
		print (X_new[:,0:selected_pca].shape)
	
	while (sil_diff >= 0):
		# Get clusters
		(model) = get_k_mean(X_new[:, 0:selected_pca],k_clusters=number_of_clusters)
		sil.append(silhouette_score(X_new[:, 0:selected_pca], model.labels_))
		model_list.append(model)
		if(number_of_clusters > 2):
			sil_diff = sil[-1] - sil[-2]
			#print "number_of_clusters is {}".format(number_of_clusters) 
			#print "sil[-1] is {} and sil[-2] is {}. The sil_diff is {}  ".format(sil[-1],sil[-2],sil_diff)
		number_of_clusters += 1
	
	return(number_of_clusters-2,X_new,selected_pca,model_list[-2])

# def get_subclusters_with_fuzz(data,columns=['dBZ','ZDR','RhoHV','PhiDP','temperature_2'], threshold=0.5,plot=True,verbose=False,debug=False):
	# """ gets an optimal number of clusters in the data
		# start with 2 clusters, increases the number of clusters while the optimal number is reached
	    # input:		data
	    # output:		number_of_clusters (int), PCA-transformed dataset, number of selected PCA, K-mean model object (sklearn)
	# """
	# #if verbose:
	# print "in get_subclusters_with_fuzz, the input data.shape is {} and the list of variables is {}".format(data.shape,columns)
	# output_dir = "./data/"
	# number_of_clusters = 2
	# fpc_diff = 0
	# fpc = []
	# sil_diff = 0
	# sil = []
	# model_list = []
	# # Get transformed data
	# X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
	# dataframe_X = pd.DataFrame(data=data, columns=columns)
	# # if(plot):
		# # scatter_matrix(dataframe_X,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
		# # savefile2 = output_dir + "original_scatters_pd_X_{}points.png".format(len(data[:,0]))
		# # #qvpp.pca_scatter(data,columns,savefile2,plot_kwds=plot_kwds)
		# # plt.savefig(savefile2 , format='png',dpi=200)
		# # plt.close()
	# #print X_std.shape
	# #print columns
	# dataframe_X_std = pd.DataFrame(data=X_std, columns=columns)
	# # if(plot):
		# # scatter_matrix(dataframe_X_std,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
		# # savefile2 = output_dir + "transformed_scatters_pd_X_std_{}points.png".format(len(data[:,0]))
		# # #qvpp.pca_scatter(dataframe_X_std,columns,savefile2,plot_kwds=plot_kwds)
		# # plt.savefig(savefile2 , format='png',dpi=200)
		# # plt.close()
	# if verbose:
		# print ("describe X_std")
		# print (dataframe_X_std.describe())
	
	# # Get PCA for X_std
	# (selected_pca, X_new) = get_PCA(X_std,evr_threshold=0.70,debug=False)
	# new_field_list = ["X_{}".format(i) for i in range(len(columns))]#["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	# # Get ICA for X_std
	# #X_new = get_ICA(X_new[:, 0:selected_pca],evr_threshold=0.70,n_components=selected_pca,debug=True)
	
	# pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
	
	# # if(plot):
		# # scatter_matrix(pd_X_new,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
		# # savefile2 = output_dir + "PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		# # plt.savefig(savefile2 , format='png',dpi=200)
		# # plt.close()
		
	# """
	# # Get ICA for X_std
	# ( X_new_ica) = get_ICA(X_new,evr_threshold=0.70,n_components=selected_pca,debug=True)
	# #new_field_list = ["X_{}".format(i) for i in range(len(columns))]#["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	# pd_X_new_ica = pd.DataFrame(data = X_new_ica, columns = new_field_list[0:selected_pca])
	
	# if(plot):
		 
		# scatter_matrix(pd_X_new_ica,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
		# savefile2 = output_dir + "ICA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		# plt.savefig(savefile2 , format='png',dpi=200)
		# plt.close()	
		
	# if verbose:
		# print ("selected_pca")
		# print (selected_pca)
		# print ("describe X_new")
		# print (pd_X_new.describe())
		# print ("X_new.shape {}; {}".format(X_new.shape, X_new[:,0:selected_pca].shape))
		# print (pd_X_new_ica.describe())
		# print ("X_new_ica.shape {}".format(pd_X_new_ica.shape))
	# """                                                                           
	
	# #selected_pca = len(pd_X_new_ica[0,:])
	# #while (fpc_diff >= 0):
	# while (sil_diff >= 0):
		# # Get clusters
		# #print (X_new[:, 0:selected_pca].shape)
		# #(model) = get_k_mean(X_new[:, 0:selected_pca],k_clusters=number_of_clusters)
		# #print ("model.labels_ {}".format(model.labels_.shape))
		# #print ("model.cluster_centers_ {}".format(model.cluster_centers_.shape))
		# #print (model.labels_[0])
		# #print (np.unique(model.labels_))
		# #data : 2d array, size (S, N) Data to be clustered. N is the number of data sets; S is the number of features within each sample vector.
		
		# # get fuzzy clustering results n_init=200)#, max_iter=500,tol=0.0001)
		# cntr, u, u0, d, jm, p, fpc_this_one  = fuzz.cmeans(np.transpose(X_new[:, 0:selected_pca]), number_of_clusters, 1.4, error=0.00001, maxiter=10000, init=None)
		# #cntr, u, u0, d, jm, p, fpc_this_one  = fuzz.cmeans(np.transpose(X_new_ica), number_of_clusters, 2, error=0.0001, maxiter=1000, init=None)
		
		# fpc.append(fpc_this_one)
		# print ("input variables {}".format(columns))
		# print("fpc {}".format(fpc))
		# # create model object for skfuzzy.cmeans results
		# fuzz_model = fpc_model(cntr, u, u0, d, jm, p)
		# print ("fuzz_model.labels_ {}".format(fuzz_model.labels_.shape))
		# print ("fuzz_model.cluster_centers_ {}".format(fuzz_model.cluster_centers_.shape))
		# # print ("u.shape {}".format(u.shape))
		# # print (u[:,0])
		# # print (fuzz_model.labels_[0])
		# # print (np.unique(fuzz_model.labels_))
		# print ("len of array with max score >= {}".format(threshold))
		# n_well_assigned = len(fuzz_model.labels_[np.max(u,axis=0)>=threshold])
		# print (n_well_assigned)
		
		# new_lables = np.copy(fuzz_model.labels_)
		# new_lables[np.max(u,axis=0)<threshold] = -1
		# #print ("new_lables")
		# #print (new_lables)
		# #print (new_lables[new_lables<0])
		# n_poor_assigned = len(new_lables[new_lables<0])
		
		# print ("len(new_lables[new_lables==-1])")
		# print (n_poor_assigned)
		
		# model_list.append(fuzz_model)
		
		# print("!!! calculating the silhouette index !!!")
		# f_sil = get_fuzz_SIL(X_new[:, 0:selected_pca], u)
		
		# #sil.append(silhouette_score(X_new[:, 0:selected_pca], fuzz_model.labels_))
		
		# sil.append(f_sil)
		# #print("crisp sil index is {}".format(sil[-1]))
		# #if(number_of_clusters > 8):
		# #	fpc_diff = fpc[-1] - fpc[-2]
		# if(number_of_clusters > 2):
			# sil_diff = sil[-1] - sil[-2]
			# if(n_poor_assigned<n_well_assigned ):
				# sil_diff = sil[-1] - sil[-2]
			# else:
				# print ("break due to n_poor_assigned<n_well_assigned number_of_clusters now is {}".format(number_of_clusters))
				# sil_diff = -1
			# #print "number_of_clusters is {}".format(number_of_clusters) 
			# #print "sil[-1] is {} and sil[-2] is {}. The sil_diff is {}  ".format(sil[-1],sil[-2],sil_diff)
		# number_of_clusters += 1
		
	# if(n_poor_assigned<n_well_assigned ):
		# pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
		# #pd_X_new = pd.DataFrame(data = X_new_ica, columns = new_field_list[0:selected_pca])
		# new_labels = np.copy(model_list[-2].labels_)
		# new_labels[np.max(model_list[-2].X,axis=1)<threshold] = -1
		
		# pd_X_new['Clusters'] = pd.Series(new_labels)
		
		# palette = sns.color_palette('deep',n_colors=(len(np.unique(new_labels))-1))
		# #print ("palette is {}".format(palette))
		
		# colors = [palette[z] if z >= 0 else (0.8, 0.8, 0.8) for z in np.unique(new_labels)]
		# print ("pd_X_new['Clusters']  are {}".format(np.unique(pd_X_new['Clusters'])))
		# print ("np.unique(new_lables) {}".format(np.unique(new_labels)))
		
		# #figurename = output_dir + "c_clustered&prob_PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		# figurename = output_dir + "c_clustered&prob_ICA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		
		# print ("plotting clusters with -1 cluster with np.unique(new_labels) {}".format(np.unique(new_labels)))
		# #(colors,labels_list) = qvpp.clustered_scatter(pd_X_new,figurename,labels_list=[],plot_kwds=plot_kwds)
		
		# labels_list = []
		# for i in np.unique(new_labels):
			# if(i==-1):
				# labels_list.append("---")
				# pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': labels_list[-1]}, regex=False)
			# else:
				# labels_list.append('cl{}'.format(i+1))
				# pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': labels_list[-1]}, regex=False)
		
		# model_list[-2].add_list(colors,labels_list)
		
		# #pd_X_new = pd_X_new.replace({'Clusters': -1}, {'Clusters': '---'}, regex=False)
		# #pd_X_new = pd_X_new.replace({'Clusters': 0}, {'Clusters': 'cl1'}, regex=False)
		# #pd_X_new = pd_X_new.replace({'Clusters': 1}, {'Clusters': 'cl2'}, regex=False)
		# # print ("new_lables_list is {}".format(new_labels_list))
		# if(plot):
			# #vars=new_field_list[0:selected_pca],  scatter_matrix(pd_X_new,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
			# g = sns.pairplot(pd_X_new,diag_kind="hist", hue='Clusters', hue_order=labels_list, palette=colors, markers=".",diag_kws={'alpha':0.7},plot_kws=plot_kwds)#'deep')#"husl")
			# savefile2 = output_dir + "c_clustered&prob_PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
			# plt.savefig(savefile2 , format='png',dpi=200)
			# plt.close()
		
		# #print ("number_of_clusters now is {} so, the returned value is {};model_list {} and we return {} with cluster_centers_ {}".format(number_of_clusters,number_of_clusters-2,model_list,model_list[-2],model_list[-2].cluster_centers_))
		
		# pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
		# pd_X_new['Clusters'] = pd.Series(model_list[-2].labels_)
		
		# print labels_list[1:]
		# figurename = output_dir + "c_clustered_PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		# (colors,labels_list) = qvpp.clustered_scatter(pd_X_new,figurename,colors=colors[1:],labels_list=labels_list[1:],plot_kwds=plot_kwds)
		
		# #for i in np.unique(model_list[-2].labels_):
		# #	pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': 'cl{}'.format(i)}, regex=False)
		# #pd_X_new = pd_X_new.replace({'Clusters': 0}, {'Clusters': 'cl1'}, regex=False)
		# #pd_X_new = pd_X_new.replace({'Clusters': 1}, {'Clusters': 'cl2'}, regex=False)
		
		# #if(plot):
			# #vars=new_field_list[0:selected_pca] scatter_matrix(pd_X_new,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
		# #	k = sns.pairplot(pd_X_new, diag_kind="hist", hue='Clusters', hue_order=new_labels_list[1:], palette=colors[1:], markers=".",diag_kws={'alpha':0.3},plot_kws=plot_kwds)#"husl")
		# #	savefile2 = output_dir + "c_clustered_PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
		# #	plt.savefig(savefile2 , format='png',dpi=200)
		# #	plt.close()
		
		# #print ("number_of_clusters now is {} so, the returned value is {};model_list {} and we return {} with cluster_centers_ {}".format(number_of_clusters,number_of_clusters-2,model_list,model_list[-2],model_list[-2].cluster_centers_))
		
		
		# X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
		# dataframe_X = pd.DataFrame(data=data, columns=columns)
		# dataframe_X['Clusters'] = pd.Series(model_list[-2].labels_)#new_labels)
		
		# figurename = output_dir + "c_clustered_original_scatters_pd_X_{}points.png".format(len(data[:,0]))
		# qvpp.clustered_scatter(dataframe_X,figurename,colors=model_list[-2].colors[1:],labels_list=model_list[-2].labels_list[1:],plot_kwds=plot_kwds)
		
		# #dataframe_X = dataframe_X.replace({'Clusters': -1}, {'Clusters': '---'}, regex=False)
		# #dataframe_X = dataframe_X.replace({'Clusters': 0}, {'Clusters': 'cl1'}, regex=False)
		# #dataframe_X = dataframe_X.replace({'Clusters': 1}, {'Clusters': 'cl2'}, regex=False)
		# # for i in np.unique(model_list[-2].labels_):
			# # if(i==-1):
				# # dataframe_X  = dataframe_X .replace({'Clusters': i}, {'Clusters': 'cl{}'.format(i+1)}, regex=False)
			# # else:
				# # dataframe_X  = dataframe_X .replace({'Clusters': i}, {'Clusters': "cl{}".format(i+1)}, regex=False)
		# # print("clusters {} labels {} colors {}  new_labels_list {}".format(np.unique(dataframe_X['Clusters']),np.unique(model_list[-2].labels_),model_list[-2].colors,model_list[-2].labels_list))
		# # if(plot):
			# # #vars=columns, scatter_matrix(dataframe_X,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
			# # l = sns.pairplot(dataframe_X, hue='Clusters', hue_order=model_list[-2].labels_list[1:], palette=model_list[-2].colors[1:],diag_kind="hist", markers=".",diag_kws={'alpha':0.3},plot_kws=plot_kwds)#"husl")
			# # savefile2 = output_dir + "c_clustered_original_scatters_pd_X_{}points.png".format(len(data[:,0]))
			# # #qvpp.pca_scatter(data,columns,savefile2,plot_kwds=plot_kwds)
			# # plt.savefig(savefile2 , format='png',dpi=200)
			# # plt.close()
			
	# return(number_of_clusters-2,X_new,selected_pca,model_list[-2])

def eigenDecomposition(A, plot = True):
    """
    :param A: Affinity matrix
    :param plot: plots the sorted eigen values for visual inspection
    :return A tuple containing:
    - the optimal number of clusters by eigengap heuristic
    - all eigen values
    - all eigen vectors
    
    This method performs the eigen decomposition on a given affinity matrix,
    following the steps recommended in the paper:
    1. Construct the normalized affinity matrix: L = D−1/2ADˆ −1/2.
    2. Find the eigenvalues and their associated eigen vectors
    3. Identify the maximum gap which corresponds to the number of clusters
    by eigengap heuristic
    
    References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf
    """
    L = csgraph.laplacian(A, normed=True)
    n_components = A.shape[0]
    
    # LM parameter : Eigenvalues with largest magnitude (eigs, eigsh), that is, largest eigenvalues in 
    # the euclidean norm of complex numbers.
    eigenvalues, eigenvectors = eigsh(L, k=n_components, which="LM", sigma=1.0, maxiter=5000)
    
    if plot:
        plt.title('Largest eigen values of input matrix')
        plt.scatter(np.arange(len(eigenvalues)), eigenvalues)
        plt.grid()
        
    # Identify the optimal number of clusters as the index corresponding
    # to the larger gap between eigen values
    index_largest_gap = np.argmax(np.diff(eigenvalues))
    nb_clusters = index_largest_gap + 1
        
    return nb_clusters, eigenvalues, eigenvectors

def vecas(data, lables, centroids):
	total_mean = np.mean(data)
	clusters = np.unique(lables)
	c = len(clusters)
	for i in clusters:
		n_i = len(labels[lables==i])
		norm2_x = norm(data[lables==1] - total_mean,2)
		beta_comp = norm2_x/n_i
		exp_comp = np.exp()
		
		
		                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
	return vecas


def get_subclusters_with_spectr(data,columns=['dBZ','ZDR','RhoHV','PhiDP','DOP'], threshold=0.5, cluster_name=None, output_dir = "./data/", plot=True,verbose=False,debug=False):
	""" gets an optimal number of clusters in the data
		start with 2 clusters, increases the number of clusters while the optimal number is reached
	    input:		data
	    output:		number_of_clusters (int), PCA-transformed dataset, number of selected PCA, K-mean model object (sklearn)
	"""
	#if verbose:
	print "in get_subclusters_with_spectr, the input data.shape is {} and the list of variables is {}".format(data.shape,columns)
	
	number_of_clusters = 2
	chs_diff = 0
	chs = []
	# Silhouette score
	sil_diff = 0
	sil = []
	# Davies-Bouldin score
	dbs_diff = 0
	dbs = []
	
	WG_diff = 0
	WG = []
	
	model_list = []
	# Get transformed data
	X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
	dataframe_X = pd.DataFrame(data=data, columns=columns)
	
	dataframe_X_std = pd.DataFrame(data=X_std, columns=columns)
	
	if verbose:
		print ("describe X_std")
		print (dataframe_X_std.describe())
	
	# Get PCA for X_std
	(selected_pca, X_new) = get_PCA(X_std,evr_threshold=0.70,debug=False)
	#X_new = get_ICA(X_std,evr_threshold=0.70,n_components=selected_pca,debug=True)
	#get_PCA(X_std,evr_threshold=0.70,debug=True)
	
	#lda = LinearDiscriminantAnalysis(n_components=selected_pca,solver='lsqr', shrinkage='auto')
	#X_r2 = lda.fit(X_new, y).transform(X_new)
	
	
	new_field_list = ["X_{}".format(i) for i in range(len(columns))]#["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	# Get ICA for X_std
	#X_new = get_ICA(X_new[:, 0:selected_pca],evr_threshold=0.70,n_components=selected_pca,debug=True)
	
	pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
		
	#affinity_matrix = getAffinityMatrix(X_new[:, 0:selected_pca], k = 7)
	#print "affinity_matrix.shape {}".format(affinity_matrix.shape)
	#number_of_clusters, _,  _ = eigenDecomposition(affinity_matrix)
	#print "eigenDecomposition number_of_clusters is {}".format(number_of_clusters)
	
	total_cond = 0 
	while ((sil_diff >= 0)):# & (X_new[:, 0:selected_pca].shape[0] > number_of_clusters)):
	#while not WG_diff:
		print "number_of_clusters {}".format(number_of_clusters)
	                
	#while ((total_cond < 2) ):#& (X_new[:, 0:selected_pca].shape[0] > number_of_clusters)):
		#print "number_of_clusters is {}".format(number_of_clusters)
		#print "X_new[:, 0:selected_pca].shape {}".format(X_new[:, 0:selected_pca].shape)
	#while (chs_diff >= 0):
		# Get clusters
		#spectral = cluster.SpectralClustering(n_clusters=number_of_clusters, eigen_solver='arpack',affinity="nearest_neighbors",assign_labels="discretize")
		
		spectral = cluster.SpectralClustering(n_clusters=number_of_clusters, assign_labels="discretize")#,affinity="nearest_neighbors",n_neighbors=20)#,affinity="nearest_neighbors",assign_labels="discretize")
		# catch warnings related to kneighbors_graph
		with warnings.catch_warnings():
			warnings.filterwarnings("ignore",
				message="the number of connected components of the " +
				"connectivity matrix is [0-9]{1,2}" +
				" > 1. Completing it to avoid stopping the tree early.",
				category=UserWarning)
			warnings.filterwarnings(
				"ignore",
				message="Graph is not fully connected, spectral embedding" +
				" may not work as expected.",
				category=UserWarning)
			spectral.fit(X_new[:, 0:selected_pca])
			
		labels = spectral.labels_.astype(np.int)	
		#fpc.append(fpc_this_one)
		
		print ("input variables {}".format(columns))
		print("labels  {}; len(labels)  {}; len(np.unique(labels))  {}".format(np.unique(labels),len(labels),len(np.unique(labels))))
		if(len(np.unique(labels))!=number_of_clusters):
			print "!!!!len(np.unique(labels))!=number_of_clusters!!!"
			sil_diff = -999
		
		# create model object for skfuzzy.cmeans results
		clf = NearestCentroid(metric='euclidean')
		clf.fit(X_new[:, 0:selected_pca], labels)
		#print(clf.centroids_)
		
		s_model = spectr_model(labels,clf.centroids_)
		
		
		#print ("spectr_model.labels_.shape {}".format(s_model.labels_.shape))
		
		model_list.append(s_model)
		if (sil_diff <> -999):
			WG.append(compute_WG(s_model.labels_, X_new[:, 0:selected_pca]))
			chs.append(calinski_harabaz_score(X_new[:, 0:selected_pca], s_model.labels_))
			sil.append(silhouette_score(X_new[:, 0:selected_pca], s_model.labels_))
			dbs.append(davies_bouldin_score(X_new[:, 0:selected_pca], s_model.labels_))
		
			
		
		#print "dbs score is {}".format(dbs)
		total_cond = 0
		if((number_of_clusters > 2) & (sil_diff <> -999)):
			WG_diff = get_WG_diff(WG)
			sil_diff = sil[-1] - sil[-2]
			if (sil_diff >= 0):total_cond +=1
			chs_diff = chs[-1] - chs[-2]
			if (chs_diff >= 0):total_cond +=1
			dbs_diff = dbs[-1] - dbs[-2]
			if (dbs_diff >= 0):total_cond +=1
		
		number_of_clusters += 1
	
	pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
	#pd_X_new = pd.DataFrame(data = X_new_ica, columns = new_field_list[0:selected_pca])
	#print "len(model_list)"
	#print len(model_list)
	if((len(model_list)>=2) & (sil_diff <> -999)):
		new_labels = np.copy(model_list[-2].labels_)
		new_centroids = np.copy(model_list[-2].cluster_centers_)
	else:
		new_labels = np.copy(model_list[-1].labels_)
		new_centroids = np.copy(model_list[-1].cluster_centers_)
	#new_labels[np.max(model_list[-2].X,axis=1)<threshold] = -1
	
	pd_X_new['Clusters'] = pd.Series(new_labels)
	
	palette = sns.color_palette('deep',n_colors=(len(np.unique(new_labels))))#-1))
	print ("palette is {} np.unique(new_labels) {}".format(len(palette),np.unique(new_labels)))
	colors = [palette[z-1] if z >= 0 else (0.8, 0.8, 0.8) for z in np.unique(new_labels)]
	
	#print ("pd_X_new['Clusters']  are {}".format(np.unique(pd_X_new['Clusters'])))
	#print ("np.unique(new_lables) {}".format(np.unique(new_labels)))
	
	#childrens_list = get_active_clusteres_list(root)
	
	#figurename = output_dir + "{}_clustered&prob_PCA_scatters_pd_X_new_{}points.png".format(len(data[:,0]))
	figurename = output_dir + "{}_clustered&prob_ICA_scatters_pd_X_new_{}points_spectr.png".format(cluster_name,len(data[:,0]))
	
	#print ("plotting clusters with -1 cluster in spectr")
	#(colors,labels_list) = qvpp.clustered_scatter(pd_X_new,figurename,labels_list=[],plot_kwds=plot_kwds)
	
	labels_list = []
	for i in np.unique(new_labels):
		if(i==-1):
			labels_list.append("---")
			pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': labels_list[-1]}, regex=False)
		elif((len(model_list)<2) | (sil_diff == -999)):
			labels_list.append('cl{}'.format(i))
			pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters':labels_list[-1]}, regex=False)
		else:
			labels_list.append('cl{}'.format(i+1))
			pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters':labels_list[-1]}, regex=False)
	
	#if((len(model_list)<2) | (sil_diff == -999)):model_list[-1].add_list(colors,labels_list)
	#else:
	model_list[-2].add_list(colors,labels_list)
	
	
	#print ("number_of_clusters now is {} so, the returned value is {};model_list {} and we return {} with cluster_centers_ {}".format(number_of_clusters,number_of_clusters-2,model_list,model_list[-2],model_list[-2].cluster_centers_))
	
	pd_X_new = pd.DataFrame(data = X_new[:, 0:selected_pca], columns = new_field_list[0:selected_pca])
	pd_X_new['Clusters'] = pd.Series(model_list[-2].labels_)
	
	figurename = output_dir + "{}_clustered_PCA_scatters_pd_X_new_{}points_spectr.png".format(cluster_name,len(data[:,0]))
	(colors,labels_list) = qvpp.clustered_scatter(pd_X_new,figurename,colors=colors,labels_list=labels_list,plot_kwds=plot_kwds,centroids=new_centroids)
	
	
	X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
	dataframe_X = pd.DataFrame(data=data, columns=columns)
	dataframe_X['Clusters'] = pd.Series(model_list[-2].labels_)#new_labels)
	
	figurename = output_dir + "{}_clustered_original_scatters_pd_X_{}points_spectr.png".format(cluster_name,len(data[:,0]))
	qvpp.clustered_scatter(dataframe_X,figurename,colors=model_list[-2].colors,labels_list=labels_list,plot_kwds=plot_kwds,centroids=None)
	
	return(number_of_clusters-2,X_new,selected_pca,model_list[-2])



def get_subclusters_with_bic(data,columns=['dBZ','ZDR','RhoHV','PhiDP','temperature_2'],verbose=False):
	""" gets an optimal number of clusters in the data
		start with 2 clusters, increases the number of clusters while the optimal number is reached
	    input:		data
	    output:		number_of_clusters (int), PCA-transformed dataset, number of selected PCA, K-mean model object (sklearn)
	"""
	if verbose:print "in get_subclusters, the input data.shape is {}".format(data.shape)
	number_of_clusters = 2
	bic_diff = 0
	bic = []
	model_list = []
	# Get transformed data
	X_std = QuantileTransformer(output_distribution='normal').fit_transform(data)
	
	dataframe_X_std = pd.DataFrame(data=X_std, columns=columns)
	
	if verbose:
		print "describe X_std"
		print dataframe_X_std.describe()
	
	# Get PCA for X_std
	(selected_pca, X_new) = get_PCA(X_std)
	#map(str, range(number + 1))["X_{}".format(i) for i in range(len(columns))]
	new_field_list = ["X_{}".format(i) for i in range(len(columns))]#["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	pd_X_new = pd.DataFrame(data = X_new, columns = new_field_list)#[0:len(columns)])
	if verbose:
		print "selected_pca"
		print selected_pca
		print "describe X_new"
		print pd_X_new.describe()
		print "X_new.shape"
		print X_new.shape
		print X_new[:,0:selected_pca].shape
		
	#sns.pairplot(X_new[:,0:selected_pca], hue='species', size=2.5)
	
	while (bic_diff <= 0):
		# Get clusters
		(model) = get_k_mean(X_new[:, 0:selected_pca],k_clusters=number_of_clusters)
		#compute_gen_bic(all_medoids, np.asarray(all_labels,dtype=np.int), all_data[:,0:last_field])
		bic.append(compute_gen_bic(model.cluster_centers_, model.labels_, X_new[:, 0:selected_pca]))#silhouette_score(X_new[:, 0:selected_pca], model.labels_))
		model_list.append(model)
		if(number_of_clusters > 2):
			bic_diff = abs(bic[-1]) - abs(bic[-2])
			#print "number_of_clusters is {}".format(number_of_clusters) 
			#print "sil[-1] is {} and sil[-2] is {}. The sil_diff is {}  ".format(sil[-1],sil[-2],sil_diff)
		number_of_clusters += 1
	
	return(number_of_clusters-2,X_new,selected_pca,model_list[-2])



def get_cluster_parent_names(clust,cluster_name_list=[]):
	#print "in get_cluster_parent_names({},{})".format(clust.name,cluster_name_list)
	parent = clust.parent
	
	if(parent==None):
		#print "parent==None"
		#print "root"
		#print [clust.name]
		return [clust.name]
	else:
		cluster_name_list = get_cluster_parent_names(parent,cluster_name_list=cluster_name_list)
		#print "back to get_cluster_parent_names in the recursion"
		cluster_name_list.append(clust.name)
		#print cluster_name_list
		return cluster_name_list
	
def get_active_clusteres_list(clust,childrens_list=[],debug=False):
	"""
	recursive function
	gets cluster id
	returns list of active clusters names
	"""
	if debug:
		print "function get_active_clusteres_list was called for clust {}, which has {} subclusters".format(clust.name,len(clust.subclusters))
		print "len(childrens_list) {} and len(clust.subclusters)".format(len(childrens_list),len(clust.subclusters))
	
	if(len(clust.subclusters)==0):
		full_list_of_names = get_cluster_parent_names(clust,cluster_name_list=[])
		childrens_list.append(".".join(full_list_of_names))
		if debug:
			print "we append childrens_list by {}".format(".".join(full_list_of_names))
			print "this is an active cluster {}".format(clust.name)
			print "len(childrens_list) is {}".format(len(childrens_list))
	else:
		if debug:print "clust {} has {} subclusters".format(clust.name,len(clust.subclusters))
		for subcluster in clust.subclusters:
			if debug:
				print "working with subclaster {}".format(subcluster.name)
				#print "len(childrens_list) {}".format(len(childrens_list))
			(childrens_list) = get_active_clusteres_list(subcluster,childrens_list=childrens_list)
			#if debug:
				#print "len(childrens_list) {}".format(len(childrens_list))
	return (childrens_list)
			

	
def get_all_children_full_name(clust,childrens_list=[]):
	#print "get_all_children_full_name"
	#print clust.name
	#print childrens_list
	if(len(clust.subclusters)==0):
		full_list_of_names = get_cluster_parent_names(clust,cluster_name_list=[])
		childrens_list.append(".".join(full_list_of_names))
		
	else:
		
		for subcluster in clust.subclusters:
			#print "subcluster.name"
			#print subcluster.name
			childrens_list = get_all_children_full_name(subcluster,childrens_list=childrens_list)
			#print "childrens_list"
			#print childrens_list
	return childrens_list	
	

		

def replace_current_cluster_by_subclusters(all_childrens_list,all_labels,all_medoids,model,current_cluster=None,medoids=[],debug=False):
	""" replaces current cluster by subclusters to check the
	    input:		data
	    output:		number_of_clusters (int), PCA-transformed dataset, number of selected PCA, K-mean model object (sklearn)
	"""
	if debug:
		print "in replace_current_cluster_by_subclusters"
		print "np.unique(all_labels)"
		print np.unique(all_labels)
		print "len(all_labels)"
		print "all_medoids"
		print all_medoids
		print "np.unique(model.labels_)"
		print np.unique(model.labels_)
		print "medoids"
		print medoids
	
	new_labels = np.full_like(model.labels_,np.nan)
	if(current_cluster is not None):
		
		if(len(all_childrens_list) == 1):
			if debug:print "len(all_childrens_list) == 1"
			to_be_replaced = 0
		else:
			to_be_replaced = np.where(np.asarray(all_childrens_list)==current_cluster)[0]
		
		if debug:print "to_be_replaced {}".fomrat(to_be_replaced)
		# list of True or False coresponding to the labels
		labels_to_replace = np.where(np.asarray(all_labels) == to_be_replaced)[0]
		if debug:
			print "len(np.asarray(all_medoids)) {}".format(len(np.asarray(all_medoids)))
			print np.asarray(all_medoids)
		# replace medoid of the current cluster by the medoid of the 
		if( len(np.asarray(all_medoids)) == 1):
			if debug:print "len(all_medoids) == 1"
			all_medoids = np.asarray(medoids)
			if debug:print all_medoids
		else:
			if debug:
				print "all_medoids"
				print all_medoids
				print np.max(np.unique(model.labels_))
				print "replacing {} medoid of all_medoids {} by the {} of medoids".format(to_be_replaced,all_medoids[to_be_replaced],np.max(np.unique(model.labels_)))
			all_medoids[to_be_replaced] = np.asarray(medoids[np.max(np.unique(model.labels_))])
			all_medoids = np.append(all_medoids,medoids[:-1],axis=0)
			
		if debug:
			print "all_medoids"
			print all_medoids
		
		
		max_label = np.max(np.unique(all_labels))+1
		if debug:
			print "max_label"
			print max_label
		
		# trick to keep the original label in the label set
		# set the values of all mode labels except of the last one to original + max_label
		# set the value of the last label equal to the one to_be_replaced (keeping label in the set)
		
		
		new_labels[model.labels_!=(np.max(np.unique(model.labels_)))] = model.labels_[model.labels_!=(np.max(np.unique(model.labels_)))] + max_label
		new_labels[model.labels_==(np.max(np.unique(model.labels_)))] = to_be_replaced
		# same showld be done with medoids
		if debug:
			print "np.max(np.unique(model.labels_))"
			print np.max(np.unique(model.labels_))
		
		
		# print model.labels_[-1]
		
		# print "len(model.labels_) is {}; len of new_labels[model.labels_!=model.labels_[-1]] is {} ".format(len(model.labels_),len(new_labels[model.labels_!=model.labels_[-1]]))
		# print "len(new_labels[model.labels_==model.labels_[-1]]) is {}".format(len(new_labels[model.labels_==model.labels_[-1]]))
		if debug:
			print "unique new_labels"
			print np.unique(new_labels)
		
		for k in np.unique(model.labels_):
			
			all_labels[labels_to_replace[model.labels_ == k]] = new_labels[model.labels_ == k]
			
		
		return(all_labels,all_medoids)
	else:
		print "the current cluster was not provided"
		return(all_labels,all_medoids)

def compute_gen_bic(medoids, labels, X,verbose=False):
    """
    Computes the BIC metric for a given clusters

    Parameters:
    -----------------------------------------
    medoids:  List of clustering centers (medoids of the clusters)
    labels:   Lables of the clusters
    X     :  multidimension np array of data points

    Returns:
    -----------------------------------------
    BIC value
    """
    if verbose:print  "in compute_gen_bic"
    # assign centers to the medoids of the clusters
    centers = medoids
   
    # check if the number of medoids corresponds to the number of cluster lables
    if(len(np.unique(labels)) != len(medoids)):
		if verbose:print "The number of medoids {} is not the same as the number of clusters {}".format(len(medoids),len(np.unique(labels)))
		return (None)
    else:
		#number of clusters
		m = len(np.unique(labels))
		if verbose:
			print "number of clusters m is {}".format(m)
			print np.unique(labels)
		
    # size of the clusters
    n = np.bincount(labels)
    
    #size of data set
    N, d = X.shape
    
    
    if(N != len(labels)):
		if verbose:
			print "The numebr of lines in the lables array is not the same as the number of lines in the data array!"
			print "N is {} and d is {}".format(N,d)
			print "len(labels) is {}".format(len(labels))
		return (None) 

    #compute variance for all clusters beforehand
    #cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[i][:]], 'euclidean')**2) for i in range(m)])

    cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[i][:]], 'sqeuclidean')) for i in range(m)])

    const_term = 0.5 * m * np.log(N) * (d+1)

    BIC = np.sum([n[i] * np.log(n[i]) -
               n[i] * np.log(N) -
             ((n[i] * d) / 2) * np.log(2*np.pi*cl_var) -
             ((n[i] - 1) * d/ 2) for i in range(m)]) - const_term
    if verbose:print "BIC to be returned from compute_gen_bic is {}".format(BIC)

    return(BIC)


def compute_WG(labels, X,verbose=False):
    """
    Computes the Wemmert-Gangarski index for a given clusters

    Parameters:
    -----------------------------------------
    medoids:  List of clustering centers (medoids of the clusters)
    labels:   Lables of the clusters
    X     :  multidimension np array of data points

    Returns:
    -----------------------------------------
    WG value
    """
    #if verbose:
    print  "in compute_WG"
    # assign centers to the medoids of the clusters
    #centers = medoids
    # create model object for skfuzzy.cmeans results
    #centers = get_mean(X, labels)#get_median(X, labels)
    #print centers
    clf = NearestCentroid(metric='euclidean')
    clf.fit(X, labels)
    centers = clf.centroids_
   
    # check if the number of medoids corresponds to the number of cluster lables
    if(len(np.unique(labels)) != len(centers)):
		if verbose:print "The number of medoids {} is not the same as the number of clusters {}".format(len(centers),len(np.unique(labels)))
		return (None)
    else:
		#number of clusters
		m = len(np.unique(labels))
		if verbose:
			print "number of clusters m is {}".format(m)
			print np.unique(labels)
		
    # size of the clusters
    n = np.bincount(labels)
    
    #size of data set
    N, d = X.shape
    
    
    if(N != len(labels)):
		if verbose:
			print "The numebr of lines in the lables array is not the same as the number of lines in the data array!"
			print "N is {} and d is {}".format(N,d)
			print "len(labels) is {}".format(len(labels))
		return (None) 

    #compute variance for all clusters beforehand
    cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[i][:]], 'sqeuclidean')) for i in range(m)])
	#sum(distance.cdist(X[np.where(labels == i)], [centers[i][:]], 'euclidean')**2)
    J=[]
    All_distances = distance.cdist(X, centers, 'sqeuclidean')
    #print "All_distances {}".format(All_distances)
    
    for i in range(m):
		
		#print X.shape
		#print len(X[np.where(labels == i)])
		#if(len(X[np.where(labels == i)]) > 0):
		this_cluster_distances = distance.cdist(X[np.where(labels == i)], centers, 'sqeuclidean')
		#print "this_cluster_distances {} and shape {}".format(this_cluster_distances,np.array(this_cluster_distances).shape)
		dist_to_other_centroids = np.delete(this_cluster_distances,[i],axis=1)
		#print "dist_to_other_centroids {} and shape {}".format(dist_to_other_centroids,np.array(dist_to_other_centroids).shape)
		#print "min_dist_to_other_centroids ="
		min_dist_to_other_centroids = np.amin(dist_to_other_centroids,axis=1)
		#print min_dist_to_other_centroids
		#print "min_dist_to_other_centroids.shape"
		#print min_dist_to_other_centroids.shape
		#print np.amin(this_cluster_distances,axis=0)
		#print np.array([distance.cdist(X[np.where(labels == i)],np.array(np.delete(centers, [i], axis=0)[:]), 'sqeuclidean')][:]).shape
		#print "(distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean') for k in [np.delete(centers, [i], axis=0)[:]]"
		#print distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean') for k in [np.delete(centers, [i], axis=0)[:]]
		#print "[(distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean'),axis=0) for k in [np.delete(centers, [i], axis=0)[:]]].shape"
		
		#print np.array(distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean') for k in [np.delete(centers, [i], axis=0)[:]]).shape
		#print "all centroids"
		#min_dist_to_other_centroids = [np.amin(distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean'),axis=0) for k in centers]
		#print len(min_dist_to_other_centroids)
		#print "all except of this one"
		#min_dist_to_other_centroids = [np.amin(distance.cdist(X[np.where(labels == i)], np.array(k), 'sqeuclidean'),axis=0) for k in [np.delete(centers, [i], axis=0)[:]]]
		#print "len(min_dist_to_other_centroids)"
		#print len(min_dist_to_other_centroids)
		#print min_dist_to_other_centroids
		#print np.array(distance.cdist(X[np.where(labels == i)], [centers[i][:]], 'sqeuclidean')).shape
		#print this_cluster_distances[:,i].shape
		R = np.divide(this_cluster_distances[:,i],min_dist_to_other_centroids)#,where=np.where(min_dist_to_other_centroids<>0))
		#print R.shape
		R_sum = sum(R)
		#print "R_sum"
		#print R_sum
		#print n[i]
		#print 1./n[i]
		
		J.append(np.max([0,(1 - (1./n[i]) * R_sum)]))
    WG = (1./N) * sum(np.multiply(np.array(J),n))
    print "WG is {}".format(WG)
   
    if verbose:print "WG to be returned from compute_WG is {}".format(WG)

    return(WG)

def check_number_of_labels(n_labels, n_samples):
    """Check that number of labels are valid.
    Parameters
    ----------
    n_labels : int
        Number of labels
    n_samples : int
        Number of samples
    """
    if not 1 < n_labels < n_samples:
        raise ValueError("Number of labels is %d. Valid values are 2 "
                         "to n_samples - 1 (inclusive)" % n_labels)

def _silhouette_reduce(D_chunk, start, labels, label_freqs):
    """Accumulate silhouette statistics for vertical chunk of X
    Parameters
    ----------
    D_chunk : shape (n_chunk_samples, n_samples)
        precomputed distances for a chunk
    start : int
        first index in chunk
    labels : array, shape (n_samples,)
        corresponding cluster labels, encoded as {0, ..., n_clusters-1}
    label_freqs : array
        distribution of cluster labels in ``labels``
    """
    # accumulate distances from each sample to each cluster
    clust_dists = np.zeros((len(D_chunk), len(label_freqs)),
                           dtype=D_chunk.dtype)
    for i in range(len(D_chunk)):
        clust_dists[i] += np.bincount(labels, weights=D_chunk[i],
                                      minlength=len(label_freqs))

    # intra_index selects intra-cluster distances within clust_dists
    intra_index = (np.arange(len(D_chunk)), labels[start:start + len(D_chunk)])
    # intra_clust_dists are averaged over cluster size outside this function
    intra_clust_dists = clust_dists[intra_index]
    # of the remaining distances we normalise and extract the minimum
    clust_dists[intra_index] = np.inf
    clust_dists /= label_freqs
    inter_clust_dists = clust_dists.min(axis=1)
    return intra_clust_dists, inter_clust_dists

def silhouette_samples(X, labels, metric='euclidean', **kwds):
	"""Compute the Silhouette Coefficient for each sample.
	The Silhouette Coefficient is a measure of how well samples are clustered
	with samples that are similar to themselves. Clustering models with a high
	Silhouette Coefficient are said to be dense, where samples in the same
	cluster are similar to each other, and well separated, where samples in
	different clusters are not very similar to each other.
	The Silhouette Coefficient is calculated using the mean intra-cluster
	distance (``a``) and the mean nearest-cluster distance (``b``) for each
	sample.  The Silhouette Coefficient for a sample is ``(b - a) / max(a,
	b)``.
	Note that Silhouette Coefficient is only defined if number of labels
	is 2 <= n_labels <= n_samples - 1.
	This function returns the Silhouette Coefficient for each sample.
	The best value is 1 and the worst value is -1. Values near 0 indicate
	overlapping clusters.
	Read more in the :ref:`User Guide <silhouette_coefficient>`.
	Parameters
	----------
	X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
			 [n_samples_a, n_features] otherwise
		Array of pairwise distances between samples, or a feature array.
	labels : array, shape = [n_samples]
			 label values for each sample
	metric : string, or callable
		The metric to use when calculating distance between instances in a
		feature array. If metric is a string, it must be one of the options
		allowed by :func:`sklearn.metrics.pairwise.pairwise_distances`. If X is
		the distance array itself, use "precomputed" as the metric.
	`**kwds` : optional keyword parameters
		Any further parameters are passed directly to the distance function.
		If using a ``scipy.spatial.distance`` metric, the parameters are still
		metric dependent. See the scipy docs for usage examples.
	Returns
	-------
	silhouette : array, shape = [n_samples]
		Silhouette Coefficient for each samples.
	References
	----------
	.. [1] `Peter J. Rousseeuw (1987). "Silhouettes: a Graphical Aid to the
	   Interpretation and Validation of Cluster Analysis". Computational
	   and Applied Mathematics 20: 53-65.
	   <https://www.sciencedirect.com/science/article/pii/0377042787901257>`_
	.. [2] `Wikipedia entry on the Silhouette Coefficient
	   <https://en.wikipedia.org/wiki/Silhouette_(clustering)>`_
	"""
	X, labels = utl.check_X_y(X, labels, accept_sparse=['csc', 'csr'])
	
	le = prp.LabelEncoder()
	labels = le.fit_transform(labels)
	n_samples = len(labels)
	label_freqs = np.bincount(labels)
	check_number_of_labels(len(le.classes_), n_samples)
	
	kwds['metric'] = metric
	
	reduce_func = functools.partial(_silhouette_reduce,
									labels=labels, label_freqs=label_freqs)
	results = zip(*pairwise_distances_chunked(X, reduce_func=reduce_func,
											  **kwds))
	intra_clust_dists, inter_clust_dists = results
	intra_clust_dists = np.concatenate(intra_clust_dists)
	inter_clust_dists = np.concatenate(inter_clust_dists)
	
	denom = (label_freqs - 1).take(labels, mode='clip')
	with np.errstate(divide="ignore", invalid="ignore"):
		intra_clust_dists /= denom
	
	sil_samples = inter_clust_dists - intra_clust_dists
	
	
	with np.errstate(divide="ignore", invalid="ignore"):
		sil_samples /= np.maximum(intra_clust_dists, inter_clust_dists)
	# nan values are for clusters of size 1, and should be 0
	
	return np.nan_to_num(sil_samples)
		
   
def get_fuzz_SIL(X, u, metric='euclidean'):
	"""
	get fuzzy Sihouette score calculated by SIL.F function of the fclust R pachage
	input:	
		X				- array with the points (n_points x n_variables) 
		u				- membership degree matrix 
		alpha=1 		- weighting cefficient 
		distance	 	- if True X is assumed to be distance matrix (default: False)
	output:	f_sil - fuzzy Silouette index value defind as in Campello&Hruschka 2006 
	"""
	
	labels = np.argmax(u,axis=0)
	sil_j = silhouette_samples(X, labels, metric=metric)#, **kwds)
	# indexes of two maximum u values per point
	u_max = u.argsort(axis=0)[-2:][::-1]
	u_pj = u_max[0]
	u_qj = u_max[1]
	
	# equasion (16)
	# sil = np.mean(sil_j)
	f_sil = np.sum(np.multiply((u[u_pj] - u[u_qj]),sil_j)) / np.sum(u[u_pj] - u[u_qj])
	
	return f_sil
	
def get_median(X, labels):
	centroids = []
	for l in np.unique(labels):
		centroid = np.median(X[labels == l,:], axis=0)
		centroids.append(centroid)
	return centroids

def get_mean(X, labels):
	centroids = []
	for l in np.unique(labels):
		centroid = np.mean(X[labels == l,:], axis=0)
		centroids.append(centroid)
	return centroids

def get_WG_diff(gen_bic):
	"""
	get WG difference for the last to values in the list.
	
	input:	list of calculated WG values
	output:	True if the difference is positive (the last split doesn't help => the split before that was optimal)
			False if the difference is negative (the last split makes the total clustering better)
	"""
	# to add calculation of BIC and return True is the number of clusters is optimal
	if(gen_bic[-1] > gen_bic[-2]):return True
	else: return False
	
def get_BIC_diff(gen_bic):
	"""
	get BIC difference for the last to values in the list
	input:	list of calculated BIC values
	output:	True if the difference is positive (the last split doesn't help => the split before that was optimal)
			False if the difference is negative (the last split makes the total clustering better)
	"""
	# to add calculation of BIC and return True is the number of clusters is optimal
	if(gen_bic[-1] < gen_bic[-2]):return True
	else: return False

def plot_BIC(gen_bic,number_of_clusteres,file_name="BIC_curve.png", output_dir="./data/",verbose=False):
	fig = plt.figure(1, figsize=(6, 6))
	ax = fig.add_subplot(111)

	#plt.xticks([0,1,2,3,4], a[:,0])
	if verbose:print "number_of_clusteres = {}".format(number_of_clusteres)
	plt.plot(number_of_clusteres,np.asarray(gen_bic).astype(float))
	
	plt.title('BIC scores through the loops')
	savefile = "{}{}".format(output_dir,file_name)
	plt.savefig(savefile , format='png',dpi=200)
	plt.close()
	if verbose:print "figure %s is ok" %savefile
	
	return savefile
	
def plot_WG(gen_WG,number_of_clusteres,file_name="WG_curve.png", output_dir="./data/",verbose=False):
	fig = plt.figure(1, figsize=(6, 6))
	ax = fig.add_subplot(111)

	#plt.xticks([0,1,2,3,4], a[:,0])
	if verbose:print "number_of_clusteres = {}".format(number_of_clusteres)
	plt.plot(number_of_clusteres,np.asarray(gen_WG).astype(float))
	
	plt.title('WG scores through the loops')
	savefile = "{}{}".format(output_dir,file_name)
	plt.savefig(savefile , format='png',dpi=200)
	plt.close()
	if verbose:print "figure %s is ok" %savefile
	
	return savefile

def get_PCA(X_std,evr_threshold=0.85,debug=False):
	"""
	Based on the transformed data get the number of PCA comoponenets and the representation of X_std in the selected components
	"""
	
	pca =  PCA()
	if(debug): 
		print ("in  get_PCA")
		print ("input data shape is {}".format(X_std.shape)) 
		print ("explained_variance_ratio threshold is {}".format(evr_threshold))
	X_new = pca.fit_transform(X_std)
	#debug = True
	if(debug):
		print "pca.explained_variance_ratio_"
		print pca.explained_variance_ratio_
		print np.sum(pca.explained_variance_ratio_)
		print "pca.components_"
		print pca.components_
	if(np.sum(pca.explained_variance_ratio_[0:(np.max(np.where(pca.explained_variance_ratio_ >= 0.10))+1)]) >= evr_threshold):
		selected_pca = np.max(np.where(pca.explained_variance_ratio_ >= 0.10)) + 1
	else:
		selected_pca = np.max(np.where(pca.explained_variance_ratio_ >= 0.10)) + 2
	if(debug):
		print "selected_pca"
		print selected_pca
		
	return (selected_pca, X_new)

def get_ICA(X_std,evr_threshold=0.85,n_components=None,debug=False):
	"""
	Based on the transformed data get the number of PCA comoponenets and the representation of X_std in the selected components
	"""
	
	ica = FastICA(n_components=n_components)
	X_new = ica.fit_transform(X_std)  # Estimate the sources
	#X_new /= X_new.std(axis=0)


	if(debug): 
		print ("in  get_ICA")
		print ("input data shape is {}".format(X_std.shape)) 
		print ("explained_variance_ratio threshold is {}".format(evr_threshold))
	#X_new = pca.fit_transform(X_std)
	#debug = True
	if(debug):
		print ("ica.mixing_.T")
		print (ica.mixing_.T)
		print ("ica.components_")
		print (ica.components_)
		print np.sum(ica.mixing_.T)
		
	# if(np.sum(pca.explained_variance_ratio_[0:(np.max(np.where(pca.explained_variance_ratio_ >= 0.10))+1)]) >= evr_threshold):
		# selected_pca = np.max(np.where(pca.explained_variance_ratio_ >= 0.10)) + 1
	# else:
		# selected_pca = np.max(np.where(pca.explained_variance_ratio_ >= 0.10)) + 2
	# if(debug):
		# print "selected_pca"
		# print selected_pca
		
	return X_new#(selected_pca, X_new)


def fuzzy_get_u(x, centroids, m):
    distances = distance.cdist(x, centroids)**2#pairwise_squared_distances(x, centroids)
    nonzero_distances = np.fmax(distances, np.finfo(np.float64).eps)
    inv_distances = np.reciprocal(nonzero_distances)**(1/(m - 1))
    return inv_distances.T/np.sum(inv_distances, axis=1)

def get_k_mean(X_new,k_clusters=2,debug=False):
	# Clusters on PCA variables
	if debug:
		print "in get_k_mean"
		print X_new.shape
		print "k_clusters = {}".format(k_clusters)
	results = []
	algorithms = {}
	algorithms['kmeans'] = cluster.KMeans(n_clusters=k_clusters, n_init=200)#, max_iter=500,tol=0.0001)
	for model in algorithms.values():
		if debug:
			print "algorithms.values()"
			print model
		#model.fit(X_new[:, 0:2]) -  for ERAD plots it was like that
		model.fit(X_new)
		results.append(list(model.labels_))
		algorithms['kmeans'].fit(X_new)#[:, 0:2])
		if debug:
			print  model.labels_
			print model.inertia_
			print "model.cluster_centers_"
			print model.cluster_centers_
			
	return (model)



# def main(args):
    # return 0

# if __name__ == '__main__':
    # import sys
    # sys.exit(main(sys.argv))
