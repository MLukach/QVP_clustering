
import copy
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import colors
import seaborn as sns

from matplotlib.colors import BoundaryNorm

import numpy as np
import pandas as pd

from pandas.plotting import scatter_matrix

from netCDF4 import num2date
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools

plot_kwds = {'alpha' : 0.25, 's' : 80, 'linewidths':0}
font = {'weight': 'normal',
        'size': 14,
        }
        #'family': 'serif',

def line_plot(data,y,columns,labels, title, cluster_label, outputdir, palette=None):
	
	fig, axes = plt.subplots(1,6,figsize=(15,5))#, sharey=True,)
	
	
	#print labels
	#print "axes"
	#print axes
	ax_ind = 0
	#print "data.describe()"
	#print data.describe()
	for i in columns:
		#print "axes[0, {}].plot(data['{}'], data['{}'], colour=palette)".format(ax_ind,y,i)
		#print data[i].as_matrix()
		#print data[y].as_matrix()
		#data.plot(x=i,y=y,kind='line',ax=axes[0, ax_ind])
		
		axes[ax_ind].plot(data[i].values(), data[y].values(), color=palette)
		axes[ax_ind].set_xlabel(labels[ax_ind], fontsize=14)
		if(ax_ind == 0):axes[0].set_ylabel(y, fontsize=14)
		else:axes[ax_ind].set_ylabel('')
		#sns.lineplot(y=y, x=i,hue='region', style='event',data=df,ax=axes[0, ax_ind])
		#print "sns.lineplot(y='{}', x='{}',hue='region', style='event',data=df)".format(y,i)
		ax_ind += 1
		
	
	# output file name
	fig.subplots_adjust(top=0.88)
	plt.setp([a.get_yticklabels() for a in axes[1:]], visible=False)
	plot_file_name="{}{}_lineplot.jpg".format(outputdir, cluster_label)
	print plot_file_name
	# save as jpeg
	#plt.ylabel(y)
	#axes[0].set_ylabel(y)
	
	#plt.tight_layout()
	fig.suptitle("{} {}".format(cluster_label,title), fontsize=16)
	plt.savefig(plot_file_name, format='png', dpi=200)
	#(pad=0.4, w_pad=0.5, h_pad=1.0)
	plt.close()
	return plot_file_name

def line_and_error_plot(data,y,columns,labels,limits, title, cluster_label, outputdir, palette=None):
	
	fig, axes = plt.subplots(1,6,figsize=(15,5))#, sharey=True,)
	
	
	#print labels
	#print "axes"
	#print axes
	ax_ind = 0
	#print "data.describe()"
	#print data.describe()
	for i in columns:
		#print "y is {} i is {} from columns {}".format(y,i,columns)
		#print "axes[0, {}].plot(data['{}'], data['{}'], colour=palette)".format(ax_ind,y,i)
		#print data
		#print data[i].as_matrix()
		#data.plot(x=i,y=y,kind='line',ax=axes[0, ax_ind])
		
		if (i is 'counts'):
			axes[ax_ind].plot(data[i].as_matrix(), data[y].as_matrix(), color=palette)
		else:
			axes[ax_ind].errorbar(data[i].as_matrix(), data[y].as_matrix(), xerr=data['{}_std'.format(i)].as_matrix(), color=palette,ecolor='#bebebe')
			axes[ax_ind].set_xlim(limits[ax_ind-1])
		axes[ax_ind].set_xlabel(labels[ax_ind], fontsize=14)
		axes[ax_ind].set_ylim([0,1500])
		if(ax_ind == 0):axes[0].set_ylabel(y, fontsize=14)
		else:axes[ax_ind].set_ylabel('')
		axes[ax_ind].set_axisbelow(True)
		axes[ax_ind].yaxis.grid(color='gray', linestyle='dashed')
		axes[ax_ind].xaxis.grid(color='gray', linestyle='dashed')

		#sns.lineplot(y=y, x=i,hue='region', style='event',data=df,ax=axes[0, ax_ind])
		#print "sns.lineplot(y='{}', x='{}',hue='region', style='event',data=df)".format(y,i)
		ax_ind += 1
		
	
	# output file name
	fig.subplots_adjust(top=0.88)
	plt.setp([a.get_yticklabels() for a in axes[1:]], visible=False)
	plot_file_name="{}{}_lineplot.jpg".format(outputdir, cluster_label)
	print plot_file_name
	# save as jpeg
	#plt.ylabel(y)
	#axes[0].set_ylabel(y)
	
	#plt.tight_layout()
	fig.suptitle("{} {}".format(cluster_label,title), fontsize=16)
	plt.savefig(plot_file_name, format='png', dpi=200)
	#(pad=0.4, w_pad=0.5, h_pad=1.0)
	plt.close()
	return plot_file_name

def line_and_error_plot_per_case(data,y,columns,labels,limits, title, cluster_label, case_label, outputdir, palette=None):
	
	fig, axes = plt.subplots(1,6,figsize=(15,5))#, sharey=True,)
	
	
	#print labels
	#print "axes"
	#print axes
	ax_ind = 0
	#print "data.describe()"
	#print data.describe()
	for i in columns:
		#print "axes[0, {}].plot(data['{}'], data['{}'], colour=palette)".format(ax_ind,y,i)
		#print data[i].as_matrix()
		#print data[y].as_matrix()
		#data.plot(x=i,y=y,kind='line',ax=axes[0, ax_ind])
		
		if (i is 'counts'):
			axes[ax_ind].plot(data[i].as_matrix(), data[y].as_matrix(), color=palette)
		else:
			axes[ax_ind].errorbar(data[i].as_matrix(), data[y].as_matrix(), xerr=data['{}_std'.format(i)].as_matrix(), color=palette,ecolor='#bebebe')
			axes[ax_ind].set_xlim(limits[ax_ind-1])
		axes[ax_ind].set_xlabel(labels[ax_ind], fontsize=14)
		axes[ax_ind].set_ylim([0,8000])
		if(ax_ind == 0):axes[0].set_ylabel(y, fontsize=14)
		else:axes[ax_ind].set_ylabel('')
		axes[ax_ind].set_axisbelow(True)
		axes[ax_ind].yaxis.grid(color='gray', linestyle='dashed')
		axes[ax_ind].xaxis.grid(color='gray', linestyle='dashed')

		#sns.lineplot(y=y, x=i,hue='region', style='event',data=df,ax=axes[0, ax_ind])
		#print "sns.lineplot(y='{}', x='{}',hue='region', style='event',data=df)".format(y,i)
		ax_ind += 1
		
	
	# output file name
	fig.subplots_adjust(top=0.88)
	plt.setp([a.get_yticklabels() for a in axes[1:]], visible=False)
	plot_file_name="{}{}_lineplot_{}.jpg".format(outputdir, cluster_label,case_label)
	print plot_file_name
	# save as jpeg
	#plt.ylabel(y)
	#axes[0].set_ylabel(y)
	
	#plt.tight_layout()
	fig.suptitle("{} {}".format(cluster_label,title), fontsize=16)
	plt.savefig(plot_file_name, format='png', dpi=200)
	#(pad=0.4, w_pad=0.5, h_pad=1.0)
	plt.close()
	return plot_file_name


def box_plot(df,variable, label, title, cluster_labels, outputdir, palette=None):
	# Reset default params
	#sns.set()

	# Set context to `"paper"`
	#sns.set_context("paper")
	#file1 = open("{}.txt".format(variable),"a")
	
	fig, ax = plt.subplots(figsize=(8,4))
	
	print "sns.boxplot(y={}, x='cluster',".format(variable)
	# if (variable == "temperature_2"):
		# if(palette is None):
			# bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
	                 # data=(df-272.15), 
	                 # width=0.5,
	                 # palette=sns.color_palette("hls", 8))#"colorblind")
		# else:
			# bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
	                 # data=(df-272.15), 
	                 # width=0.5,
	                 # palette=palette)#"colorblind")
			
	# else:
	
	if(palette is None):
		bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
                 data=df, 
                 width=0.5,
                 palette=sns.color_palette("hls", 8))#"colorblind")
	else:
		bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
                 data=df, 
                 width=0.5,
                 palette=palette)#"colorblind")
	#bplot.axes.set_title(title,fontsize=14)
	
	
	bplot.set_xlabel("Clusters", fontsize=12, fontname="Times New Roman")
 
	bplot.set_ylabel(label,fontsize=12, fontname="Times New Roman")
 
	bplot.tick_params(labelsize=12)
	bplot.yaxis.grid(True)
	
	bplot.set_xticklabels(cluster_labels,fontsize=12, fontname="Times New Roman")#rotation=0,
	#bplot.despine(left=True
	#file1.writelines(L) for L = [str1, str2, str3]

	# output file name
	plot_file_name="{}{}_boxplot_{}.jpg".format(outputdir,variable,len(cluster_labels))
	print plot_file_name
	if(variable == 'case'):
		for i in df.cluster.unique():
			#print i
			#print df[df.cluster==i].describe()
			df[df.cluster==i].describe().to_excel(r'{}_{}_cluster_stats_of{}clusters.xls'.format(plot_file_name[0:-17],i,len(df.cluster.unique())), header=True, index=True)#, mode='a')
			print r'{}_{}_cluster_stats_of{}clusters.xls'.format(plot_file_name[0:-17],i,len(df.cluster.unique()))
	# save as jpeg
	#plt.tight_layout()
	plt.title(title, fontsize=20, fontname="Times New Roman")
	bplot.figure.savefig(plot_file_name, format='png', dpi=300)
	#(pad=0.4, w_pad=0.5, h_pad=1.0)
	del bplot
	#file1.close()
	return plot_file_name

def cat_plot(df,variable, label, title, cluster_labels, outputdir, palette=None):
	# Reset default params
	#sns.set()

	# Set context to `"paper"`
	#sns.set_context("paper")
	fig, ax = plt.subplots(figsize=(10,5))
	
	print "sns.catplot(y={}, x='cluster',".format(variable)
	# if (variable == "temperature_2"):
		# if(palette is None):
			# bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
	                 # data=(df-272.15), 
	                 # width=0.5,
	                 # palette=sns.color_palette("hls", 8))#"colorblind")
		# else:
			# bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
	                 # data=(df-272.15), 
	                 # width=0.5,
	                 # palette=palette)#"colorblind")
			
	# else:
	
	#sns.catplot(x="who", y="survived", col="class",data=titanic, saturation=.5,kind="bar", ci=None, aspect=.6)
	
	if(palette is None):
		bplot=sns.catplot(y="{}".format(variable), x='cluster', 
                 data=df, 
                 kind="bar", ci=None, aspect=.6,
                 palette=sns.color_palette("hls", 13))#"colorblind")
	else:
		bplot=sns.boxplot(y="{}".format(variable), x='cluster', 
                 data=df, 
                 width=0.5,
                 kind="bar", ci=None, aspect=.6,
                 palette=palette)#"colorblind")
	#bplot.axes.set_title(title,fontsize=14)
	bplot.set_xlabel("Clusters", fontsize=18)
 
	bplot.set_ylabel(label,fontsize=18)
 
	bplot.tick_params(labelsize=12)
	bplot.yaxis.grid(True)
	
	bplot.set_xticklabels(cluster_labels, rotation=0)
	#bplot.despine(left=True)

	# output file name
	plot_file_name="{}{}_catplot.jpg".format(outputdir,variable)
	print plot_file_name
	# save as jpeg
	#plt.tight_layout()
	plt.title(title, fontsize=20)
	bplot.figure.savefig(plot_file_name, format='png', dpi=300)
	#(pad=0.4, w_pad=0.5, h_pad=1.0)
	del bplot
	return plot_file_name



def leveled_cmap(levels,cmap):
    N = len(levels)-1
    cmaplist = [cmap(i) for i in range(0,cmap.N,cmap.N/N)]
    output_cmap = colors.ListedColormap(cmaplist)
    output_bnorm = BoundaryNorm(levels, ncolors=N, clip=False)
    return(output_cmap,output_bnorm)

def leveled_discreet_cmap(levels,colors_list):
    N = len(levels)-1
    cmaplist = colors_list
    output_cmap = colors.ListedColormap(cmaplist)
    output_bnorm = BoundaryNorm(levels, ncolors=N, clip=False)
    return(output_cmap,output_bnorm)


def plot_field(qvp_file,
               field,
               Cnorm,
               Cmap,
               x_lims,
               vmin, vmax,
               count_threshold=180,variat_threshold=0.3,offset=0,scaler=1.0, ax=None,night_mask=None):
    alts = qvp_file['Height']
    field_count = field+'/Counts'
    field_mean = field+'/Means'
    field_stddev = field+'/StdDevs'
    # if field is 'temperature_2': 
		# print qvp_file[field_mean][:]
		# print (scaler*qvp_file[field_mean][:])+offset
   
    plot_field = np.where((qvp_file[field_count][:]>count_threshold),(scaler*qvp_file[field_mean][:])+offset,np.nan)
    if(night_mask is not None):plot_field = np.where((night_mask==True),plot_field,np.nan);
    
    if( ax is None):im = plt.imshow(plot_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax,cmap=Cmap, norm=Cnorm)
    else:im = ax.imshow(plot_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax,cmap=Cmap, norm=Cnorm)
    return im

def plot_cluster_field(clusters_field,
               Cmap,
               Cnorm,
               x_lims,
               vmin, vmax,alts=None,
               count_threshold=None, ax=None):
    
    #print np.nanmax(np.unique(clusters_field))
    palette = sns.color_palette('deep', int(np.nanmax(np.unique(clusters_field))) + 1)
    #print palette
    #colors = [palette[z] if z >= 0 else (0.0, 0.0, 0.0) for z in clusters_field[~np.isnan(clusters_field)].astype(int)]
    #print colors
    #cmap = colors.ListedColormap(['b','y','g','r'])
    #print Cmap
    bounds=[0,1,2,3]
    norm = colors.BoundaryNorm(bounds, Cmap.N)
    if( ax is None):im = plt.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    else:im = ax.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    return im


# def plot_cluster_field_with_faam(clusters_field,
               # Cmap,
               # Cnorm,
               # x_lims,
               # vmin, vmax,alts=None,
               # count_threshold=None, ax=None, hs=None, hvis=None):
    
    # #print np.nanmax(np.unique(clusters_field))
    # palette = sns.color_palette('deep', int(np.nanmax(np.unique(clusters_field))) + 1)
    # #print palette
    # #colors = [palette[z] if z >= 0 else (0.0, 0.0, 0.0) for z in clusters_field[~np.isnan(clusters_field)].astype(int)]
    # #print colors
    # #cmap = colors.ListedColormap(['b','y','g','r'])
    # #print Cmap
    # bounds=[0,1,2,3]
    # norm = colors.BoundaryNorm(bounds, Cmap.N)
    # print "ax !!!"
    # print ax
    
    # if(hs is not None):hs.plot(ax=ax,color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	
    # if( ax is None):im = plt.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    # else:im = ax.imshow(clusters_field,aspect='auto',origin='lower',extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),interpolation='none',vmin=vmin, vmax=vmax, cmap=Cmap)#, norm=Cnorm)
    
    # print "hs !!!!"
    # print hs
    
    	
    
    # return im


def plot_stddev_field(qvp_file,
               field,
               Cnorm,
               Cmap,
               x_lims,
               vmin, vmax,
               count_threshold=180,variat_threshold=0.3,offset=0,scaler=1.0, ax=None):
    alts = qvp_file['Height']
    field_count = field+'/Counts'
    field_mean = field+'/Means'
    field_stddev = field+'/StdDevs'
    data_count = np.empty(np.shape(qvp_file[field_count]))*np.nan
    data_count[:] = qvp_file[field_count][:]
    
    data_mean = np.empty(np.shape(qvp_file[field_mean]))*np.nan
    data_mean[:] = qvp_file[field_mean][:]
    data_stddev = np.empty(np.shape(qvp_file[field_stddev]))*np.nan
    data_stddev[:] = qvp_file[field_stddev][:]
    data_count[data_count[:] is np.inf] = np.nan
    data_mean[data_mean[:] is np.inf] = np.nan
    data_stddev[data_stddev[:] is np.inf] = np.nan
    data_count[data_count[:] is -np.inf] = np.nan
    data_mean[data_mean[:] is -np.inf] = np.nan
    data_stddev[data_stddev[:] is -np.inf] = np.nan
    data_mean[data_mean[:] == 0 ] = 0.00001
    data_stddev[data_stddev[:] == 0 ] = 0.00001
    # print "in plot_stddev_field"
    # print field
    # print "Coefficient of Variation"
    # print len(data_mean[data_mean == 0])
    nan_max_var = np.nanmax(np.divide(data_stddev[data_mean is not np.nan],data_mean[data_mean is not np.nan]))
    nan_min_var = np.nanmin(np.divide(data_stddev[data_mean is not np.nan],data_mean[data_mean is not np.nan]))
    
    # print np.nanmax(data_mean)
    # print np.nanmin(data_mean)
    nan_mean_var = (np.nanmax(data_mean)-np.nanmin(data_mean))
    sigma_mu = np.divide(data_stddev,nan_mean_var)
    sigma_mean_var = np.divide(data_stddev, nan_mean_var)
    # print len(sigma_mu[sigma_mu<0.05])
    # print np.nanmax(np.divide(data_stddev[:],(nan_mean_var)))
    #print len(np.where(((qvp_file[field_count][:]>count_threshold)
    #                    &(np.divide(qvp_file[field_stddev][:],qvp_file[field_mean][:])<variat_threshold)
    #                    &(qvp_file[field_stddev][:]/(np.nanmax(qvp_file[field_mean])-np.nanmin(qvp_file[field_mean]))<variat_threshold)))[0])
    #print qvp_file[field_stddev][np.where(((qvp_file[field_count][:]>count_threshold)
    #                    &(np.divide(qvp_file[field_stddev][:],qvp_file[field_mean][:])<variat_threshold)
    #                    &(qvp_file[field_stddev][:]/(np.nanmax(qvp_file[field_mean])-np.nanmin(qvp_file[field_mean]))<variat_threshold)))]
    
    ###&(np.divide(qvp_file[field_stddev][:],qvp_file[field_mean][:])<variat_threshold)&(qvp_file[field_stddev][:]/(np.nanmax(qvp_file[field_mean])-np.nanmin(qvp_file[field_mean]))<variat_threshold))
    if( ax is None):
		#& (sigma_mu<variat_threshold),
		im = plt.imshow(np.where((qvp_file[field_count][:]>count_threshold), 
                        (scaler*qvp_file[field_stddev][:])+offset,
                        np.nan),
                        aspect='auto',origin='lower',
                        extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),
                        interpolation='none',
                        vmin=vmin, vmax=vmax,
                        cmap=Cmap, norm=Cnorm)
        #& (sigma_mu<variat_threshold)
    else:im = ax.imshow(np.where((qvp_file[field_count][:]>count_threshold) ,
                        (scaler*qvp_file[field_stddev][:])+offset,
                        np.nan),
                        aspect='auto',origin='lower',
                        extent =(x_lims[0], x_lims[1],alts[0],alts[-1]),
                        interpolation='none',
                        vmin=vmin, vmax=vmax,
                        cmap=Cmap, norm=Cnorm)
    return im

def annotate_field_plot(cbar_label,Cnorm,hmax,date_format,xti,x0=0,x1=-1,hmin=0, ax=None, im=None):
	if( ax is None):
		
	    plt.colorbar(ticks=Cnorm.boundaries)#label=cbar_label,
	    #cbar = plt.colorbar()
	    #tick_locs = (np.arange(len(cbar_label)+1) + 0.5)*(len(cbar_label)-1)/len(cbar_label)
	    #tick_locs = (np.arange(n_clusters+1) + 0.5)#*(n_clusters)/n_clusters
	    #cbar.set_ticks(ticks_locs)
	    #cbar.set_ticklabels(cbar_label)
	    plt.ylim(hmin,hmax)
	    plt.ylabel('Height (m)', fontdict=font)
	    plt.xlabel('Time (UTC)', fontdict=font)
	    #print date_format
	    plt.gca().xaxis.set_major_formatter(date_format)
	    #print xti
	    plt.gca().xaxis.set_ticks(xti)
	    plt.grid()
	    plt.xlim(xti[x0],xti[x1])
	else:
	    divider = make_axes_locatable(ax)
	    cax = divider.append_axes('right', size='5%', pad=0.05)
	    plt.colorbar(im,cax=cax,ticks=Cnorm.boundaries)#label=cbar_label,
	    #cbar = plt.colorbar()
	    #tick_locs = (np.arange(len(cbar_label)+1) + 0.5)*(len(cbar_label)-1)/len(cbar_label)
	    ##tick_locs = (np.arange(n_clusters+1) + 0.5)#*(n_clusters)/n_clusters
	    #cbar.set_ticks(ticks_locs)
	    #cbar.set_ticklabels(cbar_label)
	    
	    ax.set_ylim(hmin,hmax)
	    ax.set_ylabel('Height (m)', fontdict=font)
	    ax.set_xlabel('Time (UTC)', fontdict=font)
	    #print date_format
	    ax.xaxis.set_major_formatter(date_format)
	    #print xti
	    ax.xaxis.set_ticks(xti)
	    ax.set_title(cbar_label, fontsize=14,fontweight= 'bold')#,verticalalignment='top')
	    ax.grid()
	    ax.set_xlim(xti[x0],xti[x1])

def annotate_cluster_plot(cbar_label,cluster_labels,hmax,date_format,xti,x0=0,x1=-1,hmin=0, ax=None, im=None):
	#print "!!!!! In annotate_cluster_plot"
	n_clusters = len(cluster_labels)
	# print "n_clusters"
	# print n_clusters
	# print "cluster_labels"
	# print cluster_labels
	
	if( ax is None):
	    #print "ax is None"
	    cbar = plt.colorbar()

# set ticks locations (not very elegant, but it works):
# - shift by 0.5
# - scale so that the last value is at the center of the last color
	    tick_locs = (np.arange(n_clusters) + 0.5)*(n_clusters-1)/n_clusters
	    #tick_locs = (np.arange(n_clusters+1) + 0.5)#*(n_clusters)/n_clusters
	    cbar.set_ticks(ticks_locs)

# set tick labels (as before)
	    cbar.set_ticklabels(cluster_labels)#np.arange(n_clusters))
	    plt.ylim(hmin,hmax)
	    plt.ylabel('Height (m)', fontdict=font)
	    plt.gca().xaxis.set_major_formatter(date_format)
	    plt.gca().xaxis.set_ticks(xti)
	    plt.grid()
	    plt.xlim(xti[x0],xti[x1])
	else:
	    #print "ax is not None"
		
	    divider = make_axes_locatable(ax)
	    cax = divider.append_axes('right', size='5%', pad=0.05)
	    cbar = plt.colorbar(im,cax=cax)

# set ticks locations (not very elegant, but it works):
# - shift by 0.5
# - scale so that the last value is at the center of the last color
	    #tick_locs = (np.arange(n_clusters) + 0.5)*(n_clusters-1)/n_clusters
	    tick_locs = (np.arange(n_clusters+1) + 0.5)
	    #print "tick_locs"
	    #print tick_locs
	    cbar.set_ticks(tick_locs)
	    #print "cluster_labels"
	    #print cluster_labels
	    cbar.set_ticklabels(cluster_labels)#np.arange(n_clusters))
	    ax.set_ylim(hmin,hmax)
	    ax.set_ylabel('Height (m)', fontdict=font)
	    ax.xaxis.set_major_formatter(date_format)
	    #print "xti"
	    #print xti
	    ax.xaxis.set_ticks(xti)
	    ax.grid()
	    ax.set_xlim(xti[x0],xti[x1])



def new_fields_scatters_plot(output_dir,X_new,selected_pca,new_field_list,plot_kwds=plot_kwds):
	figurename = output_dir + "new_fields_scatters.png"
	
	new_field_list = ["X_0","X_1","X_2","X_3","X_4","X_5","X_6"]
	
	pd_X_new = pd.DataFrame(data = X_new, columns = new_field_list)
	scatter_matrix(pd_X_new,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
	
	savefile2 = output_dir + "fields_scatters_pd_new.png"
	plt.savefig(savefile2 , format='png',dpi=200)
	plt.close()
	
	pd_X_new_main = pd.DataFrame(data = X_new[:,0:selected_pca], columns = new_field_list)
	scatter_matrix(pd_X_new_main,alpha=0.3,figsize=(12,7.27),diagonal='hist',hist_kwds={'color':'k', 'alpha':0.5, 'bins':50})
	savefile2 = output_dir + "fields_scatters_pd_new_main.png"
	plt.savefig(savefile2 , format='png',dpi=200)
	plt.close()
	
	pca_scatter(X_new,new_field_list,figurename,plot_kwds=plot_kwds)
	
	
	fig = plt.figure(1, figsize=(4, 3))
	#fig, axarr = plt.subplots(2, 2, sharex=True, figsize=(7.27, 9.69), dpi=200)
	plt.clf()
	
	return (figurename)

def pca_scatter(X_std,field_list,figurename, color=None,plot_kwds=plot_kwds):
	if ((len(field_list) - 1) % 2 == 0):
		fig3, axarr3 = plt.subplots(2, divmod((len(field_list)-1),2)[0], sharex=True, figsize=( 15,7.27), dpi=200)
		x=np.array(np.arange(len(field_list)-1).reshape(2, divmod((len(field_list)-1),2)[0]))
	else:
		fig3, axarr3 = plt.subplots(1, len(field_list)-1, sharex=True, figsize=( 15,7.27), dpi=200)
		x=np.array(np.arange(len(field_list)-1).reshape(1,len(field_list)-1))
	#print x
	for i,field in enumerate(field_list):
		if(i>0):
			#print  np.where(x==i-1)
			rows, columns = np.where(x==i-1)
			
			# print "rows, columns"
			# print rows
			# print columns
			
			if(color is not None):
				palette = sns.color_palette('deep', np.unique(color).max() + 1)
				colors = [palette[z] if z >= 0 else (0.0, 0.0, 0.0) for z in color]
				
				if ((len(field_list) - 1) % 2 == 0):axarr3[rows[0],columns[0]].scatter(X_std[:, 0], X_std[:,i], c=colors, alpha=0.3)
				else: axarr3[columns[0]].scatter(X_std[:, 0], X_std[:,i], c=colors, alpha=0.3)
			else:
				# print rows[0]
				# print columns[0]
				# print axarr3[rows[0]]
				# print X_std.shape
				if ((len(field_list) - 1) % 2 == 0):axarr3[rows[0],columns[0]].scatter(X_std[:, 0], X_std[:,i],**plot_kwds)#,alpha=0.3)
				else:axarr3[columns[0]].scatter(X_std[:, 0], X_std[:,i],**plot_kwds)#,alpha=0.3)
				
			
			if ((len(field_list) - 1) % 2 == 0):
				axarr3[rows[0],columns[0]].set_xlabel('%s'  %(field_list[0]))
				axarr3[rows[0],columns[0]].set_ylabel('%s'  %(field))
			else:
				axarr3[columns[0]].set_xlabel('%s'  %(field_list[0]))
				axarr3[columns[0]].set_ylabel('%s'  %(field))
	savefile3 = figurename
	fig3.tight_layout(h_pad=2.0)
	plt.savefig(savefile3 , format='png',dpi=200)
	plt.close()
	print ("figure %s is ok" %savefile3)


def clustered_scatter(pd_X_new,figurename,colors=None,labels_list=[],plot_kwds=plot_kwds,debug=True,centroids=None):
	"""
	Clusters are included in the pandas dataset pd_X_new. 
	If probability threshold (fuzzy clustering) is applied to the data -1 cluster is present.
	
	"""
	#print ("In clustered_scatter with {} file".format(figurename))
	if (colors is None):
		if((-1) in np.unique(pd_X_new["Clusters"])):
			#print ("(-1) in np.unique(pd_X_new['Clusters']")
			palette = sns.color_palette('deep',n_colors=(len(np.unique(pd_X_new["Clusters"]))-1))
		else:
			#print ("(-1) is not in np.unique(pd_X_new['Clusters']")
			palette = sns.color_palette('deep',n_colors=(len(np.unique(pd_X_new["Clusters"]))))
		# put grey color if the cluster is -1 use the color from the pallet otherwise
		colors = [palette[z] if z >= 0 else (0.8, 0.8, 0.8) for z in np.unique(pd_X_new["Clusters"])]
	print ("np.unique(pd_X_new['Clusters']) is {}".format(np.unique(pd_X_new["Clusters"])))
	print ("labels_list {}".format(labels_list))
		
	if (labels_list == []):
		for i in np.unique(pd_X_new["Clusters"]):
			if(i==-1):labels_list.append("---")
			else:labels_list.append('cl{}'.format(i+1))
			#pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': labels_list[-1]}, regex=False)
	print (np.unique(pd_X_new["Clusters"]))
	centroid_lables = []
	centroid_clusters=[]
	for i in np.unique(pd_X_new["Clusters"]):
		print 'Clusters i = {}'.format(i)
		print 'Clusters = {}'.format(labels_list[i])
		pd_X_new = pd_X_new.replace({'Clusters': i}, {'Clusters': labels_list[i]}, regex=False)
		centroid_lables.append("centroid_{}".format(labels_list[i]))
		centroid_clusters.append(i)
		
	#centroid_lables.append("Clusters")
	#print("hue_order {}".format(labels_list))
	#print ("np.unique(pd_X_new['Clusters']) is {}".format(np.unique(pd_X_new["Clusters"])))
	#print ("centroids.shape is {} and columns are {} ".format(centroids.shape,pd_X_new.columns))
	#pd_centroids = pd.DataFrame(centroids, columns=pd_X_new.columns[0:-1])
	#pd_centroids["label"] = centroid_lables
	#full_ds = pd.concat([pd_X_new, pd_centroids], ignore_index=True)
	
	
	#g = sns.pairplot(full_ds,diag_kind="hist", hue='Clusters', hue_order=labels_list.extend(centroid_lables), palette=colors.extend(colors), markers=".",diag_kws={'alpha':0.7},plot_kws=plot_kwds)
	#g = sns.pairplot(pd_X_new,diag_kind="hist", hue='Clusters', palette=colors,diag_kws={'alpha':0.7},plot_kws=plot_kwds)
	#centroid 
	if centroids is not None:
		centroids_array = np.c_[centroids,centroid_lables]#[0:-1]]
		print centroids_array
		pd_centroids = pd.DataFrame(centroids_array, columns=pd_X_new.columns)
		#pd_centroids["label"] = centroid_lables
		full_ds = pd.concat([pd_X_new, pd_centroids], ignore_index=True)
		markers = ["."]*len(centroid_lables)#[0:-1])
		markers.extend(["*"]*len(centroid_lables))#[0:-1]))
		size = np.repeat(20,len(centroid_lables))
		size = np.concatenate((size,np.repeat(100,len(centroid_lables))))
		print size
		#print centroid_lables[0:-1]
		#print labels_list+centroid_lables[0:-1]
		#print markers
		#print colors
		print colors+colors
	
	
		g = sns.PairGrid(full_ds, hue='Clusters', hue_order=labels_list+centroid_lables,palette=colors+colors, hue_kws={"s": size,"marker": markers})
		g.map(plt.scatter, linewidth=1, edgecolor="w")
		#g.map_offdiag(plt.scatter, linewidth=1, edgecolor="w")
		#g.map_diag(plt.hist, histtype="step", linewidth=3)
		g.add_legend()
		
		#g = sns.pairplot(pd_X_new,diag_kind="hist", hue='Clusters', hue_order=labels_list, palette=colors, markers=".",diag_kws={'alpha':0.7},plot_kws=plot_kwds)
	
		#g.data = pd_centroids
		#g.hue_vals = pd_centroids["Clusters"]
		#g.map_offdiag(plt.scatter, s=500, marker="*")
		
	else:
		g = sns.pairplot(pd_X_new,diag_kind="hist", hue='Clusters', hue_order=labels_list, palette=colors, markers=".",diag_kws={'alpha':0.7},plot_kws=plot_kwds)
	
		#print g.axes.flat
		#for ax in g.axes.flat: 
		#	plt.#setp(ax.get_xticklabels(), rotation=45)
	plt.savefig(figurename , format='png',dpi=200)
	plt.close()
	print ("figure {} is ok".format(figurename))
	return(colors,labels_list)
	

def t_contourplt(temp_data, x_lims, height_min, height_max, ax=None):
	if( ax is None):
	    cp = plt.contour(temp_data,
	                aspect='auto',origin='lower',
	                extent =(x_lims[0], x_lims[1],height_min, height_max),
	                linewidths=0.5,
	                #interpolation='none',
	                levels=[-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40],colors='k')#['k','k','k','k','k','k','k','k','k'])
	                #levels=[0,10,20],colors=['k','k','k'])
	    plt.clabel(cp, inline=True,fontsize=10,fmt = '%2.0f$^\circ$C')
	else:
		cp = ax.contour(temp_data,
	                aspect='auto',origin='lower',
	                extent =(x_lims[0], x_lims[1],height_min, height_max),
	                linewidths=0.5,
	                #interpolation='none',
	                levels=[-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40],colors='k')#['k','k','k','k','k','k','k','k','k'])
	                #levels=[0,10,20],colors=['k','k','k'])
		ax.clabel(cp, inline=True,fontsize=8,fmt = '%2.0f$^\circ$C')


def z_contourplt(qvp_file,count_threshold, x_lims, ax=None):
	if( ax is None):
	    cp = plt.contour(np.where(qvp_file['dBuZ/Counts'][:]>count_threshold,
	                    qvp_file['dBuZ/Means'][:],
	                    np.nan),
	                aspect='auto',origin='lower',
	                extent =(x_lims[0], x_lims[1],qvp_file['Height'][0],qvp_file['Height'][-1]),
	                #interpolation='none',
	                levels=[0,8,16,24,32,40],colors=['k','k','k','k','k','k'])
	                #levels=[0,10,20],colors=['k','k','k'])
	    plt.clabel(cp, inline=False,fontsize=10,fmt = '%2.0f')
	else:
		cp = ax.contour(np.where(qvp_file['dBuZ/Counts'][:]>count_threshold,
	                    qvp_file['dBuZ/Means'][:],
	                    np.nan),
	                aspect='auto',origin='lower',
	                extent =(x_lims[0], x_lims[1],qvp_file['Height'][0],qvp_file['Height'][-1]),
	                #interpolation='none',
	                levels=[0,8,16,24,32,40],colors=['k','k','k','k','k','k'])
	                #levels=[0,10,20],colors=['k','k','k'])
		ax.clabel(cp, inline=False,fontsize=10,fmt = '%2.0f')

def plot_faam_data(faam_file,field,visib=None):
	#print visib
	#line = plt.plot(faam_file['Time'][:], faam_file[field][:], color='black')#'0.75')
	if (visib is not None):
		for i,time_int in enumerate(visib):
			print "i = %d   " %(i)
			print "time_int"
			print time_int
			print faam_file['Time'][:]
			#plt.plot(faam_file['Time'][:], faam_file[field][:], color='black')
	
def plot_standard_4(qvp_file,timeticks,output_dir,case_id=None,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,verbose=False,sys_phi_list=[],night_mask=None):
    
    ## QVP xlims
    start_time = num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
    end_time = num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
    x_lims = mdates.date2num([start_time,end_time])
    time_shape = len(qvp_file['Time'])
    alt_shape = len(qvp_file['Height'])
    T = np.empty([time_shape*alt_shape], dtype=float)*np.nan
    #print T.shape
    #print np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	#                    np.nan)
    #T = np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	#                    np.nan)
	                 
    height_min = qvp_file['Height'][0]
    height_max = qvp_file['Height'][-1]
    
    mean_16colors = ("#4A6FE3","#657FE1","#7B8EE1","#8F9DE1","#A2ACE2","#B5BBE3","#C7CBE3","#D9DAE3","#E4D8DA","#E6C4C9","#E6AFB9","#E49AA9","#E28699","#DE7089","#D95979","#D33F6A")
    std_16colors = ("#11C638","#4ACA59","#6ACF73","#83D38A","#9AD79F","#B0DAB3","#C4DEC6","#D8E1D9","#E5DDD9","#EAD3C5","#EDC9B0","#EFBF9A","#F1B583","#F1AB69","#F1A149","#EF9708")
    
    zdr_colors = sns.color_palette("hls", 15)
    zdrmap,zdrnorm = leveled_discreet_cmap([-1.5,-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],zdr_colors)
   
   
    rhohv_colors = sns.color_palette("hls", 11)
    rhomap,rhonorm = leveled_discreet_cmap([0.7,0.8,0.9,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1],rhohv_colors)
    
    kdp_colors = sns.color_palette("hls", 15)
    kdpmap,kdpnorm = leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
    
    phi_colors = sns.color_palette("hls", 16)
    phimap,phinorm =leveled_discreet_cmap([-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10],phi_colors)
    
    z_colors = sns.color_palette("hls", 12)
    Zmap,Znorm = leveled_discreet_cmap([-20,-10,0,4,8,12,16,20,24,28,32,36,40],z_colors)
    
    v_colors = sns.color_palette("hls", 16)
    Vmap,Vnorm = leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
    
    
    temp_colors = sns.color_palette("hls", 15)
    tempmap,tempnorm = leveled_discreet_cmap(np.arange(10,-20,-2),temp_colors)
	
    
    #sys_phi = 109
    zdr_cal = 0#+2.1
    
    #date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    date_format = mdates.DateFormatter('%b %d %Y\n%H:%M')
    
    fig, ax_arr = plt.subplots(2, 2, sharey=False,figsize=(20,10))
    #fig.suptitle("Dual-pol observations on {}".format(start_time.strftime('%Y/%m/%d')), size=16)
    fig.subplots_adjust(top=0.7)

    #fig = plt.figure(figsize=(20,14))
    
    if x1>len(timeticks)-1:
        x1 = -1
    
    ## Plot Reflectivity - Subplot 1 ----------------------------------------------
    #plt.subplot(321)
    
    
    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'dBuZ',Znorm,Zmap,x_lims,-20,40,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
    annotate_field_plot(r'$\mathrm{Z_{H}}$ (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
    
    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])
    #t_contourplt(temp_data, x_lims, height_min, height_max, ax=None)
    
    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-1.5,2,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0],night_mask=night_mask)
    annotate_field_plot(r'$\mathrm{Z_{DR}}$ (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
    
    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,0])

    ## Plot RhoHV - Subplot 2 ----------------------------------------------------
    #plt.subplot(322)
    
    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    #im = plot_field(qvp_file,'RhoHVu',rhonorm,rhomap,x_lims,0.7,1,count_threshold=count_thr,ax=ax_arr[0,1])
    
    im = plot_field(qvp_file,'RhoHV',rhonorm,rhomap,x_lims,0.7,1,count_threshold=count_thr,ax=ax_arr[0,1],night_mask=night_mask)
    annotate_field_plot(r'$\mathrm{\rho_{HV}}$',rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
    
    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[0,1])
    

	## Plot KDP - Subplot 4 ----------------------------------------------------
	#plt.subplot(324)

    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    #im = plot_field(qvp_file,'PhiDP',phinorm,phimap,x_lims,-2,12,offset=-sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1])
    im = plot_field(qvp_file,'uPhiDP',phinorm,phimap,x_lims,-10,15,offset= -sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
    
    annotate_field_plot(r'Phase shift from signal processor (deg)',phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
    
    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])
    
    plt.tight_layout()
    #if savefile:
    if not os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d'))):os.makedirs("{}{}".format(output_dir,start_time.strftime('%Y%m%d')))
	
    savefilename = '{}{}/case{}_qvp_4var_plot.png'.format(output_dir,start_time.strftime('%Y%m%d'),case_id)
    plt.savefig(savefilename , format='png',dpi=200)
    if verbose:print savefilename
    #plt.show()
    plt.close()


def plot_first_4(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=None,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,verbose=False,sys_phi_list=[],night_mask=None):
    print "in plot_first_4"
    print "night_mask"
    print night_mask
    ## QVP xlims
    start_time = num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
    end_time = num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
    x_lims = mdates.date2num([start_time,end_time])
    time_shape = len(qvp_file['Time'])
    alt_shape = len(qvp_file['Height'])
    T = np.empty([time_shape*alt_shape], dtype=float)*np.nan
    
    #print "In plot_first_4 T.shape is"
    #print T.shape
    
    #print np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	#                    np.nan)
   # T = np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	 #                   np.nan)
	                 
    height_min = qvp_file['Height'][0]
    height_max = qvp_file['Height'][-1]
    
    mean_16colors = ("#4A6FE3","#657FE1","#7B8EE1","#8F9DE1","#A2ACE2","#B5BBE3","#C7CBE3","#D9DAE3","#E4D8DA","#E6C4C9","#E6AFB9","#E49AA9","#E28699","#DE7089","#D95979","#D33F6A")
    std_16colors = ("#11C638","#4ACA59","#6ACF73","#83D38A","#9AD79F","#B0DAB3","#C4DEC6","#D8E1D9","#E5DDD9","#EAD3C5","#EDC9B0","#EFBF9A","#F1B583","#F1AB69","#F1A149","#EF9708")
    
    
    zdr_colors = sns.color_palette("hls", 20)
    zdrmap,zdrnorm = leveled_discreet_cmap([-2,-1.5,-1.0,-0.5,0,0.5,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0],zdr_colors)
	
   
   
    rhohv_colors = sns.color_palette("hls", 12)
    rhomap,rhonorm =leveled_discreet_cmap([0.6,0.63,0.66,0.70,0.73,0.76,0.80,0.83,0.86,0.9,0.93,0.96,1.0],rhohv_colors)#
	
    kdp_colors = sns.color_palette("hls", 15)
    kdpmap,kdpnorm = leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
	
    phi_colors = sns.color_palette("hls", 11)
    phimap,phinorm =leveled_discreet_cmap([110,115,120,125,130,135,140,145,150,155,160],phi_colors)
	
    z_colors = sns.color_palette("hls", 10)
    Zmap,Znorm = leveled_discreet_cmap([-15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15],z_colors)
	
    v_colors = sns.color_palette("hls", 16)
    Vmap,Vnorm = leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
    
    
    temp_colors = sns.color_palette("hls", 15)
    tempmap,tempnorm = leveled_discreet_cmap(np.arange(10,-20,-2),temp_colors)
	
    
    #sys_phi = 109
    zdr_cal = 0#+2.1
    
    #date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    #date_format = mdates.DateFormatter('%b %d\n%H:%M')
    date_format = mdates.DateFormatter('%m-%d-%Y\n%H:%M')
    #date_format = mdates.DateFormatter('%H:%M')
    
    fig, ax_arr = plt.subplots(2, 2,sharex=True, sharey=True,figsize=(20,10))
    #fig.suptitle("Dual-pol observations on {}".format(start_time.strftime('%Y/%m/%d')), size=16)
    fig.subplots_adjust(top=0.88)

    #fig = plt.figure(figsize=(20,14))
    
    if x1>len(timeticks)-1:x1 = -1
    
    ## Plot Reflectivity - Subplot 1 ----------------------------------------------
    if('dBZ' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'dBZ',Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{Z_{H}}$ (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])#'Reflectivity (dBZ)'
    elif('dBZv' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'dBZv',Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{Z_{V}}$ (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])#'Reflectivity (dBZ)'
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[0],Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[0],Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])
    
    if('ZDR' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-2,8,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0],night_mask=night_mask)
	    annotate_field_plot( r'$\mathbf{Z_{DR}}$ (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,0])#'ZDR
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[1],zdrnorm,zdrmap,x_lims,-2,8,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[1],zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,0])
	
    ## Plot RhoHV - Subplot 2 ----------------------------------------------------
 
    if('RhoHV' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'RhoHV',rhonorm,rhomap,x_lims,0.6,1,count_threshold=count_thr,ax=ax_arr[0,1],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{\rho_{HV}}$',rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[0,1])
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,qvp_file,field_list[2],rhonorm,rhomap,x_lims,0.6,1,count_threshold=count_thr,ax=ax_arr[0,1],night_mask=night_mask)
	    annotate_field_plot(qvp_file,field_lables_list[2],rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[0,1])
    

	## Plot KDP - Subplot 4 ----------------------------------------------------
    if('KDP' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'KDP_UKMO',kdpnorm,kdpmap,x_lims,-0.3,0.5,count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{K_{DP} (^\circ km^{-1})}$',kdpnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])#[0,1])
    elif('PhiDP' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'uPhiDP',phinorm,phimap,x_lims,-10,15,offset= -sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(r'Phase shift from signal processor (deg)',phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[3],phinorm,phimap,x_lims,-10,15,offset= -sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[3],phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])
      
    #fig.suptitle(start_time.strftime('%m-%d-%Y'), fontsize=16, fontweight = 'bold')
    plt.tight_layout()
    #if savefile:
    if not os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d'))):
		print "path exist is {}".format(os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d'))))
		print "{}{}".format(output_dir,start_time.strftime('%Y%m%d'))
		os.makedirs("{}{}".format(output_dir,start_time.strftime('%Y%m%d')))
    
    print "writing to"
    print "{}{}".format(output_dir,start_time.strftime('%Y%m%d'))
    print os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d')))
    savefilename = '{}{}/case{}_qvp_first4var_plot.png'.format(output_dir,start_time.strftime('%Y%m%d'),case_id)
    print savefilename
    plt.savefig(savefilename , format='png',dpi=300)
    if verbose:print savefilename
    #plt.show()
    plt.close()

def plot_first_4_in_timeseries(qvp_file,timeticks,output_dir,field_lables_list,field_list,case_id=None,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None,verbose=False,sys_phi_list=[],night_mask=None):
    print "in plot_first_4in_timeseries"
    print "night_mask"
    print night_mask
    ## QVP xlims
    start_time = num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
    end_time = num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
    x_lims = mdates.date2num([start_time,end_time])
    time_shape = len(qvp_file['Time'])
    alt_shape = len(qvp_file['Height'])
    T = np.empty([time_shape*alt_shape], dtype=float)*np.nan
    
    #print "In plot_first_4 T.shape is"
    #print T.shape
    
    #print np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	#                    np.nan)
   # T = np.where(qvp_file['temperature_2/Counts'][:]>count_thr,
	#                    (qvp_file['temperature_2/Means'][:] - 273.15),
	 #                   np.nan)
	                 
    height_min = qvp_file['Height'][0]
    height_max = qvp_file['Height'][-1]
    
    mean_16colors = ("#4A6FE3","#657FE1","#7B8EE1","#8F9DE1","#A2ACE2","#B5BBE3","#C7CBE3","#D9DAE3","#E4D8DA","#E6C4C9","#E6AFB9","#E49AA9","#E28699","#DE7089","#D95979","#D33F6A")
    std_16colors = ("#11C638","#4ACA59","#6ACF73","#83D38A","#9AD79F","#B0DAB3","#C4DEC6","#D8E1D9","#E5DDD9","#EAD3C5","#EDC9B0","#EFBF9A","#F1B583","#F1AB69","#F1A149","#EF9708")
    
    
    zdr_colors = sns.color_palette("hls", 20)
    zdrmap,zdrnorm = leveled_discreet_cmap([-2,-1.5,-1.0,-0.5,0,0.5,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0],zdr_colors)
	
   
   
    rhohv_colors = sns.color_palette("hls", 12)
    rhomap,rhonorm =leveled_discreet_cmap([0.6,0.63,0.66,0.70,0.73,0.76,0.80,0.83,0.86,0.9,0.93,0.96,1.0],rhohv_colors)#
	
    kdp_colors = sns.color_palette("hls", 15)
    kdpmap,kdpnorm = leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
	
    phi_colors = sns.color_palette("hls", 11)
    phimap,phinorm =leveled_discreet_cmap([110,115,120,125,130,135,140,145,150,155,160],phi_colors)
	
    z_colors = sns.color_palette("hls", 10)
    Zmap,Znorm = leveled_discreet_cmap([-15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15],z_colors)
	
    v_colors = sns.color_palette("hls", 16)
    Vmap,Vnorm = leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
    
    
    temp_colors = sns.color_palette("hls", 15)
    tempmap,tempnorm = leveled_discreet_cmap(np.arange(10,-20,-2),temp_colors)
	
    
    #sys_phi = 109
    zdr_cal = 0#+2.1
    
    #date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    #date_format = mdates.DateFormatter('%b %d\n%H:%M')
    date_format = mdates.DateFormatter('%m-%d-%Y\n%H:%M')
    #date_format = mdates.DateFormatter('%H:%M')
    
    fig, ax_arr = plt.subplots(2, 2,sharex=True, sharey=True,figsize=(20,10))
    #fig.suptitle("Dual-pol observations on {}".format(start_time.strftime('%Y/%m/%d')), size=16)
    fig.subplots_adjust(top=0.88)

    #fig = plt.figure(figsize=(20,14))
    
    if x1>len(timeticks)-1:x1 = -1
    ## Plot Reflectivity - Subplot 1 ----------------------------------------------
    if('dBZ' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'dBZ',Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{Z_{H}}$ (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])#'Reflectivity (dBZ)'
    elif('dBZv' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'dBZv',Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{Z_{V}}$ (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])#'Reflectivity (dBZ)'
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[0],Znorm,Zmap,x_lims,-15,15,count_threshold=count_thr,ax=ax_arr[0,0],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[0],Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
	    #t_contourplt(T, x_lims, height_min, height_max,ax=ax_arr[0,0])
    
    if('ZDR' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-2,8,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0],night_mask=night_mask)
	    annotate_field_plot( r'$\mathbf{Z_{DR}}$ (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,0])#'ZDR
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[1],zdrnorm,zdrmap,x_lims,-2,8,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[1],zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,0])
	
    ## Plot RhoHV - Subplot 2 ----------------------------------------------------
 
    if('RhoHV' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'RhoHV',rhonorm,rhomap,x_lims,0.6,1,count_threshold=count_thr,ax=ax_arr[0,1],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{\rho_{HV}}$',rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[0,1])
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,qvp_file,field_list[2],rhonorm,rhomap,x_lims,0.6,1,count_threshold=count_thr,ax=ax_arr[0,1],night_mask=night_mask)
	    annotate_field_plot(qvp_file,field_lables_list[2],rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[0,1])
    

	## Plot KDP - Subplot 4 ----------------------------------------------------
    if('KDP' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'KDP_UKMO',kdpnorm,kdpmap,x_lims,-0.3,0.5,count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(r'$\mathbf{K_{DP} (^\circ km^{-1})}$',kdpnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])#[0,1])
    elif('PhiDP' in field_lables_list):
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'uPhiDP',phinorm,phimap,x_lims,-10,15,offset= -sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(r'Phase shift from signal processor (deg)',phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])
    else:
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,field_list[3],phinorm,phimap,x_lims,-10,15,offset= -sys_phi_list[-1], count_threshold=count_thr,ax=ax_arr[1,1],night_mask=night_mask)
	    annotate_field_plot(field_lables_list[3],phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    #t_contourplt(T, x_lims, height_min,height_max,ax=ax_arr[1,1])
      
    #fig.suptitle(start_time.strftime('%m-%d-%Y'), fontsize=16, fontweight = 'bold')
    plt.tight_layout()
    #if savefile:
    if not os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d'))):
		print "path exist is {}".format(os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d'))))
		print "{}{}".format(output_dir,start_time.strftime('%Y%m%d'))
		os.makedirs("{}{}".format(output_dir,start_time.strftime('%Y%m%d')))
    
    print "writing to"
    print "{}{}".format(output_dir,start_time.strftime('%Y%m%d'))
    print os.path.exists("{}{}".format(output_dir,start_time.strftime('%Y%m%d')))
    savefilename = '{}{}/case{}_qvp_first4var_plot.png'.format(output_dir,start_time.strftime('%Y%m%d'),case_id)
    print savefilename
    plt.savefig(savefilename , format='png',dpi=300)
    if verbose:print savefilename
    #plt.show()
    plt.close()


def Plot_standard_6(qvp_file,ver_file,timeticks,output_dir,faam_file=None,savefile=False,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None):
    visib=None
    ## QVP xlims
    start_time = num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
    end_time = num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
    x_lims = mdates.date2num([start_time,end_time])
    if (ver_file is not None):
		## Vertical xlims
		vstart_time = num2date(ver_file['Time'][0],ver_file['Time'].units,qvp_file['Time'].calendar)
		vend_time = num2date(ver_file['Time'][-1],ver_file['Time'].units,qvp_file['Time'].calendar)
		vx_lims = mdates.date2num([vstart_time,vend_time])
    
    
    mean_16colors = ("#4A6FE3","#657FE1","#7B8EE1","#8F9DE1","#A2ACE2","#B5BBE3","#C7CBE3","#D9DAE3","#E4D8DA","#E6C4C9","#E6AFB9","#E49AA9","#E28699","#DE7089","#D95979","#D33F6A")
    std_16colors = ("#11C638","#4ACA59","#6ACF73","#83D38A","#9AD79F","#B0DAB3","#C4DEC6","#D8E1D9","#E5DDD9","#EAD3C5","#EDC9B0","#EFBF9A","#F1B583","#F1AB69","#F1A149","#EF9708")
    
    zdr_colors = sns.color_palette("hls", 15)
    zdrmap,zdrnorm = leveled_discreet_cmap([-1.5,-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],zdr_colors)
   
   
    rhohv_colors = sns.color_palette("hls", 11)
    rhomap,rhonorm = leveled_discreet_cmap([0.7,0.8,0.9,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1],rhohv_colors)
    
    kdp_colors = sns.color_palette("hls", 15)
    kdpmap,kdpnorm = leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
    
    phi_colors = sns.color_palette("hls", 11)
    phimap,phinorm =leveled_discreet_cmap([-2,-1,0,1,2,3,4,5,6,8,10,12],phi_colors)
    
    z_colors = sns.color_palette("hls", 12)
    Zmap,Znorm = leveled_discreet_cmap([-20,-10,0,4,8,12,16,20,24,28,32,36,40],z_colors)
    
    v_colors = sns.color_palette("hls", 16)
    Vmap,Vnorm = leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
    
    
    temp_colors = sns.color_palette("hls", 15)
    tempmap,tempnorm = leveled_discreet_cmap(np.arange(10,-20,-2),temp_colors)
	
    
    sys_phi = 109
    zdr_cal = 0#+2.1
    
    #date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    date_format = mdates.DateFormatter('%b %d\n%H:%M')
    
    fig, ax_arr = plt.subplots(3, 2, sharey=False,figsize=(20,14))
    #fig = plt.figure(figsize=(20,14))
    
    if x1>len(timeticks)-1:
        x1 = -1
    
    ## Plot Reflectivity - Subplot 1 ----------------------------------------------
    #plt.subplot(321)
    
    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'dBuZ',Znorm,Zmap,x_lims,-20,40,count_threshold=count_thr,ax=ax_arr[0,0])
    annotate_field_plot('Reflectivity (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[0,0])
    #if(faam_file != None):
    plot_faam_data(faam_file,'ALT_GIN',visib=visib)

    # ## Plot ZDR - Subplot 3 -------------------------------------------------------
    # plt.subplot(323)
    # plot_field(qvp_file,'ZDRu',zdrnorm,zdrmap,x_lims,-1,3,offset=zdr_cal,count_threshold=count_thr)
    # annotate_field_plot('ZDR (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin)
    # z_contourplt(qvp_file, count_thr, x_lims)
    
    ## Plot ZDR - Subplot 3 -------------------------------------------------------
    #plt.subplot(323)
    
    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-2,8,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0])
    annotate_field_plot('ZDR (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[1,0])

    ## Plot RhoHV - Subplot 2 ----------------------------------------------------
    #plt.subplot(322)
    
    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'RhoHVu',rhonorm,rhomap,x_lims,0.6,1,count_threshold=count_thr,ax=ax_arr[0,1])
    annotate_field_plot('CC',rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[0,1])
    
    
    if (ver_file is not None):

		## Plot KDP - Subplot 4 ----------------------------------------------------
		#plt.subplot(324)
    
	    if(hs is not None):hs.plot(ax=ax_arr[1,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[1,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(qvp_file,'KDP_UKMO',kdpnorm,kdpmap,x_lims,-0.3,0.5,count_threshold=count_thr,ax=ax_arr[1,1])
	    annotate_field_plot(r'$\mathrm{K_{DP}}$ (deg/km)',kdpnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,1], im=im)
	    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[1,1])
	
	    ## Plot Fall speed - Subplot 5 ----------------------------------------------------
	    #plt.subplot(325)
	    if(hs is not None):hs.plot(ax=ax_arr[2,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
	    if(hvis is not None):hvis.plot(ax=ax_arr[2,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
	    im = plot_field(ver_file,'Vu',Vnorm,Vmap,vx_lims,-0.5,6,scaler=-1,count_threshold=count_thr,ax=ax_arr[2,0])
	    annotate_field_plot('Fall speed unfiltered (m/s)',Vnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[2,0], im=im)
	    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[2,0])
	
    ## Plot uPhiDP - Subplot 6 ----------------------------------------------------
    #plt.subplot(326)
    
    if(hs is not None):hs.plot(ax=ax_arr[2,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[2,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_field(qvp_file,'uPhiDP',phinorm,phimap,x_lims,-2,12,offset= -sys_phi, count_threshold=count_thr,ax=ax_arr[2,1])
    annotate_field_plot(r'Phase shift from signal processor (deg)',phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[2,1], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[2,1])


    plt.tight_layout()
    if savefile:
		savefilename = '%sqvp_plot.png' %output_dir
		plt.savefig(savefilename , format='png',dpi=200)
    #plt.show()
    plt.close()


def Plot_standard_6_stddev(qvp_file,ver_file,timeticks,output_dir,faam_file=None,savefile=False,hmax=10000,hmin=0,x0=0,x1=-1,count_thr=220,hs=None,hvis=None):
	
    ## QVP xlims
    start_time = num2date(qvp_file['Time'][0],qvp_file['Time'].units,qvp_file['Time'].calendar)
    end_time = num2date(qvp_file['Time'][-1],qvp_file['Time'].units,qvp_file['Time'].calendar)
    x_lims = mdates.date2num([start_time,end_time])
    
    if (ver_file is not None):
		## Vertical xlims
		vstart_time = num2date(ver_file['Time'][0],ver_file['Time'].units,qvp_file['Time'].calendar)
		vend_time = num2date(ver_file['Time'][-1],ver_file['Time'].units,qvp_file['Time'].calendar)
		vx_lims = mdates.date2num([vstart_time,vend_time])
    
    
    mean_16colors = ("#4A6FE3","#657FE1","#7B8EE1","#8F9DE1","#A2ACE2","#B5BBE3","#C7CBE3","#D9DAE3","#E4D8DA","#E6C4C9","#E6AFB9","#E49AA9","#E28699","#DE7089","#D95979","#D33F6A")
    std_16colors = ("#11C638","#4ACA59","#6ACF73","#83D38A","#9AD79F","#B0DAB3","#C4DEC6","#D8E1D9","#E5DDD9","#EAD3C5","#EDC9B0","#EFBF9A","#F1B583","#F1AB69","#F1A149","#EF9708")
    
    #zdr_colors = sns.color_palette("hls", 15)
    #zdrmap,zdrnorm = leveled_discreet_cmap([-1.5,-1,-0.5,0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],zdr_colors)
    std_zdr_colors = sns.color_palette("hls", 12)
    zdrmap,zdrnorm = leveled_discreet_cmap([0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2],std_zdr_colors)
	
    rhohv_colors = sns.color_palette("hls", 11)
	#rhomap,rhonorm = qvpp.leveled_discreet_cmap([0.7,0.8,0.9,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.995,1],rhohv_colors)
    rhomap,rhonorm = leveled_discreet_cmap([0,0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.22,0.24,0.26,0.28],rhohv_colors)
	
    
    kdp_colors = sns.color_palette("hls", 15)
    std_kdp_colors = sns.color_palette("hls", 12)
	#kdpmap,kdpnorm = qvpp.leveled_discreet_cmap([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6],kdp_colors)
    kdpmap,kdpnorm = leveled_discreet_cmap([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1],std_kdp_colors)
	
    phi_colors = sns.color_palette("hls", 10)
    #phimap,phinorm =leveled_discreet_cmap([-2,-1,0,1,2,3,4,5,6,8,10,12],phi_colors)
    phimap,phinorm =leveled_discreet_cmap([0,1,2,3,4,5,6,7,8,9,10],phi_colors)
    
    z_colors = sns.color_palette("hls", 12)
    Zmap,Znorm = leveled_discreet_cmap([-20,-10,0,4,8,12,16,20,24,28,32,36,40],z_colors)
    
    v_colors = sns.color_palette("hls", 16)
    Vmap,Vnorm = leveled_discreet_cmap([-4,0,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4,3.0,4.0,6.0,8.0],v_colors)
    
    sys_phi = 109
    zdr_cal = 0#+2.1
    
    #date_format = mdates.DateFormatter('%b %d\n%H:%M:%S')
    date_format = mdates.DateFormatter('%b %d\n%H:%M')
    
    fig, ax_arr = plt.subplots(3, 2, sharey=False,figsize=(20,14))
    #fig = plt.figure(figsize=(20,14))
    
    if x1>len(timeticks)-1:
        x1 = -1
    
    ## Plot Reflectivity - Subplot 1 ----------------------------------------------
    #plt.subplot(321)
    if(hs is not None):hs.plot(ax=ax_arr[0,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_stddev_field(qvp_file,'dBuZ',Znorm,Zmap,x_lims,-20,40,count_threshold=count_thr,ax=ax_arr[0,0])
    #plot_field(qvp_file,'dBuZ',Znorm,Zmap,x_lims,-20,40,count_threshold=count_thr,ax=ax_arr[0,0])
    annotate_field_plot('StdDev of Reflectivity (dBZ)',Znorm,hmax,date_format,timeticks,x0=x0,x1=x1,ax=ax_arr[0,0], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[0,0])
    if(faam_file != None):plot_faam_data(faam_file,'ALT_GIN')

    # ## Plot ZDR - Subplot 3 -------------------------------------------------------
    # plt.subplot(323)
    # plot_field(qvp_file,'ZDRu',zdrnorm,zdrmap,x_lims,-1,3,offset=zdr_cal,count_threshold=count_thr)
    # annotate_field_plot('ZDR (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin)
    # z_contourplt(qvp_file, count_thr, x_lims)
    
    ## Plot ZDR - Subplot 3 -------------------------------------------------------
    #plt.subplot(323)
    
    if(hs is not None):hs.plot(ax=ax_arr[1,0],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[1,0],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_stddev_field(qvp_file,'ZDR',zdrnorm,zdrmap,x_lims,-1.5,2,offset=zdr_cal,count_threshold=count_thr,ax=ax_arr[1,0])
    annotate_field_plot('StdDev of ZDR (dB)',zdrnorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[1,0], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[1,0])

    ## Plot RhoHV - Subplot 2 ----------------------------------------------------
    #plt.subplot(322)
    
    if(hs is not None):hs.plot(ax=ax_arr[0,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[0,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_stddev_field(qvp_file,'RhoHVu',rhonorm,rhomap,x_lims,0.,0.3,count_threshold=count_thr,ax=ax_arr[0,1])
    annotate_field_plot('StdDev of CC',rhonorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[0,1], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[0,1])

    
		## Plot KDP - Subplot 4 ----------------------------------------------------
	  
    ## Plot uPhiDP - Subplot 6 ----------------------------------------------------
    #plt.subplot(326)
    
    if(hs is not None):hs.plot(ax=ax_arr[2,1],color = '0.25', alpha=0.5, lw=3,label='FAAM flight height', legend=True)
    if(hvis is not None):hvis.plot(ax=ax_arr[2,1],color = 'black', style='*', lw=10,label='FAAM is inside QVP coverage', legend=True)
    im = plot_stddev_field(qvp_file,'uPhiDP',phinorm,phimap,x_lims,0, 10,count_threshold=count_thr,ax=ax_arr[2,1])
    annotate_field_plot(r'StdDev of Phase shift from signal processor (deg)',phinorm,hmax,date_format,timeticks,x0=x0,x1=x1,hmin=hmin,ax=ax_arr[2,1], im=im)
    z_contourplt(qvp_file, count_thr, x_lims,ax=ax_arr[2,1])


    plt.tight_layout()
    if savefile:
		savefilename = '%sqvp_stddev_plot.png' %output_dir
		plt.savefig(savefilename , format='png',dpi=200)
    #plt.show()
    plt.close()

