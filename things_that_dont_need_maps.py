#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:54:34 2020

@author: margaret

Num Waves, Seasonality, and other things that don't need to be plotted on a map
also plots things that don't need maps, like max intensity and wave distance

9 March 2021: set a color palette so I don't have to keep doing that manually

18 November 2021: added printout of text for the mean and standard deviation, changed boxplots to violinplots

5 September 2022: changed file paths for updated files (newly added years, updated TRACK to TEW filtering)
    changed ifs to whiles for fringe cases where there are multiple zero-length tracks at the end of a file
    
16 March 2023:
- all the following changes were made to adhere to Clim. Dyn. publishing requirements
- changed fig sizes
- removed fig titles
- added a/b part labels
- changed output path to PublicationVersions folder
- increased resolution to 600 dpi
- changed output file type to tiff
"""
# # # # # # # # #
# import packages
# # # # # # # # #

#general stuff
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

# netCDF-specific stuff
import netCDF4
from netCDF4 import Dataset, num2date

# # # # # # # # # # # # #
# set up things for stuff
# # # # # # # # # # # # #

level = int(sys.argv[1])

#list of files
NHfiles = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.NH.NoTC.nc'))
#NHfiles = sorted(glob.glob(f'/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredNoTC/{level}hPa/post_filter.*.NH.NoTC.nc'))
SHfiles = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.SH.NoTC.nc'))
#SHfiles = sorted(glob.glob(f'/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredNoTC/{level}hPa/post_filter.*.SH.NoTC.nc'))
#NHf = '/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/NewFiles/TRACK_netCDF/850hPa/post_filter.1981.nc'
#SHf = '/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/NewFiles/TRACK_netCDF/850hPa/post_filter.1981.SH.nc'

# =============================================================================
# NHfiles = []
# SHfiles = []
# for i in range(20):
#     NHfiles.append(NHf)
#     SHfiles.append(SHf)
# =============================================================================

#will want some arrays to hold relevant information
#global data
peryear= []
perDJF = []
perMAM = []
perJJA = []
perSON = []

#NH Data
peryear_NH = []
perDJF_NH = []
perMAM_NH = []
perJJA_NH = []
perSON_NH = []

#SH Data
peryear_SH = []
perDJF_SH = []
perMAM_SH = []
perJJA_SH = []
perSON_SH = []

#Stuff without a seasonal breakdown:
    #max int
    #wave distance
    #wave speed
max_int  = []
max_int_NH = []
max_int_SH = []

net_dist = []
net_dist_NH = []
net_dist_SH = []

avg_spd  = []
avg_spd_NH  = []
avg_spd_SH  = []

def seasons(dates):
    '''
    A function for getting seasonal counts for data that spans a year
    ----------
    
    Parameters
    ----------
    dates : cf date

    Returns
    -------
    DJF : count of things in Dec, Jan, Feb
    MAM : count of things in Mar, Apr, May
    JJA : count of things in Jun, Jul, Aug
    SON : count of thigns in Sep, Oct, Nov

    '''
    DJF = []
    MAM = []
    JJA = []
    SON = []
    for date in dates:
        if 3 <= date.month <= 5:
            MAM.append(date)
        elif 6 <= date.month <= 8:
            JJA.append(date)
        elif 9 <= date.month <= 11:
            SON.append(date)
        else:
            DJF.append(date)
    
    return DJF,MAM,JJA,SON
            
# Get a count of the tracks per year, both globally and per hemisphere
# Use the seasons function to get tracks per season, globally and per hemisphere
for NHfile,SHfile in zip(NHfiles,SHfiles):
    #might help to open the file first
    NHtracks = Dataset(NHfile, mode='r')
    SHtracks = Dataset(SHfile, mode='r')
       
    #I DO NOT WANT ALL OF THE TRACKS I WANT TO ONLY READ IN THE RELEVANT THINGS

    #things not getting seasonal breakdown
    #max intensity
    max_int.extend(NHtracks['max_int'][:])
    max_int.extend(SHtracks['max_int'][:])
    max_int_NH.extend(NHtracks['max_int'][:])
    max_int_SH.extend(SHtracks['max_int'][:])
    
    #nan the data for tracks that matched from their first points
    for i in range(len(max_int_NH)):
        if max_int_NH[i] == 0.0:
            max_int_NH[i] = np.NaN

    for i in range(len(max_int_SH)):
        if max_int_SH[i] == 0.0:
            max_int_SH[i] = np.NaN
            
    for i in range(len(max_int)):
        if max_int[i] == 0.0:
            max_int[i] = np.NaN
        
    #total distance traveled
    net_dist.extend(NHtracks['net_dist_2d'][:])
    net_dist.extend(SHtracks['net_dist_2d'][:])
    net_dist_NH.extend(NHtracks['net_dist_2d'][:])
    net_dist_SH.extend(SHtracks['net_dist_2d'][:])
    
    #average wave speed
    avg_spd.extend(NHtracks['avg_spd'][:])
    avg_spd.extend(SHtracks['avg_spd'][:])
    avg_spd_NH.extend(NHtracks['avg_spd'][:])
    avg_spd_SH.extend(SHtracks['avg_spd'][:])
    
    #things for seasonal breakdowns
    #NH stuff
    NHfirst = NHtracks.variables['FIRST_PT'][:]
    while NHfirst[-1] == len(NHtracks['longitude']):
        NHfirst = NHfirst[0:-1]
    
    NHtimes = NHtracks.variables['time'][NHfirst]
    NHdateunits = NHtracks.variables['time'].units
    NHdates = num2date(NHtimes, NHdateunits)
    
    #SH stuff
    SHfirst = SHtracks.variables['FIRST_PT'][:]
    while SHfirst[-1] == len(SHtracks['longitude']):
        SHfirst = SHfirst[0:-1]
    SHtimes = SHtracks.variables['time'][SHfirst]
    SHdateunits = SHtracks.variables['time'].units
    SHdates= num2date(SHtimes, SHdateunits)
        
    #get those seasonal bits        
    NH_DJF, NH_MAM, NH_JJA, NH_SON = seasons(NHdates)
    SH_DJF, SH_MAM, SH_JJA, SH_SON = seasons(SHdates)
    
    #annual data
    peryear.append(len(NHfirst) + len(SHfirst))
    perDJF.append(len(NH_DJF) + len(SH_DJF))
    perMAM.append(len(NH_MAM) + len(SH_MAM))
    perJJA.append(len(NH_JJA) + len(SH_JJA))
    perSON.append(len(NH_SON) + len(SH_SON))

    #NH Data
    peryear_NH.append(len(NHfirst))
    perDJF_NH.append(len(NH_DJF))
    perMAM_NH.append(len(NH_MAM))
    perJJA_NH.append(len(NH_JJA))
    perSON_NH.append(len(NH_SON))

    #SH Data
    peryear_SH.append(len(SHfirst))
    perDJF_SH.append(len(SH_DJF))
    perMAM_SH.append(len(SH_MAM))
    perJJA_SH.append(len(SH_JJA))
    perSON_SH.append(len(SH_SON))
    
    NHtracks.close()
    SHtracks.close()
    
#Might want to save stuff??? This seems to run pretty quickly though

# Shove wavecounts into Pandas Dataframes so I can have global and hemisphere all on one plot
wavecounts = pd.DataFrame({('Global','Year'): peryear,
                          ('Global','DJF'): perDJF,
                          ('Global','MAM'): perMAM,
                          ('Global','JJA'): perJJA,
                          ('Global','SON'): perSON,
                          ('NH','Year'): peryear_NH,
                          ('NH','DJF'): perDJF_NH,
                          ('NH','MAM'): perMAM_NH,
                          ('NH','JJA'): perJJA_NH,
                          ('NH','SON'): perSON_NH,
                          ('SH','Year'): peryear_SH,
                          ('SH','DJF'): perDJF_SH,
                          ('SH','MAM'): perMAM_SH,
                          ('SH','JJA'): perJJA_SH,
                          ('SH','SON'): perSON_SH,
                          })


# # having means and standard deviations is useful
# def charStats(characteristic, name='the characteristic'):
#     mean = np.nanmean(characteristic)
#     stdev = np.nanstd(characteristic)
#     
#     print(f'the mean of {name} is {mean} and the standard deviation is {stdev}')
# 
# 
# 
# charStats(wavecounts['Global','Year'], 'global wavecounts')
#     
# charStats(max_int, 'maximum intensity')
# charStats(max_int_NH, 'maximum intensity, NH,')
# charStats(max_int_SH, 'maximum intensity, SH,')
# 
# charStats(net_dist, 'propagation distance')
# charStats(net_dist_NH, 'propagation distance, NH,')
# charStats(net_dist_SH, 'propagation distance, SH,')
# 
# charStats(avg_spd, 'average speed')
# charStats(avg_spd_NH, 'average speed, NH,')
# charStats(avg_spd_SH, 'average speed, SH,')

# BOXPLOT ALL THE THINGS!!!!
# OR ACTUALLY VIOLINPLOT!!!!

mm = 1/25.4

if level == 700:
    panel = 'a'
elif level == 850:
    panel = 'b'

# #TEW counts
# fig = plt.figure(figsize=(174*mm, 234*mm/2), dpi=600)
# ax = plt.axes()
# sns.boxplot(data=wavecounts, orient='h')
# #ax.set_title(f'Number of TEWs at {level} hPa')
# 
# if level == 700:
#     plt.text(-190, -0.5, f'{panel})', fontweight='bold')
# elif level == 850:
#     plt.text(-125, -0.5, f'{panel})', fontweight='bold')
# 
# ax.set_xlabel('Count')
# plt.savefig(f'./PublicationVersions/distributions/wavecounts.boxplot.{level}.noTC.tiff', bbox_inches='tight', dpi=600)

#three tone palette for the following
palette = ['lightcoral', 'limegreen', 'cornflowerblue']
#and lables
labels = ['Global','NH','SH']

#max intensity
fig2 = plt.figure(figsize=(174*mm, 234*mm/2), dpi=600)
ax2 = plt.axes()
ax2.ticklabel_format(axis='x', style='sci', scilimits=(-5,-5))
sns.violinplot(data=[max_int, max_int_NH, max_int_SH], orient='h', palette=palette, cut=0)
ax2.set_yticklabels(labels)
#ax2.set_title(f'TEW Intensities at {level} hPa')

if level == 700:
    plt.text(-0.000025, -0.4, f'{panel})', fontweight='bold')
elif level == 850:
    plt.text(-0.000025, -0.4, f'{panel})', fontweight='bold')

ax2.set_xlabel('CV Intensity (1/s)')
ax2.set_xlim(0,25e-5)
plt.xticks(np.arange(0, (25e-5+2.5e-5), 2.5e-5)) #still need it in scientific notation       
plt.savefig(f'./PublicationVersions/distributions/intensity.violinplot.{level}.noTC.tiff', bbox_inches='tight', dpi=600)

# #wave distance no TC
# fig3 = plt.figure(figsize=(174*mm, 234*mm/2), dpi=600)
# ax3 = plt.axes()
# sns.violinplot(data=[net_dist, net_dist_NH, net_dist_SH], orient='h', palette=palette)
# ax3.set_yticklabels(labels)
# #ax3.set_title(f'TEW Propagation Distances at {level} hPa')
# 
# if level == 700:
#     plt.text(-2250, -0.4, f'{panel})', fontweight='bold')
# elif level == 850:
#     plt.text(-2500, -0.4, f'{panel})', fontweight='bold')
# ax3.set_xlabel('Distance (km)')
# plt.savefig(f'./PublicationVersions/distributions/wavedist.violinplot.{level}.noTC.tiff', bbox_inches='tight', dpi=600)
# 
# #average speed no TC
# fig4 = plt.figure(figsize=(174*mm, 234*mm/2), dpi=600)
# ax4 = plt.axes()
# sns.violinplot(data=[avg_spd, avg_spd_NH, avg_spd_SH], orient='h', palette=palette)
# ax4.set_yticklabels(labels)
# #ax4.set_title(f'TEW Average Speed at {level} hPa')
# 
# if level == 700:
#     plt.text(-3.5, -0.4, f'{panel})', fontweight='bold')
# elif level == 850:
#     plt.text(-3, -0.4, f'{panel})', fontweight='bold')
#     
# ax4.set_xlabel('Speed (m/s)')
# plt.savefig(f'./PublicationVersions/distributions/avgspeed.violinplot.{level}.noTC.tiff', bbox_inches='tight', dpi=600)
