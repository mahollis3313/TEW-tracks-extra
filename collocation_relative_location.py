#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collocation_relative_location.py
Created on Tues Aug 29 16:14 2023

@author: margaret

code for plotting the frequency of locations of 700 hPa collocation centers
relative to their respective 850 hPa center

"""

# # # # # # # # # #
# import the things
# # # # # # # # # #
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import glob
import sys
#import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import seaborn as sns
from netCDF4 import Dataset, num2date


# # # # # # # # #
#read in the data
# # # # # # # # # 
#year and hemisphere is combined  in wildcard
#trackfiles = sorted(glob.glob(f'/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredNoTCUpdated/{level}hPa/post_filter.*.nc'))
#trackfiles = sorted(glob.glob(f'/media/margaret/Big Volume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredTracks/{level}hPa/post_filter.*.nc'))
#trackfiles = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.nc'))

#if need certain years:
NH850trackfiles = []
SH850trackfiles = []
NH700trackfiles = []
SH700trackfiles = []
for year in range(1981, 2022):
    NH850trackfiles.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/850hPa/post_filter.{year}.NH.NoTC.nc')
    SH850trackfiles.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/850hPa/post_filter.{year}.SH.NoTC.nc')
    NH700trackfiles.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/700hPa/post_filter.{year}.NH.NoTC.nc')
    SH700trackfiles.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/700hPa/post_filter.{year}.SH.NoTC.nc')
    
    

def gettracks(trackfiles, level):
    
    #arrays for things to go into
    matchlats = []
    matchlons = []
    
    if level == 850:
        lev=0
    elif level == 700:
        lev=1
    
    #looop through the files and put things into the arrays
    for trackfile in trackfiles:
        tracks = Dataset(trackfile, mode='r')
        
        firstdate = num2date(tracks['time'][0], units=tracks['time'].units)
        year = firstdate.year
        hemishort = trackfile[-10:-8]
        
        # the file for the collocated points
        matchindices = np.loadtxt(f'/home/mhollis/track_stats/CollocationUpdated/match_indices_{year}_{hemishort}.csv', delimiter=',').astype('int')
        
        #the tracks themselves
        matchlats.extend(tracks['latitude'][matchindices[:,lev]].data)
        matchlons.extend(tracks['longitude'][matchindices[:,lev]].data)
    
        tracks.close()
        
    return np.asarray(matchlats), np.asarray(matchlons)
    
matchlatsNH850, matchlonsNH850 = gettracks(NH850trackfiles, 850)
matchlatsSH850, matchlonsSH850 = gettracks(SH850trackfiles, 850)
matchlatsNH700, matchlonsNH700 = gettracks(NH700trackfiles, 700)
matchlatsSH700, matchlonsSH700 = gettracks(SH700trackfiles, 700)

NHdx = matchlonsNH700 - matchlonsNH850
NHdy = matchlatsNH700 - matchlatsNH850
SHdx = matchlonsSH700 - matchlonsSH850
SHdy = matchlatsSH700 - matchlatsSH850

#can I quickly correct for the date line in the dx vars?
def datelinecorrection(dx):
    


def plotDartboard(dx, dy, hemi):    
    fig = plt.figure(figsize=(7,5),dpi=300)
    ax = plt.axes()
    sns.histplot(x=dx, y=dy, stat='count', ax=ax, binwidth=(0.625, 0.5), cbar=True, cbar_kws=dict(shrink=0.75))
    #plt.scatter(dx, dy)
    #plt.colorbar()
    ax.set_xlim(-6,6)
    ax.set_ylim(-5,5)
    ax.set_xlabel('Distance from TEW Center (degrees)')
    ax.set_ylabel('Distance from TEW Center (degrees)')
    #ax.set_title(f'Composite TEW precipitation rate, mm/hr \n for {comptype} TEW points in the {hemitext(hemi)} hemisphere at {level} hPa')
    ax.text(-7, 4.5, f'{hemi}')
    
    ax.grid()
    
    plt.savefig(f'/home/mhollis/track_stats/collocation_dartboard_{hemi}.png', bbox_inches = "tight", dpi = 300)
    
plotDartboard(NHdx, NHdy, 'NH')
plotDartboard(SHdx, SHdy, 'SH')