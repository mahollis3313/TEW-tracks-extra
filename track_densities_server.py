#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:09:03 2020

@author: margaret

track_densities_server.py
For plotting not only genesis densities, but also the track and lysis densities
THIS VERSION READS IN MULTIPLE FILES, to create density plots for MULTIPLE years

19 Oct 2020
- added some comments, as well as this documentation at the top
- standardized the colorbars
- remvoed color for lowest level

Early Jan 2021 I think
- added max intensity to things plotted

September 2022
- updated file paths for updated track files
- changed the if for last track being zero length to a while for a case where
  a file ended with two zero length tracks
  
March 2023: The Publication Updates
- fixed overlap between colorbar and lat labels
- all the following changes were made to adhere to Clim. Dyn. publishing requirements
- changed fig sizes
- removed fig titles
- added a/b part labels
- changed output path to PublicationVersions folder
- increased resolution to 600 dpi
- changed output file type to tiff

Jan 2025: I just need a map of the counts of TEWs in each lonbin so quickly doing that
No major changes to how things are displayed or what not just outputting a png again
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

# netCDF stuff
import netCDF4
from netCDF4 import Dataset

# # # # # # # # # # #
#variables from input
# # # # # # # # # # #
level = int(sys.argv[1])
binlon = float(sys.argv[2])
binlat = float(sys.argv[3])

# # # # # # # # #
#read in the data
# # # # # # # # # 
trackfiles = sorted(glob.glob(f'/data2/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.*.nc'))
#year and hemisphere is combined  in wildcard

#just gonna loop through trackfile a few times to simulate importing multiple years of tracks
#prolly wanna set up the arrays for these things to go into first...
FirstPt = []
NumPts_Gen = []
lats = []
lons = []
first_lat = []
first_lon = []
last_lat = []
last_lon = []
max_lats = []
max_lons = []


for trackfile in trackfiles:
    tracks = Dataset(trackfile, mode='r')
    FirstPt_offset = len(FirstPt)
        
    #stuff about where the tracks start
    FirstPt.extend(tracks['FIRST_PT'][:].data+FirstPt_offset)
    NumPts_Gen.extend(tracks['NUM_PTS'][:])
    
    thisfile = tracks['FIRST_PT'][:]
    
    #the tracks themselves
    lats.extend(tracks['latitude'][:].data)
    lons.extend(tracks['longitude'][:].data)
    
    while thisfile[-1] == len(tracks['longitude']):
        thisfile = thisfile[0:-1]
        
    #genesis points
    first_lat.extend(tracks['latitude'][thisfile])
    first_lon.extend(tracks['longitude'][thisfile])
    
    #lysis points
    last_lat.extend(tracks['latitude'][thisfile-1])
    last_lon.extend(tracks['longitude'][thisfile-1])
    
    # #max intensity lat/lon
    # max_lats.extend(tracks['max_lat'][:])
    # max_lons.extend(tracks['max_lon'][:])

    tracks.close()


#make those plots

#climate dynamics specifies maximum figure sizes in mm, so need a conversion factor
mm = 1/25.4

#a/b panel labels
if level == 700:
    panel = 'a'
elif level == 850:
    panel = 'b'

# =============================================================================
# #genesis
# fig1 = plt.figure(figsize=(174*mm,234*mm),dpi=600)
# ax1 = plt.axes(projection=ccrs.PlateCarree())
# ax1.coastlines(resolution='110m')
# 
# sns.histplot(x=first_lon, y=first_lat, stat='count', ax=ax1, binwidth=(binlon, binlat), cbar=True, cbar_kws=dict(shrink=0.128, pad=0.08))
# #ax1.set_title(f'Genesis Point Counts at {level} hPa')
# ax1.set_ylim(-40, 40)
# 
# grid1 = ax1.gridlines(draw_labels=True)
# grid1.top_labels = False
# grid1.left_labels = False
# grid1.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
# grid1.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])
# 
# plt.text(-200, 30, f'{panel})')
# 
# #plt.show()
# plt.savefig(f'./PublicationVersions/genesis_density.{level}.v2.noTC.tiff', bbox_inches = "tight", dpi = 600)
# =============================================================================

# =============================================================================
# #lysis
# fig2 = plt.figure(figsize=(174*mm,234*mm),dpi=600)
# ax2 = plt.axes(projection=ccrs.PlateCarree())
# ax2.coastlines(resolution='110m')
# 
# sns.histplot(x=last_lon, y=last_lat, stat='count', ax=ax2, binwidth=(binlon, binlat), cbar=True, cbar_kws=dict(shrink=0.241, pad=0.08))
# #ax2.set_title(f'End Point Counts at {level} hPa')
# ax2.set_ylim(-75,75)
# 
# grid2 = ax2.gridlines(draw_labels=True)
# grid2.top_labels = False
# grid2.left_labels = False
# grid2.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
# grid2.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])
# 
# plt.text(-200, 65, f'{panel})')
# 
# #plt.show()
# plt.savefig(f'./PublicationVersions/lysis_density.{level}.v2.noTC.tiff', bbox_inches = "tight", dpi = 600)
# =============================================================================

#track densities
fig3 = plt.figure(figsize=(174*mm,234*mm),dpi=600)
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.coastlines(resolution='110m')

sns.histplot(x=lons, y=lats, stat='count', ax=ax3, binwidth=(binlon, binlat),
             cbar=True, cbar_kws=dict(shrink=0.241, pad=0.08))
#ax3.set_title(f'Track Point Counts at {level} hPa')
ax3.set_ylim(-75, 75)

grid3 = ax3.gridlines(draw_labels=True)
grid3.top_labels = False
grid3.left_labels = False
grid3.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid3.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.text(-200, 65, f'{panel})')

#plt.show()
plt.savefig(f'./Lonbins/track_density.{level}.v2.noTC.png', bbox_inches = "tight", dpi = 600)

# #max intensity locations
# fig4 = plt.figure(dpi=300)
# ax4 = plt.axes(projection=ccrs.PlateCarree())
# ax4.coastlines(resolution='110m')
#
# sns.histplot(x=max_lons, y=max_lats, stat='count', ax=ax4, binwidth=(binlon, binlat), cbar=True)
# ax4.set_title(f'Maximum Intensities {level}hPa')
# #plt.show()
# plt.savefig(f'max_intensity_map.{level}.v2.noTC.png', bbox_inches = "tight", dpi = 300)
