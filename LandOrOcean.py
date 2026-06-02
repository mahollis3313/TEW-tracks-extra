# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:21:02 2024

@author: Margaret

LandOrOcean.py

pretty much what it sounds like, this code will determine if a TEW point is
over land or ocean (or in a coastal grid box)

for our purposes, we will use FROCEAN (fraction of ocean)
Ocean   = FROCEAN > 0.5
Land    = FROCEAN < 0.5

Initially considered having a "coastal" range of 0.25-0.75, which would correspond with
suggested land and ocean thresholds from the NASA IMERG data page,
but it didn't have many points in it relative to the number of TEWs

it will then save the result to a new variable in the TEW track file

UPDATED 14 Apr 2025 to use IMERG mask instead of MERRA2 mask
"""

import numpy as np
import sys
import glob
# netCDF stuff
import netCDF4
from netCDF4 import Dataset
#plotting stuff
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# functions I "wrote" (found on the internet years ago)
# for figuring out the closest grid point to the track point

# great circle distance
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1 = map(np.radians, [lon1, lat1])
    lon2, lat2 = map(np.radians, [lon2, lat2])
    return 6371 * (np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)))

# and a function for getting indices for TEW point and center of the mask
# shamelessly stolen from https://stackoverflow.com/questions/30180241/numpy-get-the-column-and-row-index-of-the-minimum-value-of-a-2d-array
# but modified to give back indices as ints
def find_min_idx(x):
    k = x.argmin()
    ncol = x.shape[1]
    return int(round(k/ncol)), k%ncol

level = int(sys.argv[1])

#the ocean mask data
#oceanfile = '/mnt/c/Users/Margaret/Research/Data/MERRA2_101.const_2d_asm_Nx.00000000.nc4'
oceanfile = '/data/deluge/scratch/IMERG/GPM_IMERG_LandSeaMask.2.nc4'
oceandata= Dataset(oceanfile, mode='r')

frocean = oceandata['FROCEAN'][0,:,:]
oceanlon = oceandata['lon'][:]
oceanlat = oceandata['lat'][:]

oceanlons, oceanlats = np.meshgrid(oceanlon, oceanlat)
oceandata.close()

filelist = []
for year in range(1981,2023):
    # filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.NH.NoTC.nc')
    # filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.SH.NoTC.nc')
    
    filelist.append(f'/mnt/c/Users/Margaret/Research/Data/FixedTracks/{level}hPa/post_filter.{year}.NH.NoTC.nc')
    filelist.append(f'/mnt/c/Users/Margaret/Research/Data/FixedTracks/{level}hPa/post_filter.{year}.SH.NoTC.nc')

for trackfile in filelist:
    
    tracks = Dataset(trackfile, mode='r+')

    TrackID = tracks['TRACK_ID'][:]
    FirstPt = tracks['FIRST_PT'][:]
    NumPts  = tracks['NUM_PTS'][:]
    

    hemi  = trackfile[-10:-8]
    
    lats = tracks['latitude'][:]
    lons = tracks['longitude'][:]
    
    #do something to assign land/ocean
    
    landocean = []
    
    for tracklat, tracklon in zip (lats, lons):
        
        #find the nearest MERRA point to the TEW point
        great_circle_dist = great_circle(tracklon, tracklat, oceanlons, oceanlats)
        latindex, lonindex = find_min_idx(great_circle_dist)
        
        
        
        #assign land/ocean
        if frocean[latindex, lonindex] >= 0.5:
            landocean.append(1)
        else:
            landocean.append(0)
        
    
    #save the land/ocean list
    out_landocean = tracks.createVariable('LandOcean', 'i8', ('record'))
    
    #give the new variables some attributes
    out_landocean.standard_name = 'Land Ocean'
    out_landocean.long_name     = 'Land or Ocean'
    out_landocean.description   = 'Boolean for Land (0) or Ocean (1)'
    
    
    #now actually write the data
    tracks['LandOcean'][:] = landocean
    
    tracks.close()


# =============================================================================
# #what does our ocean/coastal/land determination even look like?
# plt.figure(figsize=(16,10), dpi=300)
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines(resolution='110m')
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# plt.contourf(oceanlon, oceanlat, frocean, levels=[0.49,0.51], colors=['brown', 'g', 'blue'], extend='both')
# plt.savefig('LandOrOcean.png')
# =============================================================================

# =============================================================================
# #are land and ocean correctly assigned?
# plt.figure(figsize=(16,10), dpi=300)
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines(resolution='110m')
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# plt.scatter(lons, lats, c=landocean)
# plt.savefig('LandOrOcean.png')
# =============================================================================
