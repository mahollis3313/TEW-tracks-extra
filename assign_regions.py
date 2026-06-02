#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 2023

@author: margaret

Code to assign regions to TEWs as a variable
Region list and sub-basins

Northern Hemisphere
- North Atlantic (45E - ~80W)
    - MDR (45E - 75W)
    - western Caribbean (75W - 80W, south of 15N)
    - GOM and eastern US (75W - ~85W, north of 15N)
- North East Pacific
    - equatorial (75W - 180, south of 15N)
    - subtropical (~85W - 180, north of 15N)
- North West Pacific (180 - 110E)
- Indian Ocean (110E - 45E)

Southern Hemisphere
- South Atlantic (45E - 75W)
- South East Pacific (75W - 180)
- South West Pacific (180 - 155E)
- Indian Ocean (155E - 45E)
"""

import numpy as np
import sys
import glob
from netCDF4 import Dataset

level = int(sys.argv[1])


def NHregion(genesislats, genesislons, NumPts):
    region = []
    for lat, lon in zip(genesislats, genesislons):
        if 45 < lon <= 110: #northern indian ocean basin
            region.append(2)
        elif 110 < lon <= (360-100): #most of the north pacific
            region.append(3)
        elif ((360-100) < lon <= (360-90)) and (0 <= lat < 20): #SE Mex and south
            region.append(3)
        elif ((360-85) < lon <= (360-85)) and (0 <= lat < 15): #most of central America
            region.append(3)
        elif ((360- 70) < lon <= (360-85)) and (0 <=lat < 10): #Panama and south
            region.append(3)
        else: #everything else is the north atlantic
            region.append(1)
    
    # region number 9 for 0-length tracks, append any that need appending at the end
    nolengths = np.where(NumPts == 0)[0]
    for i in nolengths:
        if i >= len(region):
            region.append(9)
        else:
            region[i] = 9
    
    return region
            
def SHregion(genesislats, genesislons, NumPts):
    region = []
    for lat, lon in zip(genesislats, genesislons):
        if 30 < lon <= 160: #southern indian ocean basin
            region.append(5)
        elif 160 < lon <= (360-70): #south pacific
            region.append(6)
        else: #leaving south atlantic
            region.append(4)
            
    # region number 9 for 0-length tracks, append any that need appending at the end
    nolengths = np.where(NumPts == 0)[0]
    for i in nolengths:
        if i >= len(region):
            region.append(9)
        else:
            region[i] = 9
    
    return region




# filelist = sorted(glob.glob(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.1981.NH.NoTC.nc'))

# filelist = []
# for year in range(1981,2023):
#     filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.NH.NoTC.nc')
#     filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.SH.NoTC.nc')

filelist = []
for year in range(2021,2023):
    filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.NH.NoTC.nc')
    filelist.append(f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.SH.NoTC.nc')

for trackfile in filelist:
    
    tracks = Dataset(trackfile, mode='r+')

    TrackID = tracks['TRACK_ID'][:]
    FirstPt = tracks['FIRST_PT'][:]
    NumPts  = tracks['NUM_PTS'][:]
    
    #get rid of the 0-length tracks at the end of files
    while FirstPt[-1] == len(tracks['longitude']):
        FirstPt = FirstPt[0:-1]

    hemi  = trackfile[-10:-8]
    
    genesislats = tracks['latitude'][FirstPt]
    genesislons = tracks['longitude'][FirstPt]
    
    if hemi == 'NH':
        region = NHregion(genesislats, genesislons, NumPts)
    elif hemi == 'SH':
        region = SHregion(genesislats, genesislons, NumPts)
    else:
        print(f'there has been a problem with detecting the hemisphere in {trackfile}')
        break
        
    
    #save the genesis region list
    #delete any lingering variables/dimensions
    #del tracks.variables['GenesisRegion']
    #del tracks.dimensions['string']
    
    #tracks.createDimension('text', 4)
    out_region = tracks.createVariable('genesis_region', 'i8', ('tracks'))
    
    #give the new variables some attributes
    out_region.standard_name = 'Genesis Region'
    out_region.long_name     = 'TEW genesis region'
    
    
    #now actually write the data
    tracks['genesis_region'][:] = region
    
    tracks.close()