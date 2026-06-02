#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 11 Nov 2025
@author: margaret

things_that_dont_need_map_lonbins.py

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

16 Oct 2025:
so I actually pulled most of this code from ~/general_exam/plot_dynamic_comparisons.py
because it does regions
just also ITCZ

11 Nov 2025:
another overhaul this time to get the stats by lonbin/landocean and save it out to a csv
this time the overhaul largely comes form ~/moisture_modes/plot_scatter_crit1_cumulative.py
because I need the speeds for moisture modes and n-mode
and I may as well do the other vars this originally did at the same time

12 Feb 2026:
yeah I want hemispheres now too

23 Apr 2026:
just added a quick TRACK_ID var to get track counts for updating my tables
"""

#import packages
import sys
import numpy as np
from scipy import stats
import pandas as pd
import glob
import netCDF4
from netCDF4 import Dataset, num2date

#paths to files
trackfiles850NH = []
trackfiles850SH = []
trackfiles700NH = []
trackfiles700SH = []
for year in range(2001, 2024):
#for year in range(1981, 2024):
    trackfiles850NH.append(f'/data2/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/850hPa/post_filter.{year}.NH.NoTC.nc')
    trackfiles850SH.append(f'/data2/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/850hPa/post_filter.{year}.SH.NoTC.nc')
    trackfiles700NH.append(f'/data2/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/700hPa/post_filter.{year}.NH.NoTC.nc')
    trackfiles700SH.append(f'/data2/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/700hPa/post_filter.{year}.SH.NoTC.nc')
    
def getdata(filelist):
    
    longitude = []
    landocean = []
    speeds    = []
    avgspeeds = []
    intens    = []
    maxints   = []
    dists     = []
    netdists  = []
    lifetimes = []
    #localhrs  = []
    TRACK_IDs = []
    
    for file in filelist:
        
        tracks = Dataset(file, mode='r')

        lons = tracks['longitude'][:]    
        
        landbool = tracks['LandOcean'][:]
        
        speed  = tracks['speed'][:]
        avgspd = tracks['avg_spd'][:]
        inten  = tracks['curvature_vorticity'][:]
        maxint = tracks['max_int'][:]
        dist   = tracks['all_dist_2d'][:]
        netdist = tracks['net_dist_2d'][:]
        lifetime = tracks['NUM_PTS'][:]
        lifetime = lifetime/8
        TRACK_ID = tracks['TRACK_ID'][:]
        
        #local = tracks['local_solar'][:]
        
        longitude.extend(lons)
        landocean.extend(landbool)
        speeds.extend(speed)
        avgspeeds.extend(avgspd)
        intens.extend(inten)
        maxints.extend(maxint)
        dists.extend(dist)
        netdists.extend(netdist)
        lifetimes.extend(lifetime)
        #localhrs.extend(local)   
        TRACK_IDs.extend(TRACK_ID)   
    
    #for now only returning the things that are the same length as longitude
    #we can work out adding the dataframes and regions for the other stuff later if needed
    #return np.asarray(longitude), np.asarray(landocean), np.asarray(avgspeeds), np.asarray(maxints), np.asarray(netdists), np.asarray(lifetimes)
    return np.asarray(longitude), np.asarray(landocean), np.asarray(speeds), np.asarray(intens), np.asarray(TRACK_IDs)

lon850NH, landocean850NH, speed850NH, int850NH, IDs850NH = getdata(trackfiles850NH)
lon850SH, landocean850SH, speed850SH, int850SH, IDs850SH = getdata(trackfiles850SH)
lon700NH, landocean700NH, speed700NH, int700NH, IDs700NH = getdata(trackfiles700NH)
lon700SH, landocean700SH, speed700SH, int700SH, IDs700SH = getdata(trackfiles700SH)

#get land/ocean indices to make land/ocean dataframes
#this would also be what to do to separate out other things like ITCZ
land850NH  = np.where(landocean850NH==0)[0]
ocean850NH = np.where(landocean850NH==1)[0]
land850SH  = np.where(landocean850SH==0)[0]
ocean850SH = np.where(landocean850SH==1)[0]

land700NH  = np.where(landocean700NH==0)[0]
ocean700NH = np.where(landocean700NH==1)[0]
land700SH  = np.where(landocean700SH==0)[0]
ocean700SH = np.where(landocean700SH==1)[0]

#and then set up the actual dataframes
waves850NH = pd.DataFrame(np.asarray([lon850NH, speed850NH, int850NH]).T, columns=['longitude', 'speed', 'maxint'])
waves850SH = pd.DataFrame(np.asarray([lon850SH, speed850SH, int850SH]).T, columns=['longitude', 'speed', 'maxint'])
waves700NH = pd.DataFrame(np.asarray([lon700NH, speed700NH, int700NH]).T, columns=['longitude', 'speed', 'maxint'])
waves700SH = pd.DataFrame(np.asarray([lon700SH, speed700SH, int700SH]).T, columns=['longitude', 'speed', 'maxint'])

#then some smaller dataframes of just the land/ocean
landwaves850NH  = pd.DataFrame(np.asarray([lon850NH[land850NH], speed850NH[land850NH], int850NH[land850NH]]).T, columns=['longitude', 'speed', 'maxint'])
oceanwaves850NH = pd.DataFrame(np.asarray([lon850NH[ocean850NH], speed850NH[ocean850NH], int850NH[ocean850NH]]).T, columns=['longitude', 'speed', 'maxint'])
landwaves700NH  = pd.DataFrame(np.asarray([lon700NH[land700NH], speed700NH[land700NH], int700NH[land700NH]]).T, columns=['longitude', 'speed', 'maxint'])
oceanwaves700NH = pd.DataFrame(np.asarray([lon700NH[ocean700NH], speed700NH[ocean700NH], int700NH[ocean700NH]]).T, columns=['longitude', 'speed', 'maxint'])

landwaves850SH  = pd.DataFrame(np.asarray([lon850SH[land850SH], speed850SH[land850SH], int850SH[land850SH]]).T, columns=['longitude', 'speed', 'maxint'])
oceanwaves850SH = pd.DataFrame(np.asarray([lon850SH[ocean850SH], speed850SH[ocean850SH], int850SH[ocean850SH]]).T, columns=['longitude', 'speed', 'maxint'])
landwaves700SH  = pd.DataFrame(np.asarray([lon700SH[land700SH], speed700SH[land700SH], int700SH[land700SH]]).T, columns=['longitude', 'speed', 'maxint'])
oceanwaves700SH = pd.DataFrame(np.asarray([lon700SH[ocean700SH], speed700SH[ocean700SH], int700SH[ocean700SH]]).T, columns=['longitude', 'speed', 'maxint'])

#add a lonbin index column
waves850NH['lonbin'] = waves850NH['longitude'] // 15
waves700NH['lonbin'] = waves700NH['longitude'] // 15

waves850SH['lonbin'] = waves850SH['longitude'] // 15
waves700SH['lonbin'] = waves700SH['longitude'] // 15

landwaves850NH['lonbin']  = landwaves850NH['longitude'] // 15
oceanwaves850NH['lonbin'] = oceanwaves850NH['longitude'] // 15
landwaves700NH['lonbin']  = landwaves700NH['longitude'] // 15
oceanwaves700NH['lonbin'] = oceanwaves700NH['longitude'] // 15

landwaves850SH['lonbin']  = landwaves850SH['longitude'] // 15
oceanwaves850SH['lonbin'] = oceanwaves850SH['longitude'] // 15
landwaves700SH['lonbin']  = landwaves700SH['longitude'] // 15
oceanwaves700SH['lonbin'] = oceanwaves700SH['longitude'] // 15

#group by lonbin
lonbin700NH = waves700NH.groupby(['lonbin']).indices
lonbin850NH = waves850NH.groupby(['lonbin']).indices

lonbin700SH = waves700SH.groupby(['lonbin']).indices
lonbin850SH = waves850SH.groupby(['lonbin']).indices

landlonbin850NH = landwaves850NH.groupby(['lonbin']).indices
oceanlonbin850NH = oceanwaves850NH.groupby(['lonbin']).indices
landlonbin700NH = landwaves700NH.groupby(['lonbin']).indices
oceanlonbin700NH = oceanwaves700NH.groupby(['lonbin']).indices

landlonbin850SH = landwaves850SH.groupby(['lonbin']).indices
oceanlonbin850SH = oceanwaves850SH.groupby(['lonbin']).indices
landlonbin700SH = landwaves700SH.groupby(['lonbin']).indices
oceanlonbin700SH = oceanwaves700SH.groupby(['lonbin']).indices

#we don't particularly need to make plots right now
#so for now we're going to focus on just the statistical output
#but go back to the moisture mode criteria heatmaps if plots are ever needed
#it does need to be a function though because we need to do this for all of the variables

def keydescribe(varname, hemi):
    if hemi == 'NH':
        waves700 = waves700NH
        waves850 = waves850NH
        
        landwaves700 = landwaves700NH
        landwaves850 = landwaves850NH
        
        oceanwaves700 = oceanwaves700NH
        oceanwaves850 = oceanwaves850NH
        
        lonbin700 = lonbin700NH
        lonbin850 = lonbin850NH

        landlonbin700 = landlonbin700NH        
        landlonbin850 = landlonbin850NH
        
        oceanlonbin700 = oceanlonbin700NH
        oceanlonbin850 = oceanlonbin850NH
        
    elif hemi=='SH':
        waves700 = waves700SH
        waves850 = waves850SH
        
        landwaves700 = landwaves700SH
        landwaves850 = landwaves850SH
        
        oceanwaves700 = oceanwaves700SH
        oceanwaves850 = oceanwaves850SH
        
        lonbin700 = lonbin700SH
        lonbin850 = lonbin850SH

        landlonbin700 = landlonbin700SH        
        landlonbin850 = landlonbin850SH
        
        oceanlonbin700 = oceanlonbin700SH
        oceanlonbin850 = oceanlonbin850SH

    stats700 = np.zeros((25,6))
    stats850 = np.zeros((25,6))
    
    landstats700 = np.zeros((25,6))
    oceanstats700 = np.zeros((25,6))
    
    landstats850 = np.zeros((25,6))
    oceanstats850 = np.zeros((25,6))
    
    #global stats
    globalstats700 = stats.describe(waves700[f'{varname}'], nan_policy='omit')
    globalstats850 = stats.describe(waves850[f'{varname}'], nan_policy='omit')
    
    globallandstats700 = stats.describe(landwaves700[f'{varname}'], nan_policy='omit')
    globallandstats850 = stats.describe(landwaves850[f'{varname}'], nan_policy='omit')
    
    globaloceanstats700 = stats.describe(oceanwaves700[f'{varname}'], nan_policy='omit')
    globaloceanstats850 = stats.describe(oceanwaves850[f'{varname}'], nan_policy='omit')
    
    #and go in the last row of the arrays
    stats700[-1,:] = [-1, globalstats700.nobs, globalstats700.minmax[0], globalstats700.minmax[0], globalstats700.mean, globalstats700.variance]
    stats850[-1,:] = [-1, globalstats850.nobs, globalstats850.minmax[0], globalstats850.minmax[0], globalstats850.mean, globalstats850.variance]
    
    landstats700[-1,:] = [-1, globallandstats700.nobs, globallandstats700.minmax[0], globallandstats700.minmax[1], globallandstats700.mean, globallandstats700.variance]
    landstats850[-1,:] = [-1, globallandstats850.nobs, globallandstats850.minmax[0], globallandstats850.minmax[1], globallandstats850.mean, globallandstats850.variance]
    
    oceanstats700[-1,:] = [-1, globaloceanstats700.nobs, globaloceanstats700.minmax[0], globaloceanstats700.minmax[1], globaloceanstats700.mean, globaloceanstats700.variance]
    oceanstats850[-1,:] = [-1, globaloceanstats850.nobs, globaloceanstats850.minmax[0], globaloceanstats850.minmax[1], globaloceanstats850.mean, globaloceanstats850.variance]
    
    for lonbin in range(24):
        
        #all-terrain stats exist in each lonbin
        stat700 = stats.describe(waves700.iloc[lonbin700[lonbin]][f'{varname}'], nan_policy='omit')
        stat850 = stats.describe(waves850.iloc[lonbin850[lonbin]][f'{varname}'], nan_policy='omit')
        
        stats700[lonbin,:] = [lonbin*15, stat700.nobs, stat700.minmax[0], stat700.minmax[1], stat700.mean, stat700.variance]
        stats850[lonbin,:] = [lonbin*15, stat850.nobs, stat850.minmax[0], stat850.minmax[1], stat850.mean, stat850.variance]
        
        #but specific terrains may not
        if lonbin in landlonbin700:
            lonlandstat700 = stats.describe(landwaves700.iloc[landlonbin700[lonbin]][f'{varname}'], nan_policy='omit')
            landstats700[lonbin,:] = [lonbin*15, lonlandstat700.nobs, lonlandstat700.minmax[0], lonlandstat700.minmax[1], lonlandstat700.mean, lonlandstat700.variance]

        else:
            print('empty group -- continuing')
        
        if lonbin in oceanlonbin700:
            lonoceanstat700 = stats.describe(oceanwaves700.iloc[oceanlonbin700[lonbin]][f'{varname}'], nan_policy='omit')
            oceanstats700[lonbin,:] = [lonbin*15, lonoceanstat700.nobs, lonoceanstat700.minmax[0], lonoceanstat700.minmax[1], lonoceanstat700.mean, lonoceanstat700.variance]
        
        else:
            print('empty group -- continuing')
            
        if lonbin in landlonbin850:
            lonlandstat850 = stats.describe(landwaves850.iloc[landlonbin850[lonbin]][f'{varname}'], nan_policy='omit')
            landstats850[lonbin,:] = [lonbin*15, lonlandstat850.nobs, lonlandstat850.minmax[0], lonlandstat850.minmax[1], lonlandstat850.mean, lonlandstat850.variance]

        else:
            print('empty group -- continuing')
        
        if lonbin in oceanlonbin850:
            lonoceanstat850 = stats.describe(oceanwaves850.iloc[oceanlonbin850[lonbin]][f'{varname}'], nan_policy='omit')
            oceanstats850[lonbin,:] = [lonbin*15, lonoceanstat850.nobs, lonoceanstat850.minmax[0], lonoceanstat850.minmax[1], lonoceanstat850.mean, lonoceanstat850.variance]
        
        else:
            print('empty group -- continuing')
            
    return stats700, stats850, landstats700, landstats850, oceanstats700, oceanstats850

avgspeed700NH, avgspeed850NH, landspeed700NH, landspeed850NH, oceanspeed700NH, oceanspeed850NH = keydescribe('speed', 'NH')
avgspeed700SH, avgspeed850SH, landspeed700SH, landspeed850SH, oceanspeed700SH, oceanspeed850SH = keydescribe('speed', 'SH')
intensity700NH, intensity850NH, landint700NH, landint850NH, oceanint700NH, oceanint850NH       = keydescribe('maxint', 'NH') 
intensity700SH, intensity850SH, landint700SH, landint850SH, oceanint700SH, oceanint850SH       = keydescribe('maxint', 'SH') 
#distance700, distance850, landdist700, landdist850, oceandist700, oceandist850     = keydescribe('netdist')
#lifetime700, lifetime850, landlife700, landlife850, oceanlife700, oceanlife850     = keydescribe('lifetime')

#now write these out to csv files with headers
#speed
np.savetxt('./descriptive_stats/avg_speed_700NH.csv', avgspeed700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850NH.csv', avgspeed850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/avg_speed_700SH.csv', avgspeed700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850SH.csv', avgspeed850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/avg_speed_700NH_land.csv', landspeed700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850NH_land.csv', landspeed850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/avg_speed_700SH_land.csv', landspeed700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850SH_land.csv', landspeed850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/avg_speed_700NH_ocean.csv', oceanspeed700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850NH_ocean.csv', oceanspeed850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/avg_speed_700SH_ocean.csv', oceanspeed700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/avg_speed_850SH_ocean.csv', oceanspeed850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

#intensity
np.savetxt('./descriptive_stats/intensity_700NH.csv', intensity700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850NH.csv', intensity850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/intensity_700SH.csv', intensity700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850SH.csv', intensity850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/intensity_700NH_land.csv', landint700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850NH_land.csv', landint850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/intensity_700SH_land.csv', landint700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850SH_land.csv', landint850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/intensity_700NH_ocean.csv', oceanint700NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850NH_ocean.csv', oceanint850NH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

np.savetxt('./descriptive_stats/intensity_700SH_ocean.csv', oceanint700SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
np.savetxt('./descriptive_stats/intensity_850SH_ocean.csv', oceanint850SH, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')

# #distance
# np.savetxt('./descriptive_stats/distance_700.csv', distance700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/distance_850.csv', distance850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
# np.savetxt('./descriptive_stats/distance_700_land.csv', landdist700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/distance_850_land.csv', landdist850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
# np.savetxt('./descriptive_stats/distance_700_ocean.csv', oceandist700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/distance_850_ocean.csv', oceandist850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
# #lifetime
# np.savetxt('./descriptive_stats/lifetime_700.csv', lifetime700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/lifetime_850.csv', lifetime850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
# np.savetxt('./descriptive_stats/lifetime_700_land.csv', landlife700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/lifetime_850_land.csv', landlife850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
# np.savetxt('./descriptive_stats/lifetime_700_ocean.csv', oceanlife700, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# np.savetxt('./descriptive_stats/lifetime_850_ocean.csv', oceanlife850, delimiter=',', header='lonbin, nobs, min, max, mean, variance', comments='')
# 
#     