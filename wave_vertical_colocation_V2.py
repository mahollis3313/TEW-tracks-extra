#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 13:49:04 2021

@author: margaret

wave_vertical_colocation.py

Comapres wave tracks at 850 and 700 hPa

Something about if this creates new files or modifies the existing ones

"""
# packages to import
import sys
import numpy as np
import numpy.ma as ma
from math import radians, degrees, sin, cos, asin, acos, sqrt

# netCDF stuff
from netCDF4 import Dataset, num2date

# other functions
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1 = map(np.radians, [lon1, lat1])
    lon2, lat2 = map(np.radians, [lon2, lat2])
    return 6371 * (np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)))

#a function that finds minimums in arrays
#shamelessly stolen from https://stackoverflow.com/questions/30180241/numpy-get-the-column-and-row-index-of-the-minimum-value-of-a-2d-array
#but modified to give back indices as ints
def find_min_idx(x):
    k = x.argmin()
    ncol = x.shape[1]
    return int(round(k/ncol)), k%ncol

# # # # # # # # # # # # # #
# User-supplied information
# # # # # # # # # # # # # #

year = int(sys.argv[1])
hemi = int(sys.argv[2])

# # # # # # # # # # # # # # # # # # # #
# opening files and accessing variables
# # # # # # # # # # # # # # # # # # # #

if hemi == 1:
    hemitext = 'northern'
    hemishort = 'NH'
elif hemi == 2:
    hemitext = 'southern'
    hemishort = 'SH'
else:
    print('Hemisphere is not valid.')
    sys.exit()

# The 850 hPa tracks
trackfile_850 = f'/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredNoTCUpdated/850hPa/post_filter.{year}.{hemishort}.NoTC.nc' #f'/data3/TRACK/TRACKOutput/MERRA2/FilteredTracks/850hPa/post_filter.{year}.{hemishort}.nc'
alltracks_850 = Dataset(trackfile_850, mode='r')

# the variables we need, and time converted to datetime
TrackID_850 = alltracks_850['TRACK_ID'][:]
numrecs_850 = alltracks_850.dimensions['record'].size
firstpt_850 = alltracks_850['FIRST_PT'][:]
index_850   = alltracks_850['index'][:]

TrackTimes_850 = alltracks_850['time'][:]
TrackTimeUnits_850 = alltracks_850.variables['time'].units
dates850 = num2date(TrackTimes_850, TrackTimeUnits_850)

lats850 = alltracks_850['latitude'][:]
lons850 = alltracks_850['longitude'][:]



# Then the 700 hPa tracks
trackfile_700 = f'/media/margaret/BigVolume/Data/NASA_PMM_TEW/TRACK_stuff/FilteredNoTCUpdated/700hPa/post_filter.{year}.{hemishort}.NoTC.nc' #'/data3/TRACK/TRACKOutput/MERRA2/FilteredTracks/700hPa/post_filter.{year}.{hemishort}.nc'
alltracks_700 = Dataset(trackfile_700, mode='r')

#the variables
TrackID_700 = alltracks_700['TRACK_ID'][:]
numrecs_700 = alltracks_700.dimensions['record'].size
firstpt_700 = alltracks_700['FIRST_PT'][:]
index_700   = alltracks_700['index'][:]

TrackTimes_700 = alltracks_700['time'][:]
TrackTimeUnits_700 = alltracks_700.variables['time'].units
dates700 = num2date(TrackTimes_700, TrackTimeUnits_700)

lats700 = alltracks_700['latitude'][:]
lons700 = alltracks_700['longitude'][:]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Go through filtered tracks and make a lot of meshgrids to compare
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# the meshgrids to compare the dates, lats, and lons

datemesh700, datemesh850 = np.meshgrid(dates700, dates850)

latmesh700, latmesh850 = np.meshgrid(lats700, lats850)

lonmesh700, lonmesh850 = np.meshgrid(lons700, lons850)

# We have a 500km radius for the lat/lon matches
# so we need something to allow for that before doing the truth values
distances = great_circle(lonmesh700, latmesh700, lonmesh850, latmesh850)


# do we have matches?
# times should be exact, while distance has a 500km radius
datematch = (datemesh700 == datemesh850)
distmatch = (distances <= 500)

# shoutout to stackoverflow for pointing me to np.logical_and()
# Which will return True for when the two arrays in question are True at the same locations
# https://stackoverflow.com/questions/44578571/intersect-two-boolean-arrays-for-true 
matchmesh = np.logical_and(datematch, distmatch)


# then just do some sums to get T/F for the indices where there's a match!
# crap this loses if a point in one track matches to multiple in another track
# which is itself an interesting question
# that I should not be answering at 2am
# but am going to try to anyways

matchindices850 = np.sum(matchmesh, axis=1)
matchindices700 = np.sum(matchmesh, axis=0)

# remove sus matches
# go through each list and find where the sum is greater than 1
# then find the indices in the row where the sus points are
# compare the distances and see which one is closer
# make that be the one that gets saved





for i in range(len(matchindices850)):
    #sus rows will have a sum greater than 1
    if matchindices850[i] > 1:
        #for curiosity's sake, I want to know how many sus points there are
        #print(f'{matchindices850[i]} points matched to a single point')
        
        #grab the row from matchmesh
        sus = matchmesh[i,:]
        
        #then use np.where to get the indices where True is present
        susindex = np.where(sus)[0]
        
        #use *those* indices to get the sus distances
        susdist = distances[i,susindex]
        
        #the shortest distance is the one we want to keep
        #let's get its index with np.argmax() courtesy of method 2 on
        #https://www.geeksforgeeks.org/how-to-find-the-maximum-and-minimum-value-in-numpy-1d-array/
        keeperindex = np.argmin(susdist)
        
        keeper = susindex[keeperindex]
        
        #and now we can make a new row with only one matching point
        replacement = np.ones(len(sus)) * False
        replacement[keeper] = True
        
        matchmesh[i,:] = replacement
        
#remake sums
print('resetting match indices')
matchindices850 = np.sum(matchmesh, axis=1)
matchindices700 = np.sum(matchmesh, axis=0)
        
#now do it again for the 700 hPa list
#I didn't write a function solely because I was too lazy to 

for j in range(len(matchindices700)):
    #sus rows will have a sum greater than 1
    if matchindices700[j] > 1:
        #for curiosity's sake, I want to know how many sus points there are
        #print(f'{matchindices700[j]} points matched to a single point')
        
        #grab the row from matchmesh
        sus = matchmesh[:,j]
        
        #then use np.where to get the indices where True is present
        susindex = np.where(sus)[0]
        
        #use *those* indices to get the sus distances
        susdist = distances[susindex,j]
        
        #the shortest distance is the one we want to keep
        #let's get its index with np.argmax() courtesy of method 2 on
        #https://www.geeksforgeeks.org/how-to-find-the-maximum-and-minimum-value-in-numpy-1d-array/
        keeperindex = np.argmin(susdist)
        
        keeper = susindex[keeperindex]
        
        #and now we can make a new row with only one matching point
        replacement = np.ones(len(sus)) * False
        replacement[keeper] = True
        
        matchmesh[:,j] = replacement


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Now just use matchindices to make the files of matched and unmatched TEWs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# get the indices where matched, not matched
matches = np.where(matchmesh)
match850 = matches[0]
match700 = matches[1]

allindices850 = np.arange(0, len(lats850))
allindices700 = np.arange(0, len(lats700))
nomatch850 = np.delete(allindices850, match850)
nomatch700 = np.delete(allindices700, match700)

# save some csv files because they're easy to read, write, and share
np.savetxt(f"./CollocationUpdated/match_indices_{year}_{hemishort}.csv", np.array((match850, match700)).transpose(), header=f'Paired matched indices for TEW points in the {hemitext} hemisphere for {year} \n 850 hPa, 700 hPa', delimiter=',')
np.savetxt(f'./CollocationUpdated/nomatch_indices_850hPa_{year}_{hemishort}.csv', nomatch850, header=f'Indices for TEW points in the {hemitext} hemisphere for {year} at 850 hPa which did not match to a TEW at 700 hPa', delimiter=',')
np.savetxt(f'./CollocationUpdated/nomatch_indices_700hPa_{year}_{hemishort}.csv', nomatch700, header=f'Indices for TEW points in the {hemitext} hemisphere for {year} at 700 hPa which did not match to a TEW at 850 hPa', delimiter=',')




# # # # # # # # #
# close the files
# # # # # # # # #
alltracks_850.close()
alltracks_700.close()