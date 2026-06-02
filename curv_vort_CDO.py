##### BASIC DESCRIPTION #####
# Code for calculating curvature vorticity on a pressure level from gridded wind components
# By Margaret Hollis
# Core calculations written Fall 2019
# Updated Jan 2020 to process a full year at a time and write to a netCDF file

import glob
import sys
import numpy as np
from netCDF4 import num2date, Dataset, date2num
import netCDF4
import metpy
import metpy.calc as mpcalc
import subprocess

### THINGS THAT NEED TO BE UPDATED DEPENDING ON FILES, LEVEL, ETC. ###

#what year and level are we running today?
year = int(sys.argv[1])
print(f'The provided year is {year}. If this is not the year you wish to run, please quit now')
level = int(sys.argv[2])
dataset = 'MERRA2'

#where's the data to be processed located at?
#this should be a list of filepaths that can be iterated through
#filepaths = sorted(glob.glob('/data/deluge/reanalysis/REANALYSIS/MERRA2/3D_PL/native/'+str(year)+'/*/*.nc4'))
filepaths = sorted(glob.glob(f'/data/deluge/reanalysis/REANALYSIS/MERRA2/3D_PL/native/{year}/*/*.nc4'))

#where to save the output to?
#this should be a string and should end with '/' so that files can be saved into that folder
outpath = '/data2/Data_Processed/MERRA2-Margaret/'

# # Shouldn't need to change things after this, but read assumptions below to be sure ###
# 1 Indexing follows this order:
#   [time, level, lat, lon]
# 2 Levels are consistent from file to file (index of level of interest does not change)
# 3 The lat lon grid is consistent from file to file
# 4 The number of times in the file is the same from file to file (day lengths are the same)
#   Basically this code assumes that things are consistent from file to file

#counter is mostly just for setting up the output file
#but being able to print "progress reports" is useful too
count = 0

for filepath in filepaths:
	#First off, the stuff that needs to be done once per file
	file = netCDF4.Dataset(filepath, mode='r')
	#print(f'reached file {filepath}')
	
	if count == 0:
	
		#Find the index of the level of interest from the list of levels
		levels = file.variables['lev'][:]
		
		levindex = 0
		for lev in range(len(levels)):
			if int(levels[lev]) == level:
				levindex = lev
				print(f'found the index of {level} hPa level')
				break
			elif lev+1 == len(levels):
				print('Error: The given level does not exist in the dataset')
				levindex = np.nan
			else:
				continue
	
	#pull out U, V wind components at all times and the appropriate level
	#Order of indices may need to be changed if data isn't in the assumed order, listed above
	U = file.variables['U'][:,levindex,:,:]
	V = file.variables['V'][:,levindex,:,:]

	#time stuff
	times = file.variables['time'][:]
	timeUnits = file.variables['time'].units
	date_time = num2date(times,units=timeUnits)

	#before we do any actual calculations, set up output file
	#will also read in lat lon lists here, since those aren't changing
	if count == 0:
		#lat lon lists
		lat = file.variables['lat'][:]
		lon = file.variables['lon'][:]
		
		#the original x and y dimensions will be useful later on with the wrapping
		x = len(lon)
		y = len(lat)
        
		#create the file, open it		
		outfile = f'{outpath}CurvVort/cur_vort.{dataset}.{level}.{year}.nc'
		cvort_file = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
		
		#create dimensions
		#follows file.createDimension("dim name", var size or None for unlimited)
		out_time = cvort_file.createDimension("time", None)
		out_lon = cvort_file.createDimension("lon", x)
		out_lat = cvort_file.createDimension("lat", y)

		#create variables
		#follows file.createVariable("var name", datatype, dimension (leave off if creating scalar))
		out_cv = cvort_file.createVariable("cv", "f8", ("time", "lat", "lon",))
		out_times = cvort_file.createVariable("time", "i8", ("time",))
		out_lats = cvort_file.createVariable("lat", "d", ("lat"))
		out_lons = cvort_file.createVariable("lon", "d", ("lon"))
		out_lev = cvort_file.createVariable("level", "d")

		#variable attributes
		#units, long name
		out_cv.units = "1/s"
		out_cv.long_name = "curvature_vorticity"

		#the units for the time outputs need to be consistent
		#but the numbers will get *really* big if we just use the "minutes since" date units
		#so we'll have the output units be hours since start of file
		out_timeUnits = f'hours since {date_time[0]}'
		out_times.units = out_timeUnits
		out_times.long_name = "time"

		out_lats.units = "degrees_north"
		out_lats.long_name = "latitude"

		out_lons.units = "degrees_east"
		out_lons.long_name = "longitude"

		out_lev.units = "hPa"
		out_lev.long_name = "level"

		#global attributes
		cvort_file.description = f'{dataset} curvature vorticity at {level} hPa for {year}'
		
		#we can go ahead and write the level, lats, and lons now.
		out_lev.assignValue(level)
		out_lons[:] = lon
		out_lats[:] = lat
		
		#the number of times in the file is also useful
		numtimes = len(times)
	
	#close the data file, because we've gotten out what we need
	file.close()

	#wrapping of lat/lons only needs to be done once as long as data grid is consistent
	if count == 0:
		#wrap and set up meshgrid now that the output file has been set up
		#wrap two indices in each direction
		#don't wrap lat because that's not physical
		lon = np.pad(lon, 2, 'wrap')
		lons, lats = np.meshgrid(lon, lat)
	
	#now for loops for the hours
	for h in range(len(times)):
		#arrays are [time, level, lat, lon] and levels are done bottom up
		u = U[h,:,:]
		v = V[h,:,:]
		
		#because we want centered derivatives everywhere, let's wrap all of our arrays
		#so that there's a bit of overlap at the 180/-180 line
        #only doing a zonal wrap because I don't currently care about the poles
		u = u.take(range(-2,x+2), axis=1, mode='wrap')
		v = v.take(range(-2,x+2), axis=1, mode='wrap')

		#metpy wants things to have units, so let's attach them to u, v
		u = metpy.units.masked_array(u, 'meter/second')
		v = metpy.units.masked_array(v, 'meter/second')
        
		#will need to define dx, dy for the metpy module
		#mpcalc.lat_lon_grid_deltas will return dx, dy in meters
		#even when given data on degree grids
		dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

		#start off with relative vorticity (relvort)
		#if this thing throws errors,
		#mpcalc.vorticity needs the units to be attached to the arrays

		relvort = mpcalc.vorticity(u, v, dx=dx, dy=dy) #arguments are u array, v array, dx, dy

		#then calculate shear vorticity (shrvort)
		#finally figured out directional derivatives using this page:
		#http://mathworld.wolfram.com/DirectionalDerivative.html
		
		wind = np.sqrt(u**2 + v**2)

		dVdx = mpcalc.first_derivative(wind, axis=1, delta=dx)
		dVdy = mpcalc.first_derivative(wind, axis=0, delta=dy)

		#and now we can *finally* calculate shear vorticity
		shrvort = -((dVdx * (-v/wind)) + (dVdy * (u/wind)))

		#finally do difference, since curvature vort = rel vort - shear vort
		curvort = relvort - shrvort

		#for saving, will want to bring curvature vort back to going just from -180 to 180
		#lons were already saved in original dimensions back when setting up the file
		curvort = curvort[:,2:x+2]
		
		#and write the curvature vorticity and its time to the netCDF file
		
		out_times[(numtimes*count+h)] = date2num(date_time[h], out_timeUnits)
		out_cv[(numtimes*count+h),:,:] = np.ma.array(curvort.data, mask = curvort.mask)


	#progress update
	if (count+1)%10 == 0:
		print(f'successfully made it through file number {count}')
		
	count += 1
	
#close curvature vort file when done with everything else
print(f'Successfully completed processing {year}. Now closing output file and starting CDO processing')
cvort_file.close()

subprocess.run(['cdo','setrtomiss,1.0e+10,1.0e+38',
                outfile,
                f'{outpath}setmiss/cur_vort.{dataset}.{level}.{year}.new.nc'])

subprocess.run(['cdo','fillmiss',
                f'{outpath}setmiss/cur_vort.{dataset}.{level}.{year}.new.nc',
                f'{outpath}filled/cur_vort.{dataset}.{level}.{year}.filled.nc'])