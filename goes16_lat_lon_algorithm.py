#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os

# navigate to directory with .nc data files (raw radiance files in my case)
os.chdir('./rad_nc_files/')
full_direc = os.listdir()
nc_files = [ii for ii in full_direc if ii.endswith('.nc')]
nc_indx = 5
g16_data_file = nc_files[nc_indx] # select .nc file
print(nc_files[nc_indx]) # print file name

# designate dataset
g16nc = Dataset(g16_data_file, 'r')

# GOES-R projection info and retrieving relevant constants
proj_info = g16nc.variables['goes_imager_projection']
lon_origin = proj_info.longitude_of_projection_origin
H = proj_info.perspective_point_height+proj_info.semi_major_axis
r_eq = proj_info.semi_major_axis
r_pol = proj_info.semi_minor_axis

# Data info
lat_rad_1d = g16nc.variables['x'][:]
lon_rad_1d = g16nc.variables['y'][:]

# close file when finished
g16nc.close()
g16nc = None

# create meshgrid filled with radian angles
lat_rad,lon_rad = np.meshgrid(lat_rad_1d,lon_rad_1d)

# lat/lon calc routine from satellite radian angle vectors

lambda_0 = (lon_origin*np.pi)/180.0

a_var = np.power(np.sin(lat_rad),2.0) + (np.power(np.cos(lat_rad),2.0)*(np.power(np.cos(lon_rad),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(lon_rad),2.0))))
b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
c_var = (H**2.0)-(r_eq**2.0)

r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)

s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
s_y = - r_s*np.sin(lat_rad)
s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)

lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)


# print test coordinates
print('{} N, {} W'.format(lat[318,1849],abs(lon[318,1849])))
