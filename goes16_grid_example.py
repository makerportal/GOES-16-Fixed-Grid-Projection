#!/usr/bin/python

from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import os

def lat_lon_reproj(nc_folder,nc_indx):
    os.chdir(nc_folder)
    full_direc = os.listdir()
    nc_files = [ii for ii in full_direc if ii.endswith('.nc')]
    g16_data_file = nc_files[nc_indx] # select .nc file
    print(nc_files[nc_indx]) # print file name

    # designate dataset
    g16nc = Dataset(g16_data_file, 'r')
    var_names = [ii for ii in g16nc.variables]
    var_name = var_names[0]
    try:
        band_id = g16nc.variables['band_id'][:]
        band_id = ' (Band: {},'.format(band_id[0])
        band_wavelength = g16nc.variables['band_wavelength']
        band_wavelength_units = band_wavelength.units
        band_wavelength_units = '{})'.format(band_wavelength_units)
        band_wavelength = ' {0:.2f} '.format(band_wavelength[:][0])
        print('Band ID: {}'.format(band_id))
        print('Band Wavelength: {} {}'.format(band_wavelength,band_wavelength_units))
    except:
        band_id = ''
        band_wavelength = ''
        band_wavelength_units = ''

    # GOES-R projection info and retrieving relevant constants
    proj_info = g16nc.variables['goes_imager_projection']
    lon_origin = proj_info.longitude_of_projection_origin
    H = proj_info.perspective_point_height+proj_info.semi_major_axis
    r_eq = proj_info.semi_major_axis
    r_pol = proj_info.semi_minor_axis

    # grid info
    lat_rad_1d = g16nc.variables['x'][:]
    lon_rad_1d = g16nc.variables['y'][:]

    # data info    
    data = g16nc.variables[var_name][:]
    data_units = g16nc.variables[var_name].units
    data_time_grab = ((g16nc.time_coverage_end).replace('T',' ')).replace('Z','')
    data_long_name = g16nc.variables[var_name].long_name
    
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

    return lon,lat,data,data_units,data_time_grab,data_long_name,band_id,band_wavelength,band_wavelength_units,var_name

nc_folder = './rad_nc_files/' # define folder where .nc files are located

file_indx = 15 # be sure to pick the correct file. Make sure the file is not too big either,
# some of the bands create large files (stick to band 7-16)

# main data grab from function above
lon,lat,data,data_units,data_time_grab,data_long_name,band_id,band_wavelength,band_units,var_name = lat_lon_reproj(nc_folder,file_indx)

bbox = [np.min(lon),np.min(lat),np.max(lon),np.max(lat)] # set bounds for plotting

# figure routine for visualization
fig = plt.figure(figsize=(9,4),dpi=200)

n_add = 0
m = Basemap(llcrnrlon=bbox[0]-n_add,llcrnrlat=bbox[1]-n_add,urcrnrlon=bbox[2]+n_add,urcrnrlat=bbox[3]+n_add,resolution='i', projection='cyl')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.25)
m.pcolormesh(lon.data, lat.data, data, latlon=True)
parallels = np.linspace(np.min(lat),np.max(lat),5.)
m.drawparallels(parallels,labels=[True,False,False,False])
meridians = np.linspace(np.min(lon),np.max(lon),5.)
m.drawmeridians(meridians,labels=[False,False,False,True])
cb = m.colorbar()

data_units = ((data_units.replace('-','^{-')).replace('1','1}')).replace('2','2}')
plt.rc('text', usetex=True)
cb.set_label(r'%s $ \left[ \mathrm{%s} \right] $'% (var_name,data_units))
plt.title('{0}{2}{3}{4} on {1}'.format(data_long_name,data_time_grab,band_id,band_wavelength,band_units))
plt.tight_layout()

#plt.savefig('goes_16_demo.png',dpi=200,transparent=True) # uncomment to save figure
plt.show()
