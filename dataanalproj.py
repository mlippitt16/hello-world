#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:40:30 2023

@author: maureenlippitt
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
#import metpy.calc as mp



#Open the data files
ds = xr.open_dataset('CAM_FV3L30_HS_month_12_3D.nc')
ds_m = xr.open_dataset('CAM_FV3L30_moist_HS_month_12_3D.nc')
#Open the U variable, selecting for levels 0-29. 
U = ds['U'].isel(lev = slice(0, 29))
U_m = ds_m['U'].isel(lev = slice(0, 29))
#Select for just one time value to eliminate the time dimension
U = U.isel(time = 0)
U_m = U_m.isel(time = 0)

'''
#I DONT UNDERSTAND HOW TO FIND GEOSTROPHIC WIND
'''
#Calculate geostrophic wind
def geostrophic(U=U):
    #Pull out arrays from U
    lat = U['lat']
    lon = U['lon']
    lev = U['lev']
    
    #Define constants
    g = 9.81 #m/s
    f = 2 * 7.2921159e-5 * np.sin(np.radians(lat))  # Coriolis parameter (1/s)
    a = 6371000.0  # Earth's radius in meters
   
    # Calculate the geostrophic wind
    
    
    
    return geostrophic_wind














'''
QUESTION TWO
'''

#2A
#Define important constants
R = 287.05  # Specific gas constant for dry air
P0 = ds['P0'] / 100  # Reference pressure in hPa
Cp = 1004 # Specific heat capacity for dry air (J/(kg*K))
T = ds['T']  #Temperature data
P = ds['lev'] #Pressure data in hPa

#Potential temperature equation
#This creates a 3D array of potential temperatures  
pot_temp = T * (P0 / P) ** (R / Cp)
#Remove the time variable by selecting one timestamp
pot_temp = pot_temp.isel(time=0)
#Average over all longitudes to get the zonal mean
zonal_mean = pot_temp.mean(dim='lon')

#PLOTTING
# Define custom contour levels as specified in problem
contour_levels = np.concatenate([np.arange(260, 350, 10), np.arange(400, 700, 100)])
# Create a contour plot
fig, ax = plt.subplots(figsize = (12,5))
plt.contourf(zonal_mean['lat'], zonal_mean['lev'], zonal_mean, levels=contour_levels, extend='both', cmap='rainbow')
plt.colorbar(label='Potential Temperature (K)') #Plot colorbar
plt.title('Zonal-Mean Potential Temperature') #Plot title
plt.xlabel('Latitude (degrees)') #Plot x axis label
plt.ylabel('Pressure (hPa)') #Plot y axis label
plt.gca().invert_yaxis() # Reverse the y-axis
# Add contour lines for the specified contour levels
contours = plt.contour(zonal_mean['lat'], zonal_mean['lev'], zonal_mean, levels=contour_levels, colors='k')
plt.clabel(contours, inline=True, fmt='%1.0f K', fontsize=10)
# Show the plot
plt.show()


#2D
# Define custom contour levels for the potential temperature.
contour_levels = [260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600, 650, 700]
# Create a contour plot with custom contour levels and a non-equidistant color scheme
fig, ax = plt.subplots(figsize = (12,5))
plt.contourf(zonal_mean['lat'], np.log10(zonal_mean['lev']), zonal_mean, levels=contour_levels, extend='both', cmap='rainbow')
plt.colorbar(label='Potential Temperature (K)')
plt.title('Zonal-Mean Potential Temperature')
plt.xlabel('Latitude (degrees)')
plt.ylabel('log10(Pressure) (log10(hPa)')
plt.gca().invert_yaxis() # Reverse the y-axis
# Add contour lines for the specified contour levels
contours = plt.contour(zonal_mean['lat'], np.log10(zonal_mean['lev']), zonal_mean, levels=contour_levels, colors='k')
plt.clabel(contours, inline=True, fmt='%1.0f K', fontsize=10)
# Show the plot
plt.show()











'''
QUESTION THREE
'''

















