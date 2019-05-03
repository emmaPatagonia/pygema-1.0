 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) Diego Gonzalez-Vidal <diegogonzalezvidal@gmail.com>
#
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with This program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import mpl_toolkits
mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, PathPatch, Ellipse
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.core import UTCDateTime, read, Stream
from obspy.imaging.beachball import beach, MomentTensor, mt2axes, plot_mt, mt2plane
from scipy.io.netcdf import NetCDFFile
from matplotlib.colors import LinearSegmentedColormap
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging


from pygema.plot.GMTcptcity import gmtColormap
from pygema.plot.shading import set_shade, hillshade 
from pygema.read.parameters import find_pygema_parent_directory


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def reverse_colourmap(cmap, name = 'my_cmap_r'):
  reverse = []
  k = []   
  for key in cmap._segmentdata:    
    k.append(key)
    channel = cmap._segmentdata[key]
    data = []
    for t in channel:                    
      data.append((1-t[0],t[2],t[1]))            

    reverse.append(sorted(data))    

  LinearL = dict(zip(k,reverse))
  my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
  return my_cmap_r


def add_dem_srtm1_topo(map, netcdffile, wlon, elon, slat, nlat):
  pygema_path = find_pygema_parent_directory()
  netcdf = NetCDFFile(netcdffile)
  topoin = netcdf.variables['dem30_no0_biobio_ProjectRast'][:]
  lon = netcdf.variables['lon'][:]
  lat = netcdf.variables['lat'][:]
  plon = np.where( (lon>=wlon) & (lon<=elon) )[0]
  plat = np.where( (lat>=slat) & (lat<=nlat) )[0]
  elev = topoin[plat].T[plon].T
  elev[elev==-32767] = np.ma.masked
  elev[elev<0] = np.ma.masked
  left, bottom = map(wlon,slat)
  right, top = map(elon,nlat)
  colormap = pygema_path+'/src/ib12.cpt' 
  cdict = gmtColormap(colormap)
  cmap1 = LinearSegmentedColormap('my_colormap',cdict,256)
  cmap1 = reverse_colourmap(cmap1)
  cmap1.set_bad('w',alpha=0.)
  intensity = hillshade(elev,scale=20,azdeg=290.0,altdeg=45.0)
  rgb = set_shade(elev,intensity=intensity,cmap=cmap1) 
  map.imshow(rgb, cmap=cmap1, zorder = 1, origin='upper',clip_on=True, extent=[left, right, bottom, top])
  colormap = pygema_path+'/src/DEM_print.cpt'
  cdict = gmtColormap(colormap) 
  cmap2 = LinearSegmentedColormap('my_colormap',cdict,256)
  cmap2.set_bad('w',alpha=0.)
  #map.imshow(elev, cmap=cmap2, zorder = 1,origin='upper',alpha=0.3,clip_on=True, extent=[left, right, bottom, top])
  return map