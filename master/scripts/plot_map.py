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
#mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, PathPatch, Ellipse
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.geodetics.base import kilometer2degrees
from obspy.core import UTCDateTime, read, Stream
from obspy.imaging.beachball import beach, MomentTensor, mt2axes, plot_mt, mt2plane
from scipy.io.netcdf import NetCDFFile
from matplotlib.colors import LinearSegmentedColormap
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
GEMA_PATH = "%s/GEMA/PyGEMA" % (os.getenv("HOME"))
sys.path.append( GEMA_PATH )

from pygema.read.parameters import load_station_metadata
from pygema.core.mysqlDB import select_event_list, select_events_manual_loc

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

"""
dlon = 0.3
dlat = 0.25
elon = -71.3+dlon
wlon = -71.3-dlon
nlat = -37.9+dlat
slat = -37.9-dlat
"""
elon = -71.0
wlon = -71.75
nlat = -37.7
slat = -38.6


fig = plt.figure()
ax1 = plt.axes([0.1, 0.1, 0.9, 0.9])
map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='f', area_thresh=1000000, epsg=4269)
map.drawparallels(np.arange(slat,nlat,0.1), dashes=[1,2], labels=[1,1,0,0], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
map.drawmeridians(np.arange(wlon,elon,0.1), dashes=[1,2], labels=[0,0,1,1], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
#map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)
map.drawcoastlines(linewidth=0.8,zorder=5)
map.drawcountries(linewidth=0.5, linestyle='-',zorder=2)
#map.fillcontinents(color='0.3', lake_color='steelblue',alpha=0.4,zorder=1)

xpixels = 800
dpi = 300
map.arcgisimage(service='World_Shaded_Relief', xpixels=xpixels, dpi=dpi, verbose=False, zorder=1)
#img = map.arcgisimage(service='NatGeo_World_Map', xpixels=xpixels, dpi=dpi, verbose=False, zorder=0)
#img = map.arcgisimage(service='World_Topo_Map', xpixels=xpixels, dpi=dpi, verbose=False, zorder=0)
img = map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=xpixels, dpi=dpi, verbose=False, zorder=0)
#img = map.arcgisimage(service='ESRI_StreetMap_World_2D', xpixels=xpixels, dpi=dpi, verbose=False, zorder=0)
img.set_alpha(0.6)

"""
# add recent seismicity
events_list = select_events_manual_loc( UTCDateTime(2019,4,5), UTCDateTime().now(), table="LOC", max_gap=270, max_depth=50)
for this_event in events_list:
  xi,yi = map(this_event[1], this_event[2])
  evdx_deg = kilometer2degrees(this_event[9])
  evdy_deg = kilometer2degrees(this_event[10])
  xi_new, yi_new = map(this_event[1]+evdx_deg, this_event[2]+evdy_deg)
  xi_err = abs(xi - xi_new)
  yi_err = abs(yi - yi_new)

  ax1.scatter(xi,yi,color='red', s=120, marker='*', zorder=20, alpha=1, clip_on=True, lw=0.7, edgecolors='k')
  text = "\n%s\n%s\n Ml = %.1f" % (this_event[0].strftime("%Y-%m-%d"), this_event[0].strftime("%H:%M:%S"), this_event[8] )
  #ax1.annotate( text, (xi,yi), fontsize=7, ha='center', va='top', fontweight='bold' )
  ax1.errorbar(xi, yi, xerr=xi_err, yerr=yi_err, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=19)


# add historic seismicity
catalog_pygema = select_events_manual_loc( UTCDateTime(1970, 1,1), UTCDateTime(), table="LOC", max_gap=270, max_depth=50)
for event in catalog_pygema:
  evtime = date2num(event[0].datetime)
  evlon = event[1]
  evlat = event[2]
  evdep = event[3]
  evmag = event[8]
  evdx = event[9]
  evdy = event[10]
  evdz = event[11]

  xi,yi = map(evlon, evlat)
  evdx_deg = kilometer2degrees(evdx)
  evdy_deg = kilometer2degrees(evdy)
  xi_new, yi_new = map(evlon+evdx_deg, evlat+evdy_deg)
  xi_err = abs(xi - xi_new)
  yi_err = abs(yi - yi_new)
  ax1.scatter(xi,yi,color='C1', s=evmag**2*2.5, marker='o', zorder=15, alpha=0.7, clip_on=True, lw=0.3, edgecolors='k')
  #ax1.errorbar(xi, yi, xerr=xi_err, yerr=yi_err, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=14)
"""

si,sj = map(-71.611683, -37.910428)
#ax1.annotate('Embalse\nPangue',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
map.scatter(si, sj, c='C1', linewidths=0.3, edgecolors='k', alpha=1, zorder=15, marker='X', s=40, clip_on=True)

si,sj = map(-71.475571, -38.046040)
#ax1.annotate('Embalse\nRalco',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
map.scatter(si, sj, c='C1', linewidths=0.3, edgecolors='k', alpha=1, zorder=15, marker='X', s=40, clip_on=True)


# add volcanoes
vnames = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
xvolc,yvolc = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
xvolc,yvolc = map(xvolc,yvolc)
ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.7,lw=0.,ms=6.,zorder=10,alpha=1, clip_on=True)

"""
for x,y,vname in zip(xvolc,yvolc,vnames):
  if vname == 'Copahue' or vname == 'Tolhuaca':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
  if vname == 'Callaqui' or vname == 'Trolon':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
  if vname == 'Lonquimay':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
"""

# add stations
networks, stations, stlons, stlats, stalts = load_station_metadata()
for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
  if net=='GM':
    x,y = map(float(lon),float(lat))
    ax1.scatter(x,y,marker='s',color='c', s=30, zorder=2000+2, alpha=1, clip_on=True, lw=0.5, edgecolors='k')
    #plt.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
  #elif net=='C1':
  #  x,y = map(float(lon),float(lat))
  #  ax1.scatter(x,y,marker='d',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=2, edgecolors='C0')
  #  ax1.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)


# ADD faults sielfeld et al. 2019
shp = map.readshapefile(GEMA_PATH+'/src/shapes/gerd/fallas_sielfeld_etal_2019_mod', 'fallas', drawbounds=False)
types = np.unique(np.array([ info['Name']  for info, shape in zip(map.fallas_info, map.fallas) ]))
#print("types of faults =", types)
for info, shape in zip(map.fallas_info, map.fallas):
  if info['Name'] == 'LOFS':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='yellow', alpha=1., linestyle='-', linewidth=0.8, zorder=9, clip_on=True)
  elif info['Name'] == 'ATF':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='yellow', alpha=1., linestyle='--', linewidth=0.8, zorder=9, clip_on=True)



s1 = plt.scatter([], [], c='c', linewidths=0.5, edgecolors='k', alpha=1, zorder=1000, marker='s', s=25)
s2, = plt.plot([], [], c='None', markeredgecolor='C3', markeredgewidth=1.4, lw=0., alpha=1, zorder=1000, marker='^', ms=5)
f1, = plt.plot([-100,-99], [-100,-98], c='yellow', linewidth=1.2, alpha=1, zorder=1000, ls='-') # darkgoldenrod
f2, = plt.plot([-100,-99], [-100,-98], c='yellow', linewidth=1.2, alpha=1, zorder=1000, ls='--') # darkgoldenrod
old_events = plt.scatter([],[], marker='o',color='C1', s=3.0**2*2.5, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
actual_event = plt.scatter([],[], marker='*',color='red', s=90, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')


leg = plt.legend([s2, s1, f1, f2], 
                 [ 'Volcán', 'Estación sismológica', 'Sistema de Fallas\nLiquiñe-Ofqui', 'Falla Biobio-Aluminé' ], 
                 fontsize=6, ncol=1, frameon=True, fancybox=True, shadow=False, framealpha=0.3, loc=4 ) 
leg.set_zorder(1000)


figname = "figs/map.jpg" 
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close("all")
















