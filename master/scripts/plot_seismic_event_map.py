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


max_gap = 270
max_depth = 50

events_list = select_events_manual_loc( UTCDateTime(1970, 1,1), UTCDateTime().now(), table="LOC", max_gap=max_gap, max_depth=max_depth)

events_list = np.array(events_list)
events_list = events_list[np.argsort(events_list[:, 0])]

# print events found on the screen
count = 1; count_list = []
print("\n\n")
for event in events_list:
  count_list.append(str(count))
  pattern = "[%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f %s " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S.%f"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
  print(pattern)
  count += 1

print("(from PyGEMA database: Table = LOC; max_gap = %.1f deg; max_depth = %.1f km)" % (max_gap, max_depth) )

#print( "[%i] TODOS!!" % (count) )
flag = input("\n+ Type the seismic event that you want to plot: ")
while not flag in count_list:
  flag = input("+ Type the seismic event that you want to plot: ")
  if flag in count_list:
    break

this_event = events_list[int(flag)-1]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

dlon = 0.3
dlat = 0.3

elon = this_event[1]+dlon
wlon = this_event[1]-dlon
nlat = this_event[2]+dlat
slat = this_event[2]-dlat

#elon = -71.183+dlon
#wlon = -71.183-dlon
#nlat = -37.856+dlat
#slat = -37.856-dlat

fig = plt.figure()
ax1 = plt.axes([0.1, 0.1, 0.9, 0.9])

map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='c', area_thresh=1000000, epsg=4269)
map.drawparallels(np.arange(slat,nlat,dlat/5.), dashes=[1,2], labels=[1,1,0,0], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
map.drawmeridians(np.arange(wlon,elon,dlon/5.), dashes=[1,2], labels=[0,0,1,1], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
#map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)

map.drawcoastlines(linewidth=0.8,zorder=5)
#map.drawcountries(linewidth=0.5, linestyle=':',zorder=2)
map.fillcontinents(color='w', lake_color='steelblue',alpha=0.,zorder=1)
#map.arcgisimage(service='World_Shaded_Relief', xpixels=1920, dpi=300, verbose=False, zorder=0)
map.arcgisimage(service='NatGeo_World_Map', xpixels=1920, dpi=300, verbose=False, zorder=0)
#map.arcgisimage(service='World_Topo_Map', xpixels=1920, dpi=300, verbose=False, zorder=0)
#map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=1920, dpi=300, verbose=False, zorder=0)
#map.arcgisimage(service='ESRI_StreetMap_World_2D', xpixels=1920, dpi=300, verbose=False, zorder=0)

# add current event
xi,yi = map(this_event[1], this_event[2])
evdx_deg = kilometer2degrees(this_event[9])
evdy_deg = kilometer2degrees(this_event[10])
xi_new, yi_new = map(this_event[1]+evdx_deg, this_event[2]+evdy_deg)
xi_err = abs(xi - xi_new)
yi_err = abs(yi - yi_new)

ax1.scatter(xi,yi,color='red', s=120, marker='*', zorder=20, alpha=1, clip_on=True, lw=0.7, edgecolors='k')
text = "%s\n" % (this_event[0].strftime("%Y-%m-%d %H:%M:%S UTC") )
ax1.annotate( text, (xi,yi), fontsize=6, ha='left', va='bottom', fontweight='bold' )
ax1.errorbar(xi, yi, xerr=xi_err, yerr=yi_err, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=19)



# add old seismicity
catalog_pygema = select_events_manual_loc( UTCDateTime(1970, 1,1), this_event[0], table="LOC", max_gap=max_gap, max_depth=max_depth)
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


si,sj = map(-71.611683, -37.910428)
ax1.annotate('Embalse\nPangue',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
map.scatter(si, sj, c='k', linewidths=0.5, edgecolors='k', alpha=1, zorder=15, marker='X', s=15, clip_on=True)

si,sj = map(-71.475571, -38.046040)
ax1.annotate('Embalse\nRalco',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
map.scatter(si, sj, c='k', linewidths=0.5, edgecolors='k', alpha=1, zorder=15, marker='X', s=15, clip_on=True)


# add volcanoes
vnames = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
xvolc,yvolc = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
xvolc,yvolc = map(xvolc,yvolc)
ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.1,lw=0.,ms=5.,zorder=10,alpha=1, clip_on=True)

for x,y,vname in zip(xvolc,yvolc,vnames):
  if vname == 'Copahue' or vname == 'Tolhuaca':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
  if vname == 'Callaqui' or vname == 'Trolon':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
  if vname == 'Lonquimay':
    ax1.annotate('V. '+vname, (x, y), color='k', weight='normal', fontsize=6, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')


# add stations
networks, stations, stlons, stlats, stalts = load_station_metadata()
for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
  if net=='GM':
    x,y = map(float(lon),float(lat))
    ax1.scatter(x,y,marker='s',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=2, edgecolors='C0')
    plt.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
  elif net=='C1':
    x,y = map(float(lon),float(lat))
    ax1.scatter(x,y,marker='d',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=2, edgecolors='C0')
    ax1.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)


# ADD faults sielfeld et al. 2019
shp = map.readshapefile(GEMA_PATH+'/src/shapes/gerd/fallas_sielfeld_etal_2019_mod', 'fallas', drawbounds=False)
types = np.unique(np.array([ info['Name']  for info, shape in zip(map.fallas_info, map.fallas) ]))
#print("types of faults =", types)
for info, shape in zip(map.fallas_info, map.fallas):
  if info['Name'] == 'LOFS':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='k', alpha=1., linestyle='-', linewidth=0.5, zorder=9, clip_on=True)
  elif info['Name'] == 'ATF':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='k', alpha=1., linestyle='--', linewidth=0.5, zorder=9, clip_on=True)


s1 = plt.scatter([], [], c='None', linewidths=1.1, edgecolors='C0', alpha=1, zorder=1000, marker='s', s=15)
s2, = plt.plot([], [], c='None', markeredgecolor='C3', markeredgewidth=1.1, lw=0., alpha=1, zorder=1000, marker='^', ms=5)
f1, = plt.plot([-100,-99], [-100,-98], c='k', linewidth=.75, alpha=1, zorder=1000, ls='-') # darkgoldenrod
f2, = plt.plot([-100,-99], [-100,-98], c='k', linewidth=.75, alpha=1, zorder=1000, ls='--') # darkgoldenrod
old_events = plt.scatter([],[], marker='o',color='C1', s=3.0**2*2.5, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
actual_event = plt.scatter([],[], marker='*',color='red', s=90, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')


leg = plt.legend([actual_event, old_events, s1, s2, f1, f2], 
                 ['Sismicidad reciente', 'Sismicidad histórica', 'Estación sismológica', 'Volcán Holoceno', 'Sistema de Fallas\nLiquiñe-Ofqui', 'Falla Biobio-Aluminé' ], 
                 fontsize=6, ncol=1, frameon=True, fancybox=True, shadow=False, framealpha=0.8, loc=2 ) 
leg.set_zorder(1000)



figname = "figs/event_map.jpg" 
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close("all")



