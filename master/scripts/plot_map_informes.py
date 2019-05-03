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


# create map
elon = -70.8
wlon = -71.9
nlat = -37.65
slat = -38.55

fig = plt.figure()
ax1 = plt.axes([0.1, 0.3, 0.9, 0.6])
map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='f', area_thresh=1000000, epsg=4269)
map.drawparallels(np.arange(slat,nlat,0.1), dashes=[1,2], labels=[1,0,0,0], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
map.drawmeridians(np.arange(wlon,elon,0.1), dashes=[1,2], labels=[0,0,1,0], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
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


# add places
si,sj = map(-71.611683, -37.910428)
ax1.annotate('Pangue ',xy=(si,sj), ha='right', va='top', color='k',zorder=15,fontsize=6,fontweight='bold',xytext=(-0,-0),textcoords='offset points', fontstyle='normal', clip_on=True)
map.scatter(si, sj, c='C1', linewidths=0.3, edgecolors='k', alpha=1, zorder=15, marker='X', s=40, clip_on=True)

si,sj = map(-71.475571, -38.046040)
ax1.annotate('Ralco ',xy=(si,sj), ha='right', va='top', color='k',zorder=15,fontsize=6,fontweight='bold',xytext=(-0,-0),textcoords='offset points', fontstyle='normal', clip_on=True)
map.scatter(si, sj, c='C1', linewidths=0.3, edgecolors='k', alpha=1, zorder=15, marker='X', s=40, clip_on=True)

# add volcanoes
vnames = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
xvolc,yvolc = np.loadtxt(GEMA_PATH+'/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
xvolc,yvolc = map(xvolc,yvolc)
ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.7,lw=0.,ms=6.,zorder=10,alpha=1, clip_on=True)

for x,y,vname in zip(xvolc,yvolc,vnames):
  if vname == 'Tolhuaca':
    ax1.annotate('V. '+vname+'  ', (x, y), color='k', weight='bold', fontsize=5.5, ha='right', va='center', clip_on=True, zorder=2000+3, fontstyle='normal')
  elif vname == 'Copahue':
    ax1.annotate('   V. '+vname, (x, y), color='k', weight='bold', fontsize=5.5, ha='left', va='top', clip_on=True, zorder=2000+3, fontstyle='normal')
  elif vname == 'Callaqui':
    ax1.annotate(' V. '+vname, (x, y), color='k', weight='bold', fontsize=5.5, ha='left', va='bottom', clip_on=True, zorder=2000+3, fontstyle='normal')
  elif vname == 'Trolon':
    ax1.annotate('V. '+vname+' ', (x, y), color='k', weight='bold', fontsize=5.5, ha='right', va='bottom', clip_on=True, zorder=2000+3, fontstyle='normal')
  elif vname == 'Lonquimay':
    ax1.annotate('  V. '+vname, (x, y), color='k', weight='bold', fontsize=5.5, ha='left', va='bottom', clip_on=True, zorder=2000+3, fontstyle='normal')

# add stations
networks, stations, stlons, stlats, stalts = load_station_metadata()
for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
  if net=='GM':
    x,y = map(float(lon),float(lat))
    ax1.scatter(x,y,marker='s',color='c', s=30, zorder=2000+2, alpha=1, clip_on=True, lw=0.5, edgecolors='k')
    #plt.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
  elif net=='C1':
    x,y = map(float(lon),float(lat))
    ax1.scatter(x,y,marker='d',color='None', s=20, zorder=2000+2, alpha=1, clip_on=True, lw=1.7, edgecolors='C0')
    #ax1.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)

# ADD faults sielfeld et al. 2019
shp = map.readshapefile(GEMA_PATH+'/src/shapes/gerd/fallas_sielfeld_etal_2019_mod', 'fallas', drawbounds=False)
types = np.unique(np.array([ info['Name']  for info, shape in zip(map.fallas_info, map.fallas) ]))
#print("types of faults =", types)
for info, shape in zip(map.fallas_info, map.fallas):
  if info['Name'] == 'LOFS':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='yellow', alpha=1., linestyle='-', linewidth=0.75, zorder=9, clip_on=True)
  elif info['Name'] == 'ATF':
    x, y = zip(*shape) 
    map.plot(x, y, marker=None, color='yellow', alpha=1., linestyle='--', linewidth=0.75, zorder=9, clip_on=True)


# add seismicity gema
max_gap = 360
max_depth = 250
events_list = select_events_manual_loc( UTCDateTime(1970, 1,1), UTCDateTime(), table="LOC", max_gap=max_gap, max_depth=max_depth)

zmin = 0
zmax = 25
cmap = plt.cm.gnuplot
nlevels = 20
bounds = np.linspace(zmin, zmax,nlevels, endpoint=True)
levels = np.linspace(zmin, zmax,nlevels, endpoint=True)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

x = []; y = []; z = []; s = []; xerr = []; yerr = []; zerr = []
for event in events_list:
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

  x.append(xi); y.append(yi); z.append(evdep); s.append( evmag**3*1.5 ); xerr.append(xi_err); yerr.append(yi_err); zerr.append(evdz)

ax1.scatter(x,y,c=z, cmap=cmap, s=s, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
#ax1.errorbar(x, y, xerr=xerr, yerr=yerr, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=19)


# add seismicity and focal mechanism from cgmt and seielfeld

catalog_gcmt = np.loadtxt(GEMA_PATH+"/src/global_cmt.txt", usecols=(0,1,2,3,4,5,6,7,8,9) )

fm_xi = []; fm_yi = []; fm_colors = []; fm_lons = []; fm_lats = []; fm_depths = []; fm_mags = []; fm_strikes = []; fm_dips = []; fm_rakes = []
evlon = catalog_gcmt[0]
evlat = catalog_gcmt[1]
xi, yi = map(evlon, evlat)
evdep = catalog_gcmt[2]
evmag = 5.5 #catalog_gcmt[3]
mt = MomentTensor(catalog_gcmt[3:9], catalog_gcmt[9] )
plane = mt2plane(mt)
strike = plane.strike
dip = plane.dip
rake = plane.rake
if rake >= 180: rake = rake - 360
fm_xi.append(xi)
fm_yi.append(yi)
fm_lons.append(evlon)
fm_lats.append(evlat)
fm_depths.append(evdep)
fm_mags.append(evmag)
fm_strikes.append(strike)
fm_dips.append(dip)
fm_rakes.append(rake)
fm_colors.append('purple')

catalog_sielfeld2019 = np.loadtxt(GEMA_PATH+"/src/sielfeld_2019.txt", usecols=(2, 1, 3, 9, 4, 5, 6) )
for event in catalog_sielfeld2019:
  evlon = event[0]
  evlat = event[1]
  xi, yi = map(evlon, evlat)
  evdep = event[2]
  evmag = event[3]
  strike = event[4]
  dip = event[5]
  rake = event[6]
  if rake >= 180: rake = rake - 360

  fm_xi.append(xi)
  fm_yi.append(yi)
  fm_lons.append(evlon)
  fm_lats.append(evlat)
  fm_depths.append(evdep)
  fm_mags.append(evmag)
  fm_strikes.append(strike)
  fm_dips.append(dip)
  fm_rakes.append(rake)
  fm_colors.append('c')

ax1.scatter(fm_xi,fm_yi,c=fm_depths, cmap=cmap, s=np.array(fm_mags)**3*1.5, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
ax1.scatter(fm_xi,fm_yi,color='k', s=2, marker='o', zorder=21, alpha=0.8, clip_on=True, lw=0., edgecolors='k')

"""
lon1 = -71.72
lat1 = -37.92
increment_lon = -0.03
increment_lat = -0.1
for xi,yi,strike,dip,rake,color,mag in zip(fm_xi, fm_yi, fm_strikes, fm_dips, fm_rakes, fm_colors, fm_mags):
  print(xi,yi,strike,dip,rake,color,mag)
  x, y = map(lon1,lat1)
  xi, yi = map(xi, yi)
  b1 = beach([strike,dip,rake], xy=(x,y), width=6000, linewidth=0.5, facecolor=color, alpha=.9)      
  b1.set_zorder(10)
  b1.set_clip_on(True)
  ax1.add_collection(b1)
  ax1.plot([x,xi],[y,yi], lw=0.4, color='0.2', zorder=9 )
  ax1.annotate('M = %.1f'%(mag), (x,y), textcoords='offset points', size=4, xytext=(-9, 7), zorder=10000, weight='normal', clip_on=True)
  lon1+=increment_lon
  lat1+=increment_lat
"""


# add legend
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# add axis at bottom

ax2 = plt.axes([0.275, 0.1, 0.55, 0.17])
ax2.minorticks_on()
ax2.tick_params(axis='both', which='major', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in')
ax2.tick_params(axis='both', which='minor', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in')
ax2.tick_params(axis='x', which='major', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
ax2.tick_params(axis='x', which='minor', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
ax2.spines['right'].set_visible(True)
ax2.spines['top'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax2.spines['bottom'].set_visible(True)
ax2.set_xlim(wlon, elon)
ax2.set_ylim(zmax, zmin)
labels = [ r"%.1f$^{\circ}$" % (item) for item in ax2.get_xticks().tolist() ]
ax2.set_xticklabels(labels)
ax2.set_ylabel("Prof. (km)", fontsize=6)

# add seismic events pygema
x = []; y = []; z = []; s = []; xerr = []; yerr = []; zerr = []
for event in events_list:
  evlon = event[1]
  evlat = event[2]
  evdep = event[3]
  evmag = event[8]
  evdx = event[9]
  evdy = event[10]
  evdz = event[11]

  evdx_deg = kilometer2degrees(evdx)
  evdy_deg = kilometer2degrees(evdy)
  xi_new, yi_new = evlon+evdx_deg, evlat+evdy_deg
  xi_err = abs(evlon - xi_new)
  yi_err = abs(evlat - yi_new)

  x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 ); xerr.append(xi_err); yerr.append(yi_err); zerr.append(evdz)

ax2.scatter(x,z,c=z, cmap=cmap, s=s, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
#ax2.errorbar(x, z, xerr=xerr, yerr=zerr, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=19)

# add seismic events sielfeld et al 2019 and gcmt
x = []; y = []; z = []; s = []
for event in catalog_sielfeld2019:
  evlon = event[0]
  evlat = event[1]
  evdep = event[2]
  evmag = event[3]
  x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 )

evlon = catalog_gcmt[0]
evlat = catalog_gcmt[1]
evdep = catalog_gcmt[2]
evmag = 5.5 #catalog_gcmt[3]
x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 )

ax2.scatter(x,z,c=z, cmap=cmap, s=s, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
ax2.scatter(x,z,color='k', s=1.5, marker='o', zorder=21, alpha=0.8, clip_on=True, lw=0., edgecolors='k')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# add axis on right

ax3 = plt.axes([0.85, 0.3, 0.15, 0.6])
ax3.minorticks_on()
ax3.yaxis.tick_right()
ax3.xaxis.tick_top()
ax3.tick_params(axis='both', which='major', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in')
ax3.tick_params(axis='both', which='minor', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in')
ax3.tick_params(axis='x', which='major', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
ax3.tick_params(axis='x', which='minor', labelsize=6, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
ax3.spines['right'].set_visible(True)
ax3.spines['top'].set_visible(True)
ax3.spines['left'].set_visible(True)
ax3.spines['bottom'].set_visible(True)
ax3.set_xlim(zmin, zmax)
ax3.set_ylim(slat, nlat)
labels = [ r"%.1f$^{\circ}$" % (item) for item in ax3.get_yticks().tolist() ]
ax3.set_yticklabels(labels)
ax3.set_xlabel("Prof. (km)", fontsize=6)
ax3.xaxis.set_label_position('top')


# add seismic events pygema
x = []; y = []; z = []; s = []; xerr = []; yerr = []; zerr = []
for event in events_list:
  evlon = event[1]
  evlat = event[2]
  evdep = event[3]
  evmag = event[8]
  evdx = event[9]
  evdy = event[10]
  evdz = event[11]

  evdx_deg = kilometer2degrees(evdx)
  evdy_deg = kilometer2degrees(evdy)
  xi_new, yi_new = evlon+evdx_deg, evlat+evdy_deg
  xi_err = abs(evlon - xi_new)
  yi_err = abs(evlat - yi_new)

  x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 ); xerr.append(xi_err); yerr.append(yi_err); zerr.append(evdz)

ax3.scatter(z,y,c=z, cmap=cmap, s=s, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
#ax3.errorbar(z, y, xerr=zerr, yerr=yerr, elinewidth=0.4, ecolor='0.4', capsize=1.7, capthick=0.4, linewidth=0, zorder=19)

# add seismic events sielfeld et al 2019 and gcmt
x = []; y = []; z = []; s = []
for event in catalog_sielfeld2019:
  evlon = event[0]
  evlat = event[1]
  evdep = event[2]
  evmag = event[3]
  x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 )

evlon = catalog_gcmt[0]
evlat = catalog_gcmt[1]
evdep = catalog_gcmt[2]
evmag = 5.5 #catalog_gcmt[3]
x.append(evlon); y.append(evlat); z.append(evdep); s.append( evmag**3*1.5 )

ax3.scatter(z,y,c=z, cmap=cmap, s=s, marker='o', norm=norm, zorder=20, alpha=0.8, clip_on=True, lw=0.3, edgecolors='k')
ax3.scatter(z,y,color='k', s=1.5, marker='o', zorder=21, alpha=0.8, clip_on=True, lw=0., edgecolors='k')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# add outline legend
ml2 = plt.scatter([],[], marker='o',color='w', s=2.0**3*1.5, zorder=2000+2, alpha=1, clip_on=True, lw=1, edgecolors='k')
ml3 = plt.scatter([],[], marker='o',color='w', s=3.0**3*1.5, zorder=2000+2, alpha=1, clip_on=True, lw=1, edgecolors='k')
ml4 = plt.scatter([],[], marker='o',color='w', s=4.0**3*1.5, zorder=2000+2, alpha=1, clip_on=True, lw=1, edgecolors='k')
ml5 = plt.scatter([],[], marker='o',color='w', s=5.0**3*1.5, zorder=2000+2, alpha=1, clip_on=True, lw=1, edgecolors='k')

leg = plt.legend([ml5, ml4, ml3, ml2, s1, s2, f1, f2], 
                 ['5.0', '4.0', '3.0', '2.0' ], title=r"Mag. (M$_l$)", title_fontsize=6, labelspacing=1.2, 
                 fontsize=6, ncol=1, frameon=True, fancybox=True, shadow=False, framealpha=0.2, bbox_to_anchor=(0.55, -0.03) ) 
leg.set_zorder(1000)
frame = leg.get_frame()
frame.set_facecolor('w')

# add colorbar tomogrpahic inversion
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
cax = fig.add_axes([0.96, 0.108, 0.017, 0.14])
cb1 = plt.colorbar(sm, format='%i', extend='neither', norm=norm, spacing='proportional', orientation='vertical', cax=cax, ticks = np.linspace(zmin, zmax, 5)) 
cb1.ax.invert_yaxis()
cb1.set_label('Prof. (km)', size=6)
for j in cb1.ax.get_yticklabels(): j.set_fontsize(6)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# export figure

figname = "figs/mapa_informe.jpg" 
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close("all")



















