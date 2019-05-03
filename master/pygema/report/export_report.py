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
from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import Rectangle, Polygon, PathPatch, Ellipse
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.geodetics.base import kilometer2degrees
from obspy.core import UTCDateTime, read, Stream
from obspy.imaging.beachball import beach, MomentTensor, mt2axes, plot_mt, mt2plane
from obspy.geodetics.base import calc_vincenty_inverse
from scipy.io.netcdf import NetCDFFile
from matplotlib.colors import LinearSegmentedColormap
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

#from pygema.report.build_latex import build_report_automatic_event
from pygema.read.parameters import load_station_metadata
from pygema.read.parameters import find_pygema_parent_directory
from pygema.read.seiscomp3 import get_streams_seiscomp3



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def plot_map(evlon, evlat, evdep, dlon=0.2, dlat=0.2, add_topo=False, show_plot=True, save_plot=False, savedir="pygema/report/figs", format='jpg', dpi=150):
  pygema_path = find_pygema_parent_directory()
  utc_now = UTCDateTime().now()


  elon = evlon+dlon
  wlon = evlon-dlon
  nlat = evlat+dlat
  slat = evlat-dlat
  zmin = 0
  zmax = 25

  fig = plt.figure()
  ax1 = plt.axes([0.1, 0.1, 0.9, 0.9])
  #ax1.set_title("last update: %s UTC" % (utc_now.strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right')

  map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='c', area_thresh=1000000, epsg=4269)
  map.drawparallels(np.arange(slat,nlat,dlat/5.), dashes=[1,2], labels=[1,1,0,0], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
  map.drawmeridians(np.arange(wlon,elon,dlon/5.), dashes=[1,2], labels=[0,0,1,1], linewidth=0.004,fontsize=6, zorder=100, labelstyle='+/-')
  #map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)

  map.drawcoastlines(linewidth=0.8,zorder=5)
  #map.drawcountries(linewidth=0.5, linestyle=':',zorder=2)
  map.fillcontinents(color='w', lake_color='steelblue',alpha=0.,zorder=1)
  if add_topo:
    #map.arcgisimage(service='World_Shaded_Relief', xpixels=1920, dpi=dpi, verbose=False, zorder=0)
    map.arcgisimage(service='NatGeo_World_Map', xpixels=1920, dpi=dpi, verbose=False, zorder=0)
    #map.arcgisimage(service='World_Topo_Map', xpixels=1920, dpi=dpi, verbose=False, zorder=0)
    #map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=1920, dpi=dpi, verbose=False, zorder=0)
    #map.arcgisimage(service='ESRI_StreetMap_World_2D', xpixels=1920, dpi=dpi, verbose=False, zorder=0)

  # add event
  xi,yi = map(evlon, evlat)
  ax1.scatter(xi,yi,color='red', s=100, marker='*', zorder=20, alpha=1, clip_on=True, lw=0.7, edgecolors='k')


  # add city
  #si,sj = map(-71.637212, -37.882078)
  #ax1.annotate('Comuna\nRalco',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
  #map.scatter(si, sj, c='k', linewidths=0.5, edgecolors='k', alpha=1, zorder=15, marker='X', s=15, clip_on=True)

  si,sj = map(-71.611683, -37.910428)
  ax1.annotate('E. Pangue',xy=(si,sj), color='C2',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
  map.scatter(si, sj, c='k', linewidths=0.6, edgecolors='k', alpha=1, zorder=15, marker='X', s=20, clip_on=True)

  si,sj = map(-71.475571, -38.046040)
  ax1.annotate('E. Ralco',xy=(si,sj), color='k',zorder=15,fontsize=6,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic', clip_on=True)
  map.scatter(si, sj, c='k', linewidths=0.6, edgecolors='C2', alpha=1, zorder=15, marker='X', s=20, clip_on=True)


  # add volcanoes
  vnames = np.loadtxt(pygema_path+'/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
  xvolc,yvolc = np.loadtxt(pygema_path+'/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
  xvolc,yvolc = map(xvolc,yvolc)
  ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='r',markeredgewidth=1.1,lw=0.,ms=5.,zorder=10,alpha=1, clip_on=True)

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
      ax1.scatter(x,y,marker='s',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=3, edgecolors='blue')
      #plt.annotate(stat, (x, y), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
    elif net=='C1':
      x,y = map(float(lon),float(lat))
      ax1.scatter(x,y,marker='s',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=3, edgecolors='blue')
      #ax1.annotate(stat, (x, y), weight='bold', fontsize=4, ha='left', va='center', clip_on=True, zorder=2000+3)


  # ADD faults sielfeld et al. 2019
  shp = map.readshapefile(pygema_path+'/src/shapes/gerd/fallas_sielfeld_etal_2019_mod', 'fallas', drawbounds=False)
  types = np.unique(np.array([ info['Name']  for info, shape in zip(map.fallas_info, map.fallas) ]))
  #print("types of faults =", types)
  for info, shape in zip(map.fallas_info, map.fallas):
    if info['Name'] == 'LOFS':
      x, y = zip(*shape) 
      map.plot(x, y, marker=None, color='k', alpha=1., linestyle='-', linewidth=0.5, zorder=9, clip_on=True)
    elif info['Name'] == 'ATF':
      x, y = zip(*shape) 
      map.plot(x, y, marker=None, color='k', alpha=1., linestyle='--', linewidth=0.5, zorder=9, clip_on=True)

  """
  # add legend
  s0 = plt.scatter([], [],color='red', s=40, marker='*', zorder=20, alpha=1, clip_on=True, lw=0.5, edgecolors='k')
  s1 = plt.scatter([], [], c='None', linewidths=1.1, edgecolors='blue', alpha=1, zorder=1000, marker='s', s=15)
  s2, = plt.plot([], [], c='None', markeredgecolor='r', markeredgewidth=1.1, lw=0., alpha=1, zorder=1000, marker='^', ms=5)
  f1, = plt.plot([-100,-99], [-100,-98], c='k', linewidth=.75, alpha=1, zorder=1000, ls='-') # darkgoldenrod
  f2, = plt.plot([-100,-99], [-100,-98], c='k', linewidth=.75, alpha=1, zorder=1000, ls='--') # darkgoldenrod
  nada = plt.scatter([], [], c='None', linewidths=0, edgecolors='None', alpha=1, zorder=1000, marker='s', s=15)
  leg = plt.legend([s0, s1, s2, f1, f2], 
                   ['Evento localizado automáticamente', 'Estación sismológica', 'Volcán Holoceno', 'S.F. Liquiñe-Ofqui', 'Falla Biobio-Aluminé' ], 
                   fontsize=6, ncol=1, frameon=True, fancybox=True, shadow=False, framealpha=0.8, loc=2 ) 
  leg.set_zorder(1000)
  """

  if save_plot:
    if savedir is None:
      outdir = "figs"
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/map.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, transparent=False, bbox_inches='tight') # 
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')




def plot_waveforms(evtime, evlon, evlat, freqmin=3, freqmax=10, show_plot=False, save_plot=True, savedir="pygema/report/figs", format='jpg', dpi=300):
  starttime = evtime - 15
  endtime = evtime + 90

  networks, stations, stlons, stlats, stalts = load_station_metadata()
  st, gaps = get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)
  for tr in st:
    ind = np.where( tr.stats.station == stations )[0][0]
    tr.stats.distance = calc_vincenty_inverse( stlats[ind], stlons[ind], evlat, evlon )[0]
    tr.detrend('demean')
    tr.detrend('linear')
    #tr.taper(max_percentage=0.015,type='hann')
    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2)

  fig = plt.figure()
  st.plot(type='section', method='full', orientation='vertical', time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig )
  ax = fig.axes[0]
  fig.suptitle("Origin Time: %s" % (evtime.strftime("%Y-%m-%d %H:%M:%S")), fontsize=10)
  transform = blended_transform_factory(ax.transData, ax.transAxes)
  for tr in st:
    ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, fontsize=6, rotation=45, va="bottom", ha="center", transform=transform, zorder=10)

  if save_plot:
    if savedir is None:
      outdir = "figs"
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/waveforms.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, transparent=False, bbox_inches='tight') # 
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')

