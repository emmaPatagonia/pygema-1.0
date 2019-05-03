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
from pygema.plot.elevation import reverse_colourmap, add_dem_srtm1_topo
from pygema.read.parameters import load_station_metadata
from pygema.read.seismic_catalog import read_iris_events,  get_iris_events
from pygema.read.parameters import find_pygema_parent_directory
from pygema.core.mysqlDB import select_event_list


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def plot_events_time_serie(evtimes, starttime, endtime, evmags=None, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):

  if (endtime - starttime > 86400.) and (endtime - starttime <= 86400.*2) :
    date_format = DateFormatter('%Y-%m-%d\n%H:%M:%S')
  elif (endtime - starttime <= 86400.):
    date_format = DateFormatter('%H:%M:%S')
  else:
    date_format = DateFormatter('%Y-%m-%d')

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  origin_times = sorted(evtimes)
  nevents = []; utcs = []
  this_day = UTCDateTime(starttime.strftime("%Y-%m-%d"))
  while this_day <= endtime:
    inds = np.where( (np.array(origin_times)>=this_day) & (np.array(origin_times)<(this_day+86400.) ) )[0]
    num = len(inds)
    utcs.append( date2num(this_day.datetime) )
    nevents.append(num)
    this_day+=86400

  fig = plt.figure(dpi=dpi)
  gs = mpl.gridspec.GridSpec(2, 1, hspace=0) 
  ax1 = fig.add_subplot(gs[0])

  ax1.set_title("last update: %s UTC" % (endtime.strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right')
  ax1.minorticks_on()
  ax1.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.spines['right'].set_visible(True)
  ax1.spines['top'].set_visible(True)
  ax1.spines['left'].set_visible(True)
  ax1.spines['bottom'].set_visible(True)
  ax1.set_ylabel("# events per day")
  ax1.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  ax1.xaxis.set_major_formatter(date_format)
  markerline, stemlines, baseline = ax1.stem(utcs, nevents, markerfmt='None', basefmt='None', linefmt='r-')
  plt.setp(stemlines, 'linewidth', 1.2)
  plt.setp( ax1.get_xticklabels(), visible=False)

  ax2 = fig.add_subplot(gs[1])
  ax2.minorticks_on()
  ax2.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax2.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax2.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=90)
  ax2.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=90)
  ax2.spines['right'].set_visible(True)
  ax2.spines['top'].set_visible(True)
  ax2.spines['left'].set_visible(True)
  ax2.spines['bottom'].set_visible(True)
  ax2.set_ylabel("accumulated events")
  ax2.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  ax2.xaxis.set_major_formatter(date_format)

  cumsum = np.cumsum( np.ones(len(origin_times)) )
  cumsum_times = []
  for time in origin_times:
    cumsum_times.append( date2num(time.datetime) )

  ax2.step(cumsum_times, cumsum, color='r', where='pre', zorder=2, lw=0.8)

  if save_plot:
    if savedir is None:
      outdir = "figs"
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/events_time_serie.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_events_map(wlon, elon, nlat, slat, event_file_iris='pygema/src/iris_events.txt', add_topo=False, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  fig = plt.figure(dpi=dpi)

  ax1 = plt.axes([0.1, 0.3, 0.6, 0.6])
  #ax1.set_title("last update: %s UTC" % (endtime.strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right')
  ax1.minorticks_on()
  ax1.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.spines['right'].set_visible(True)
  ax1.spines['top'].set_visible(True)
  ax1.spines['left'].set_visible(True)
  ax1.spines['bottom'].set_visible(True)

  map = Basemap(llcrnrlon=-75, urcrnrlon=-70, llcrnrlat=-40, urcrnrlat=-35, projection='mill', resolution='i', area_thresh=1000)
  if dark_background:
    map.drawparallels(np.arange(-90,90,1), labels=[1,0,0,0], linewidth=0.006,fontsize=6, zorder=100, color='w')
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,1), labels=[0,0,1,0], linewidth=0.006,fontsize=6, zorder=100, color='w')
    #map.drawmapscale(minlongitude+0.11, maxlatitude-0.04, minlongitude+0.11, maxlatitude-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='w', fillcolor1='w', fillcolor2='w', format='%d', zorder=1000+1)
  else:
    map.drawparallels(np.arange(-90,90,1), labels=[1,0,0,0], linewidth=0.006,fontsize=6, zorder=100)
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,1), labels=[0,0,1,0], linewidth=0.006,fontsize=6, zorder=100)
    #map.drawmapscale(minlongitude+0.11, maxlatitude-0.04, minlongitude+0.11, maxlatitude-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)

  map.drawcoastlines(linewidth=0.8,zorder=5)
  map.drawcountries(linewidth=0.5,zorder=5, linestyle='-')
  map.fillcontinents(color='0.8', lake_color='steelblue',zorder=5,alpha=0.7)

  # add seismicity
  evtimes, evlons, evlats, evdepths, evmags = read_iris_events(starttime=UTCDateTime(2000,1,1), endtime=UTCDateTime(2018,12,31), event_file=event_file_iris, minlatitude=-40, maxlatitude=-35, minlongitude=-75, maxlongitude=-70, mindepth=0, maxdepth=150, minmagnitude=3, maxmagnitude=10)
  x,y = map(evlons,evlats)
  ax1.scatter(x,y,marker='o',c=evdepths, cmap=plt.cm.jet_r, s=evmags**3/12, zorder=3000+2, alpha=1, clip_on=True, lw=0.2, edgecolors='k')

  # add stations
  networks, stations, stlons, stlats, stalts = load_station_metadata("pygema/src/stationALL.lst")
  dx = -3700; dy = -2800
  for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
    if net=='GM':
      x,y = map(float(lon),float(lat))
      plt.scatter(x,y,marker='s',color='None', s=15, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      #plt.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
    elif net=='C1':
      x,y = map(float(lon),float(lat))
      plt.scatter(x,y,marker='d',color='None', s=15, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      #plt.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)

  # add volcanoes
  vnames = np.loadtxt('pygema/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
  xvolc,yvolc = np.loadtxt('pygema/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
  xvolc,yvolc = map(xvolc,yvolc)
  ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.1,lw=0.,ms=4.5,zorder=10,alpha=1, clip_on=True)
  #for x,y,vname in zip(xvolc,yvolc,vnames):
    #if vname == 'Copahue' or vname == 'Callaqui' or vname == 'Lonquimay' or vname == 'Tolhuaca':
      #plt.annotate('V. '+vname, (x-3000, y+3000), weight='normal', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')

  #
  def plot_fault_type_luis(m, fault_type='falla observada', color_shp='m', lw_shp=1., ls_shp='-', alp_shp=0.8, zo_shp=2, marker_size=2.5, marker_scale=90, dp=1):
    shp = m.readshapefile('pygema/src/shapes/lucho/lofz', 'fallas', drawbounds=False)
    types = np.unique(np.array([ info['Tipo_de_fa']  for info, shape in zip(m.fallas_info, m.fallas) ]))
    for info, shape in zip(m.fallas_info, m.fallas):
      if info['Tipo_de_fa'] == fault_type:
        x, y = zip(*shape) 
        m.plot(x, y, marker=None, color=color_shp, alpha=alp_shp, linestyle=ls_shp, linewidth=lw_shp, zorder=zo_shp, clip_on=True)

    return m   

  plot_fault_type_luis(map, fault_type='falla de rumbo dextral', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)
  plot_fault_type_luis(map, fault_type='falla dextral inversa', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)

  x1,y1 = map(wlon,slat) 
  x2,y2 = map(wlon,nlat) 
  x3,y3 = map(elon,nlat) 
  x4,y4 = map(elon,slat)
  xy = [(x1,y1),(x2,y2),(x3,y3),(x4,y4)]
  poly = Polygon( xy, facecolor='None', edgecolor='k', alpha=1, zorder=11, fill=True, linewidth=1.7, clip_on=False, linestyle='solid')
  plt.gca().add_patch(poly)

  ax3 = plt.axes([0.65, 0.05, 0.35, 0.9])
  ax3.minorticks_on()
  ax3.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax3.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax3.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax3.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax3.spines['right'].set_visible(True)
  ax3.spines['top'].set_visible(True)
  ax3.spines['left'].set_visible(True)
  ax3.spines['bottom'].set_visible(True)


  map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='c', area_thresh=0.1)
  if dark_background:
    map.drawparallels(np.arange(-90,90,.1), labels=[0,1,0,0], linewidth=0.006,fontsize=6, zorder=100, color='w')
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,.2), labels=[0,0,1,1], linewidth=0.006,fontsize=6, zorder=100, color='w')
    map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='w', fillcolor1='w', fillcolor2='w', format='%d', zorder=1000+1)
  else:
    map.drawparallels(np.arange(-90,90,.1), labels=[0,1,0,0], linewidth=0.006,fontsize=6, zorder=100)
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,.2), labels=[0,0,1,1], linewidth=0.006,fontsize=6, zorder=100)
    map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)

  map.drawcoastlines(linewidth=0.8,zorder=5)
  map.drawcountries(linewidth=0.5,zorder=5, linestyle='-')
  map.fillcontinents(color='None', lake_color='steelblue',zorder=5,alpha=0.7)

  evtimes, evlons, evlats, evdepths, evmags = read_iris_events(starttime=UTCDateTime(2000,1,1), endtime=UTCDateTime(2018,12,31), event_file=event_file_iris, minlatitude=-40, maxlatitude=-35, minlongitude=-75, maxlongitude=-70, mindepth=0, maxdepth=150, minmagnitude=3, maxmagnitude=10)
  x,y = map(evlons,evlats)
  ax3.scatter(x,y,marker='o',c=evdepths, cmap=plt.cm.jet_r, s=evmags**3/3, zorder=3000+2, alpha=1, clip_on=True, lw=0.2, edgecolors='k')

  x,y=map(-71.40, -38.04)
  mt = MomentTensor([0.16, 2.15, -2.32, 0.17, -0.01, 1.16],24.)
  plane = mt2plane(mt)
  strike = plane.strike
  dip = plane.dip
  rake = plane.rake
  if rake >= 180: rake = rake - 360
  b1 = beach([strike,dip,rake], xy=(x,y), width=6000, linewidth=0.8, facecolor='k', alpha=.9)      
  b1.set_zorder(10000000000-1)
  b1.set_clip_on(True)
  ax3.add_collection(b1)
  ax3.annotate('2006-12-31 14:55\nMw = 5.5', (x,y), textcoords='offset points', size=5, xytext=(5, -15), zorder=10000, weight='normal', clip_on=True)


  # add stations
  networks, stations, stlons, stlats, stalts = load_station_metadata("pygema/src/stationALL.lst")
  dx = -3700; dy = -2800
  for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
    if net=='GM':
      x,y = map(float(lon),float(lat))
      plt.scatter(x,y,marker='s',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      #plt.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
    elif net=='C1':
      x,y = map(float(lon),float(lat))
      plt.scatter(x,y,marker='d',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      plt.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)

  # add volcanoes
  vnames = np.loadtxt('pygema/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
  xvolc,yvolc = np.loadtxt('pygema/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
  xvolc,yvolc = map(xvolc,yvolc)
  ax3.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.1,lw=0.,ms=7.,zorder=10,alpha=1, clip_on=True)
  for x,y,vname in zip(xvolc,yvolc,vnames):
    if vname == 'Copahue' or vname == 'Tolhuaca':
      plt.annotate('V. '+vname, (x-13000, y+3000), weight='normal', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
    if vname == 'Callaqui' or vname == 'Lonquimay':
      plt.annotate('V. '+vname, (x-3000, y+3000), weight='normal', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')



  # add topography shaded
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
    colormap = 'pygema/src/ib12.cpt' 
    cdict = gmtColormap(colormap)
    cmap1 = LinearSegmentedColormap('my_colormap',cdict,256)
    cmap1 = reverse_colourmap(cmap1)
    cmap1.set_bad('w',alpha=0.)
    intensity = hillshade(elev,scale=20,azdeg=290.0,altdeg=45.0)
    rgb = set_shade(elev,intensity=intensity,cmap=cmap1) 
    map.imshow(rgb, cmap=cmap1, zorder = 1, origin='upper',clip_on=True, extent=[left, right, bottom, top])
    colormap = 'pygema/src/DEM_print.cpt'
    cdict = gmtColormap(colormap) 
    cmap2 = LinearSegmentedColormap('my_colormap',cdict,256)
    cmap2.set_bad('w',alpha=0.)
    #map.imshow(elev, cmap=cmap2, zorder = 1,origin='upper',alpha=0.3,clip_on=True, extent=[left, right, bottom, top])
    return map

  netcdffile = "pygema/src/DEM_srtm1_30m.nc" 
  if add_topo:
    add_dem_srtm1_topo(map, netcdffile, wlon, elon, slat, nlat)

  #
  def plot_fault_type_luis(m, fault_type='falla observada', color_shp='m', lw_shp=1., ls_shp='-', alp_shp=0.8, zo_shp=2, marker_size=2.5, marker_scale=90, dp=1):
    shp = m.readshapefile('pygema/src/shapes/lucho/lofz', 'fallas', drawbounds=False)
    types = np.unique(np.array([ info['Tipo_de_fa']  for info, shape in zip(m.fallas_info, m.fallas) ]))
    for info, shape in zip(m.fallas_info, m.fallas):
      if info['Tipo_de_fa'] == fault_type:
        x, y = zip(*shape) 
        m.plot(x, y, marker=None, color=color_shp, alpha=alp_shp, linestyle=ls_shp, linewidth=lw_shp, zorder=zo_shp, clip_on=True)

    return m   

  plot_fault_type_luis(map, fault_type='falla de rumbo dextral', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)
  plot_fault_type_luis(map, fault_type='falla dextral inversa', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)


  # add city
  si,sj = map(-71.637212, -37.882078)
  ax3.annotate('Comuna\nRalco',xy=(si-2000,sj+1500),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)

  si,sj = map(-71.611683, -37.910428)
  ax3.annotate('Embalse\nPangue',xy=(si-3500,sj-5000),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)

  si,sj = map(-71.475571, -38.046040)
  ax3.annotate('Embalse\nRalco',xy=(si-3500,sj-5000),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)



  ax2 = plt.axes([0.2, 0.1, 0.4, 0.15])
  ax2.minorticks_on()
  ax2.tick_params(axis='both', which='major', labelsize=6, bottom='off', top='on', left='on', right='off', direction='in')
  ax2.tick_params(axis='both', which='minor', labelsize=6, bottom='off', top='on', left='on', right='off', direction='in')
  ax2.tick_params(axis='x', which='major', labelsize=6, bottom='off', top='on', left='on', right='off', direction='in', rotation=0)
  ax2.tick_params(axis='x', which='minor', labelsize=6, bottom='off', top='on', left='on', right='off', direction='in', rotation=0)
  ax2.spines['right'].set_visible(False)
  ax2.spines['top'].set_visible(True)
  ax2.spines['left'].set_visible(True)
  ax2.spines['bottom'].set_visible(False)
  ax2.set_xlim(-75, -70)
  ax2.set_ylim(120,0)
  ax2.xaxis.set_ticks_position('top')

  # add seismicity
  evtimes, evlons, evlats, evdepths, evmags = read_iris_events(starttime=UTCDateTime(2000,1,1), endtime=UTCDateTime(2018,12,31), event_file=event_file_iris, minlatitude=-40, maxlatitude=-35, minlongitude=-75, maxlongitude=-70, mindepth=0, maxdepth=150, minmagnitude=3, maxmagnitude=10)
  ax2.scatter(evlons,evdepths,marker='o',c=evdepths, cmap=plt.cm.jet_r, s=evmags**3/12, zorder=3000+2, alpha=1, clip_on=True, lw=0.3, edgecolors='k')

  x1,y1 = (wlon,20)
  x2,y2 = (wlon,0)
  x3,y3 = (elon,0) 
  x4,y4 = (elon,20)
  xy = [(x1,y1),(x2,y2),(x3,y3),(x4,y4)]
  poly = Polygon( xy, facecolor='None', edgecolor='k', alpha=1, zorder=11, fill=True, linewidth=1.7, clip_on=False, linestyle='solid')
  plt.gca().add_patch(poly)


  if save_plot:
    if savedir is None:
      outdir = "figs"
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/map_seismic_events.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')



def plot_fault_type_luis(m, fault_type='falla observada', color_shp='m', lw_shp=1., ls_shp='-', alp_shp=0.8, zo_shp=2, marker_size=2.5, marker_scale=90, dp=1):
  pygema_path = find_pygema_parent_directory()

  shp = m.readshapefile(pygema_path+'/src/shapes/lucho/lofz', 'fallas', drawbounds=False)
  types = np.unique(np.array([ info['Tipo_de_fa']  for info, shape in zip(m.fallas_info, m.fallas) ]))
  for info, shape in zip(m.fallas_info, m.fallas):
    if info['Tipo_de_fa'] == fault_type:
      x, y = zip(*shape) 
      m.plot(x, y, marker=None, color=color_shp, alpha=alp_shp, linestyle=ls_shp, linewidth=lw_shp, zorder=zo_shp, clip_on=True)

  return m   


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_map(wlon, elon, nlat, slat, add_topo=False, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):
  pygema_path = find_pygema_parent_directory()
  utc_now = UTCDateTime().now()

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  fig = plt.figure(dpi=dpi)
  ax1 = plt.axes([0.1, 0.1, 0.9, 0.9])
  ax1.set_title("last update: %s UTC" % (utc_now.strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right')
  ax1.minorticks_on()
  ax1.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax1.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax1.spines['right'].set_visible(True)
  ax1.spines['top'].set_visible(True)
  ax1.spines['left'].set_visible(True)
  ax1.spines['bottom'].set_visible(True)

  map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='c', area_thresh=1000)
  map = Basemap(llcrnrlon=wlon, urcrnrlon=elon, llcrnrlat=slat, urcrnrlat=nlat, projection='mill', resolution='c', area_thresh=1000)
  if dark_background:
    map.drawparallels(np.arange(-90,90,.1), labels=[1,1,0,0], linewidth=0.006,fontsize=6, zorder=100, textcolor='w')
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,.2), labels=[0,0,0,1], linewidth=0.006,fontsize=6, zorder=100, textcolor='w')
    map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='w', fillcolor1='w', fillcolor2='w', format='%d', zorder=1000+1)
  else:
    map.drawparallels(np.arange(-90,90,.1), labels=[1,1,0,0], linewidth=0.006,fontsize=6, zorder=100)
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,.2), labels=[0,0,0,1], linewidth=0.006,fontsize=6, zorder=100)
    map.drawmapscale(wlon+0.11, nlat-0.04, wlon+0.11, nlat-0.04, 20, barstyle='fancy', fontsize=6, yoffset=0.1 , fontcolor='k', fillcolor1='k', fillcolor2='w', format='%d', zorder=1000+1)

  map.drawcoastlines(linewidth=0.8,zorder=5)
  map.drawcountries(linewidth=0.5,zorder=5, linestyle='-')
  map.fillcontinents(color='None', lake_color='steelblue',zorder=5,alpha=0.7)


  # add topography shaded
  netcdffile = pygema_path+"/src/DEM_srtm1_30m.nc" 
  if add_topo:
    add_dem_srtm1_topo(map, netcdffile, wlon, elon, slat, nlat)

  # add stations
  networks, stations, stlons, stlats, stalts = load_station_metadata()
  dx = -3700; dy = -2800
  for net, stat, lon, lat in zip(networks, stations, stlons, stlats):
    if net=='GM':
      x,y = map(float(lon),float(lat))
      ax1.scatter(x,y,marker='s',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      #ax1.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)
    elif net=='C1':
      x,y = map(float(lon),float(lat))
      ax1.scatter(x,y,marker='d',color='None', s=25, zorder=2000+2, alpha=1, clip_on=True, lw=1.3, edgecolors='C0')
      ax1.annotate(stat, (x+dx, y+dy), weight='bold', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3)

  # add volcanoes
  vnames = np.loadtxt(pygema_path+'/src/smithsonian_volcanoes.txt', usecols=(0,), dtype='str')
  xvolc,yvolc = np.loadtxt(pygema_path+'/src/smithsonian_volcanoes.txt', usecols=(2,1)).T
  xvolc,yvolc = map(xvolc,yvolc)
  ax1.plot(xvolc,yvolc,marker='^',color='None',markeredgecolor='C3',markeredgewidth=1.1,lw=0.,ms=7.,zorder=10,alpha=1, clip_on=True)
  for x,y,vname in zip(xvolc,yvolc,vnames):
    if vname == 'Copahue' or vname == 'Tolhuaca':
      ax1.annotate('V. '+vname, (x-13000, y+3000), weight='normal', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')
    if vname == 'Callaqui' or vname == 'Lonquimay':
      ax1.annotate('V. '+vname, (x-3000, y+3000), weight='normal', fontsize=5, ha='left', va='center', clip_on=True, zorder=2000+3, fontstyle='italic')

  plot_fault_type_luis(map, fault_type='falla de rumbo dextral', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)
  plot_fault_type_luis(map, fault_type='falla dextral inversa', color_shp='k', lw_shp=.75, ls_shp='--', alp_shp=0.8, zo_shp=6, marker_size=2.5, marker_scale=90, dp=1)

  # add city
  si,sj = map(-71.637212, -37.882078)
  ax1.annotate('Comuna\nRalco',xy=(si-2000,sj+1500),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)

  si,sj = map(-71.611683, -37.910428)
  ax1.annotate('Embalse\nPangue',xy=(si-3500,sj-5000),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)

  si,sj = map(-71.475571, -38.046040)
  ax1.annotate('Embalse\nRalco',xy=(si-3500,sj-5000),zorder=200000000,fontsize=5,fontweight='normal',xytext=(-0,-0),textcoords='offset points', fontstyle='italic')
  map.scatter(si, sj, c='None', linewidths=0.9, edgecolors='purple', alpha=1, zorder=1000, marker='D', s=12)


  # add seismicity last day
  events_list = []

  starttime = utc_now - 86400.
  endtime = utc_now
  events = select_event_list(starttime, endtime, status="automatic", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="confirmed", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="manual", max_gap=270)
  for event in events:
    events_list.append(event)

  lons = []; lats = []; depths = []
  for event in events_list:
    lons.append(event[1])    
    lats.append(event[2])    
    depths.append(event[3])    

  x,y = map(lons,lats)
  ax1.scatter(x,y,color="r", lw=0.5, edgecolors='k', s=20, zorder=100)

  # add seismicity last week
  events_list = []

  starttime = utc_now - 86400.*7
  endtime = utc_now  - 86400
  events = select_event_list(starttime, endtime, status="automatic", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="confirmed", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="manual", max_gap=270)
  for event in events:
    events_list.append(event)

  lons = []; lats = []; depths = []
  for event in events_list:
    lons.append(event[1])    
    lats.append(event[2])    
    depths.append(event[3])    

  x,y = map(lons,lats)
  ax1.scatter(x,y,color="y", lw=0.5, edgecolors='k', s=20, zorder=100)


  # add seismicity last month
  events_list = []

  starttime = utc_now - 86400.*28
  endtime = utc_now  - 86400*7
  events = select_event_list(starttime, endtime, status="automatic", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="confirmed", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="manual", max_gap=270)
  for event in events:
    events_list.append(event)

  lons = []; lats = []; depths = []
  for event in events_list:
    lons.append(event[1])    
    lats.append(event[2])    
    depths.append(event[3])    

  x,y = map(lons,lats)
  ax1.scatter(x,y,color="g", lw=0.5, edgecolors='k', s=20, zorder=100)


  # add seismicity older
  events_list = []

  starttime = utc_now - 86400.*365*10
  endtime = utc_now  - 86400*28
  events = select_event_list(starttime, endtime, status="automatic", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="confirmed", max_gap=270)
  for event in events:
    events_list.append(event)

  events = select_event_list(starttime, endtime, status="manual", max_gap=270)
  for event in events:
    events_list.append(event)

  lons = []; lats = []; depths = []
  for event in events_list:
    lons.append(event[1])    
    lats.append(event[2])    
    depths.append(event[3])    

  x,y = map(lons,lats)
  ax1.scatter(x,y,color="b", lw=0.5, edgecolors='k', s=20, zorder=100)

  # add legend
  s3 = plt.scatter([], [], c='None', linewidths=1.1, edgecolors='C0', alpha=1, zorder=1000, marker='s', s=15)
  s4, = plt.plot([], [], c='None', markeredgecolor='C3', markeredgewidth=1.1, lw=0., alpha=1, zorder=1000, marker='^', ms=5)
  f1, = plt.plot([-100,-99], [-100,-98], c='0.7', linewidth=.75, alpha=1, zorder=1000, ls='--') # darkgoldenrod
  #s5 = plt.scatter([],[], marker='*',color='None', s=60, zorder=2000+2, alpha=1, clip_on=True, lw=.6, edgecolors='0.7')
  s51 = plt.scatter([],[], marker='o',color='r', s=40, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
  s52 = plt.scatter([],[], marker='o',color='y', s=40, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
  s53 = plt.scatter([],[], marker='o',color='g', s=40, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
  s54 = plt.scatter([],[], marker='o',color='b', s=40, zorder=2000+2, alpha=1, clip_on=True, lw=.5, edgecolors='k')
  s2 = plt.scatter([],[], marker='D',color='None', s=10, zorder=2000+2, alpha=1, clip_on=True, lw=.7, edgecolors='purple')


  leg = plt.legend([s51, s52, s53, s54, s3, s4, f1], ['Past day event', 'Past week event', 'Past month event', 'Older event', 'Seismic station', 'Present-day volcano', r'Liqui$\mathrm{\~n}$e-Ofqui F-Z' ], fontsize=5, frameon=True, fancybox=True, shadow=True, framealpha=0.8, loc=4) 
  leg.set_zorder(1000)


  if save_plot:
    if savedir is None:
      outdir = "figs"
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/map_seismic_events.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')





