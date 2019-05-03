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

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.core import UTCDateTime, read, Stream 
from obspy.imaging.cm import pqlx
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.read.parameters import find_pygema_parent_directory, load_warnings_rsam #load_warnings_ssam

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_rsam_and_ssam(st, gaps, rsam_utcs, rsam_values, ssam_utcs, ssam_values, ssam_freqs, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150, warnings_rsam=False):

  starttime = st[0].stats.starttime
  endtime = st[-1].stats.endtime
  network = st[0].stats.network
  station = st[0].stats.station
  rsam_times_matplotlib = [ date2num(utc.datetime) for utc in rsam_utcs ]
  ssam_times_matplotlib = [ date2num(utc.datetime) for utc in ssam_utcs ]


  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])


  #####################################################
  # init figure
  fig = plt.figure(dpi=dpi)
  fig.suptitle("%s" % (st[0].id), y=0.95, fontsize=10)
  gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[2,2,3], hspace=0.15) 
  ax = fig.add_subplot(gs[0])
  ax1 = fig.add_subplot(gs[1], sharex=ax)
  ax2 = fig.add_subplot(gs[2], sharex=ax1)
  date_format = DateFormatter('%H:%M:%S')

  for axx in [ax,ax1,ax2]:
    axx.minorticks_on()
    axx.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='on', right='on', direction='in')
    axx.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='on', right='on', direction='in')
    axx.tick_params(axis='x', which='major', labelsize=8, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
    axx.tick_params(axis='x', which='minor', labelsize=8, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
    axx.spines['right'].set_visible(True)
    axx.spines['top'].set_visible(False)
    axx.spines['left'].set_visible(True)
    axx.spines['bottom'].set_visible(False)
    axx.xaxis.set_major_formatter(date_format)
    axx.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
    axx.grid(axis='x', lw=0.7, ls=':', color='0.5')

  #####################################################
  # subplot 1

  ax.set_title("Last update: %s UTC" % (UTCDateTime().now().strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right' )
  ax.set_ylabel("Rawdata (counts)", fontsize=8)
  ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax.yaxis.offsetText.set_fontsize(6)
  plt.setp( ax.get_xticklabels(), visible=False)

  for tr in st:
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.015,type='hann')
    tr.filter("bandpass", freqmin=3, freqmax=10, corners=2)
    ax.plot(tr.times("matplotlib"), tr.data, lw=0.35, color="r")

  if len(gaps)>0:
    for gap in gaps:
      t1 = date2num(UTCDateTime(gap[4]).datetime)
      t2 = date2num(UTCDateTime(gap[5]).datetime)
      if dark_background:
        y1 = ax.get_ylim()[0]; y2 = ax.get_ylim()[1]
        ax.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='k', alpha=1, linewidth=0, zorder=100) )
        y1 = ax1.get_ylim()[0]; y2 = ax1.get_ylim()[1]
        ax1.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='k', alpha=1, linewidth=0, zorder=100) )
        y1 = ax2.get_ylim()[0]; y2 = ax2.get_ylim()[1]
        ax2.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='k', alpha=1, linewidth=0, zorder=100) )
      else:
        y1 = ax.get_ylim()[0]; y2 = ax.get_ylim()[1]
        ax.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='white', alpha=1, linewidth=0, zorder=100) )
        y1 = ax1.get_ylim()[0]; y2 = ax1.get_ylim()[1]
        ax1.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='white', alpha=1, linewidth=0, zorder=100) )
        y1 = ax2.get_ylim()[0]; y2 = ax2.get_ylim()[1]
        ax2.add_patch( Rectangle((t1, y1), t2-t1, y2-y1, fill=True, edgecolor=None, facecolor='white', alpha=1, linewidth=0, zorder=100) )



  #####################################################
  # subplot 2

  ax1.set_ylabel("RSAM", fontsize=8)
  ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax1.yaxis.offsetText.set_fontsize(6)
  plt.setp( ax1.get_xticklabels(), visible=False)

  #ax1.fill_between(rsam_times_matplotlib, np.nanmin(rsam_values), rsam_values, zorder=1, edgecolor="0.75", facecolor="0.3", alpha=1, linewidth=0.9)
  ax1.plot(rsam_times_matplotlib, rsam_values, color='r', alpha=1, linewidth=1)
  if warnings_rsam is True:
    rsam_lower_corner, rsam_upper_corner, rsam_yellow_warning, rsam_orange_warning = load_warnings_rsam(station)
    ax1.axhline(rsam_yellow_warning, color='yellow', ls='--', lw=0.8, zorder=100)
    ax1.axhline(rsam_orange_warning, color='orange', ls='--', lw=0.8, zorder=100)


  #####################################################
  # subplot 3

  ax2.set_ylabel("Frequency (Hz)", fontsize=8)
  ax2.set_ylim([0, 10])
  ax2.set_xlabel("UTC time (local time = UTC - 03:00)", fontsize=8 )

  cmap = pqlx
  nlevels = 32
  ssam_min =  np.nanmin(ssam_values)
  ssam_max =  np.nanmax(ssam_values)
  bounds = np.linspace(ssam_min,ssam_max,nlevels, endpoint=True)
  levels = np.linspace(ssam_min,ssam_max,nlevels, endpoint=True)
  norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
  x, y = np.meshgrid(ssam_times_matplotlib, ssam_freqs)
  im = ax2.contourf(x, y, ssam_values.T, cmap=cmap, norm=norm, levels=levels)

  cbar_ax = fig.add_axes([0.92, 0.14, 0.02, 0.22])
  cbar = fig.colorbar(im, cax=cbar_ax, extend='both', norm=norm, spacing='proportional', format='%i')
  cbar.ax.set_title('dB', fontsize = 8)
  for j in cbar.ax.get_yticklabels(): j.set_fontsize(8)

  ###################################################
  # export figure

  if save_plot:
    if savedir is None:
      outdir = "figs/%s" % (station)
    else:
      outdir = "%s/%s" % (savedir, station)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/%s_rsam_and_ssam.%s" % (outdir, station, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')

