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
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging 

from pygema.core.mysqlDB import select_rsam_warnings

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_warnings_rsam(networks, stations, starttime, endtime, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):

  if endtime - starttime > 86400.:
    date_format = DateFormatter('%Y-%m-%d\n%H:%M:%S')
  else:
    date_format = DateFormatter('%H:%M:%S')

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  fig = plt.figure(dpi=dpi)
  gs = mpl.gridspec.GridSpec(len(stations), 1, hspace=0.) 
  for i in range(len(stations)):
    network=networks[i]
    station = stations[i]

    if i == 0:
      ax1 = fig.add_subplot(gs[i])
      ax1.set_title("last update: %s UTC" % (endtime.strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right')
    else:
      ax1 = fig.add_subplot(gs[i], sharex=ax1)

    plt.setp(ax1.spines.values(), linewidth=.5)
    ax1.minorticks_on()
    ax1.set_ylabel(station, fontsize=6, rotation=0, labelpad=10)
    ax1.set_ylim([-1,1])
    ax1.xaxis.set_major_formatter(date_format)
    ax1.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
    plt.setp( ax1.get_yticklabels(), visible=False)
    if i == len(stations)-1:
      plt.setp( ax1.get_xticklabels(), visible=True)
      ax1.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='off', left='off', right='off', direction='in')
      ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='off', left='off', right='off', direction='in')
      ax1.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='off', left='off', right='off', direction='in', rotation=90)
      ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='off', left='off', right='off', direction='in', rotation=90)
      ax1.spines['right'].set_visible(False)
      ax1.spines['top'].set_visible(False)
      ax1.spines['left'].set_visible(False)
      ax1.spines['bottom'].set_visible(True)
    else:
      plt.setp( ax1.get_xticklabels(), visible=False)
      ax1.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
      ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
      ax1.tick_params(axis='x', which='major', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in', rotation=90)
      ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in', rotation=90)
      ax1.spines['right'].set_visible(False)
      ax1.spines['top'].set_visible(False)
      ax1.spines['left'].set_visible(False)
      ax1.spines['bottom'].set_visible(False)

    # add colors
    if network=="GM":
      # green polygon
      on = date2num(starttime.datetime); off = date2num(endtime.datetime)
      y1 = ax1.get_ylim()[0]*.1; y2 = ax1.get_ylim()[1]*1.0
      ax1.add_patch( Rectangle((on, y1), off-on, y2-y1, fill=True, edgecolor=None, facecolor='g', alpha=1, linewidth=0, zorder=1) )
      # color polygones
      colors, times, times_matplotlib = select_rsam_warnings(station, starttime, endtime)
      if len(colors)>0:
        for color,timeutc in zip(colors,times):
          on = date2num((timeutc-1.0).datetime)
          off = date2num((timeutc+1.0).datetime)
          if color=='YELLOW':
            ax1.add_patch( Rectangle((on, y1), off-on, y2-y1, fill=True, edgecolor=None, facecolor='y', alpha=1, linewidth=0, zorder=10) )
          elif color=='ORANGE':
            ax1.add_patch( Rectangle((on, y1), off-on, y2-y1, fill=True, edgecolor=None, facecolor='r', alpha=1, linewidth=0, zorder=15) )


  if save_plot:
    if savedir is None:
      outdir = "figs" 
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/warnings_rsam.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close()

  if show_plot:
    plt.show()
    plt.close()




