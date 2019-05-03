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

from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.core.mysqlDB import select_triggers_stalta, select_event_list


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def adjustFigAspect(fig,aspect=1):
  xsize,ysize = fig.get_size_inches()
  minsize = min(xsize,ysize)
  xlim = .4*minsize/xsize
  ylim = .4*minsize/ysize
  if aspect < 1:
    xlim *= aspect
  else:
    ylim /= aspect

  fig.subplots_adjust(left=.5-xlim, right=.5+xlim, bottom=.5-ylim, top=.5+ylim)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_helicorder(st, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  starttime = UTCDateTime(st[0].stats.starttime.strftime("%Y-%m-%d"))
  endtime = starttime + 86400. - st[0].stats.delta

  max_gap = 360
  max_depth = 250
  events_list = select_event_list( starttime, endtime, status='manual', table="LOC", max_gap=max_gap, max_depth=max_depth)
  events = [ ]
  for event in events_list:
    dic = {"time":event[0], "text": "Ml %.1f" % (event[-1]) }
    events.append(dic)


  fig = plt.figure(dpi=dpi)
  st.plot(type='dayplot', 
        interval=60, 
        one_tick_per_line=True, 
        show_y_UTC_label=True, tick_format="%H:%M",
        events = events,
        linewidth=0.3, method='full', starttime=starttime, endtime=endtime, fig=fig, title=None)


  fig.suptitle("%s" % (st[0].id), y=0.99, fontsize=10)
  adjustFigAspect(fig,aspect=0.75)

  ax = fig.axes[0]
  ax.set_title("last update: %s UTC" % (UTCDateTime().now().strftime("%Y/%m/%d %H:%M:%S")), fontsize=6, loc='right' )
  ax.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in')
  ax.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='in', rotation=0)
  ax.spines['right'].set_visible(True)
  ax.spines['top'].set_visible(True)
  ax.spines['left'].set_visible(True)
  ax.spines['bottom'].set_visible(True)
  ax.grid(axis='x', lw=0.7, ls=':', color='0.75')

  if dark_background:
    ax.patch.set_facecolor('k')
  else:
    ax.patch.set_facecolor('w')

  if save_plot:
    if savedir is None:
      outdir = "figs/%s" % (st[0].stats.station)
    else:
      outdir = "%s/%s" % (savedir, st[0].stats.station)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/%s_helicorder.%s" % (outdir, st[0].stats.station, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_event_record_section(this_event, st, dark_background=True, show_plot=True, save_plot=False, savedir=None, format='jpg', dpi=150):

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  fig = plt.figure(dpi=dpi)
  st.plot(type='section', method='full', fig=fig, title=None, color='b', orientation='horizontal', alpha=0.7)
  fig.suptitle("Origin time: %s (M = %.1f)\nevlat: %.4f  evlon: %.4f  evdep: %.1f km" % ( this_event[0].strftime("%Y-%m-%d %H:%M:%S.%f UTC"), this_event[8], this_event[2], this_event[1], this_event[3] ), y=0.99, fontsize=10)

  ax = fig.axes[0]
  ax.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out')
  ax.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out')
  ax.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out', rotation=0)
  ax.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out', rotation=0)

  if save_plot:
    if savedir is None:
      outdir = "figs" 
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/record_section.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_three_components(network, station, starttime, endtime, freqmin=2, freqmax=10, deconvolve=False, triggers_list=None, dark_background=False, show_plot=False, save_plot=True, savedir=None, format='jpg', dpi=150):
  st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=False, merge_method=None, remove_traces_with_gaps=False)
  for tr in st:
    tr.detrend('demean')
    tr.detrend('linear')
    #tr.taper(max_percentage=0.015,type='hann')
    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2)


  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  date_format = DateFormatter('%H:%M:%S')

  fig = plt.figure(dpi=dpi)
  gs = mpl.gridspec.GridSpec(3, 1, hspace=0.) 
  ax1 = fig.add_subplot(gs[0])
  ax2 = fig.add_subplot(gs[1], sharex=ax1)
  ax3 = fig.add_subplot(gs[2], sharex=ax2)

  ax1.set_title("%s - %s UTC" % (starttime.strftime("%Y/%m/%d %H:%M:%S"), endtime.strftime("%Y/%m/%d %H:%M:%S")), fontsize=8, loc='center' )
  ax1.minorticks_on()
  ax1.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in')
  ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in')
  ax1.tick_params(axis='x', which='major', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in', rotation=0)
  ax1.tick_params(axis='x', which='minor', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in', rotation=0)
  ax1.spines['right'].set_visible(False)
  ax1.spines['top'].set_visible(False)
  ax1.spines['left'].set_visible(True)
  ax1.spines['bottom'].set_visible(False)
  ax1.set_ylabel("vertical (counts)", fontsize=10)
  #ax1.set_ylim([np.nanmin(tr.data), np.nanmax(tr.data)])
  ax1.xaxis.set_major_formatter(date_format)
  ax1.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  plt.setp( ax1.get_xticklabels(), visible=False)
  plt.setp( ax1.get_yticklabels(), visible=True)
  ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax1.yaxis.offsetText.set_fontsize(6)
  for tr in st.select(channel='*HZ'):
    ax1.plot(tr.times("matplotlib"), tr.data, lw=0.5, color="b")

  ax2.minorticks_on()
  ax2.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in')
  ax2.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in')
  ax2.tick_params(axis='x', which='major', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in', rotation=0)
  ax2.tick_params(axis='x', which='minor', labelsize=8, bottom='off', top='off', left='on', right='off', direction='in', rotation=0)
  ax2.spines['right'].set_visible(False)
  ax2.spines['top'].set_visible(False)
  ax2.spines['left'].set_visible(True)
  ax2.spines['bottom'].set_visible(False)
  ax2.set_ylabel("northing (counts)", fontsize=10)
  #ax1.set_ylim([np.nanmin(tr.data), np.nanmax(tr.data)])
  ax2.xaxis.set_major_formatter(date_format)
  ax2.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  plt.setp( ax2.get_xticklabels(), visible=False)
  plt.setp( ax2.get_yticklabels(), visible=True)
  ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax2.yaxis.offsetText.set_fontsize(6)
  for tr in st.select(channel='*HN'):
    ax2.plot(tr.times("matplotlib"), tr.data, lw=0.5, color="b")

  ax3.minorticks_on()
  ax3.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='off', left='on', right='off', direction='in')
  ax3.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='off', left='on', right='off', direction='in')
  ax3.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='off', left='on', right='off', direction='in', rotation=0)
  ax3.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='off', left='on', right='off', direction='in', rotation=0)
  ax3.spines['right'].set_visible(False)
  ax3.spines['top'].set_visible(False)
  ax3.spines['left'].set_visible(True)
  ax3.spines['bottom'].set_visible(True)
  ax3.set_ylabel("easting (counts)", fontsize=10)
  #ax3.set_ylim([np.nanmin(tr.data), np.nanmax(tr.data)])
  ax3.xaxis.set_major_formatter(date_format)
  ax3.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  plt.setp( ax3.get_xticklabels(), visible=True)
  plt.setp( ax3.get_yticklabels(), visible=True)
  ax3.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax3.yaxis.offsetText.set_fontsize(6)
  for tr in st.select(channel='*HE'):
    ax3.plot(tr.times("matplotlib"), tr.data, lw=0.5, color="b")


  if triggers_list:
    for trigger in triggers_list:
      on = date2num(trigger[0].datetime)
      off = date2num(trigger[1].datetime)
      for ax in [ax1]:#, ax2, ax3]:
        y1 = ax.get_ylim()[0]; y2 = ax.get_ylim()[1]
        ax.axvline(on, y1, y2, lw=0.7, color='r', zorder=100)
        ax.annotate( "P", (on, y2), color='r', zorder=100, fontsize=6)


  if save_plot:
    if savedir is None:
      outdir = "figs" 
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/event_three_components_%s.%s" % (outdir, station, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def plot_triggers_stalta(networks, stations, starttime, endtime, include_trig_off=False, dark_background=False, show_plot=False, save_plot=True, savedir=None, format='jpg', dpi=150):

  if dark_background:
    plt.style.use(['dark_background'])
  else:
    plt.style.use(['default'])

  fig = plt.figure(dpi=dpi)
  adjustFigAspect(fig,aspect=1.75)
  gs = mpl.gridspec.GridSpec(len(stations), 1, hspace=0.) 
  date_format = DateFormatter('%H:%M:%S')

  for i in range(len(stations)):
    network = networks[i]
    station = stations[i]
    st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)

    if i == 0:
      ax1 = fig.add_subplot(gs[i])
      ax1.set_title("last update: %s UTC" % (UTCDateTime().now().strftime("%Y/%m/%d %H:%M:%S")), fontsize=4.5, loc='right')
    else:
      ax1 = fig.add_subplot(gs[i], sharex=ax1)

    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
    ax1.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
    ax1.tick_params(axis='x', which='major', labelsize=6, bottom='off', top='off', left='off', right='off', direction='in', rotation=0)
    ax1.tick_params(axis='x', which='minor', labelsize=6, bottom='off', top='off', left='off', right='off', direction='in', rotation=0)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_ylabel(station, fontsize=6, rotation=0, labelpad=10)
    ax1.xaxis.set_major_formatter(date_format)
    ax1.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
    plt.setp( ax1.get_yticklabels(), visible=False)
    plt.setp(ax1.spines.values(), linewidth=.5)
    ax1.grid(axis='x', lw=0.7, ls=':', color='0.75')

    if i == len(stations)-1:
      plt.setp( ax1.get_xticklabels(), visible=True)
    else:
      plt.setp( ax1.get_xticklabels(), visible=False)

    tr = st.select(station=station)
    if len(tr)>0:
      tr = tr[0]
      if not isinstance(tr.data, np.ma.masked_array):
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.005,type='hann')
        tr.filter("bandpass", freqmin=1, freqmax=10, corners=2)

      if dark_background:
        ax1.plot(tr.times("matplotlib"), tr.data, lw=0.09, color="0.85")
      else:
        ax1.plot(tr.times("matplotlib"), tr.data, lw=0.09, color="k")


    triggers_list = select_triggers_stalta(station, starttime, endtime)
    for trigger in triggers_list:
      on = date2num(trigger[0].datetime)
      off = date2num(trigger[1].datetime)
      y1 = ax1.get_ylim()[0]; y2 = ax1.get_ylim()[1]
      ax1.axvline(on, y1, y2, lw=0.5, color='r', zorder=100, clip_on=True, alpha=0.9)
      #ax1.annotate( "P", (on, y2), color='r', zorder=100, fontsize=6)
      if include_trig_off:
        ax1.axvline(off, y1, y2, lw=0.7, color='g', zorder=100)
        if dark_background:
          ax1.add_patch( Rectangle((on, y1), off-on, y2-y1, fill=True, edgecolor=None, facecolor='w', alpha=0.8, linewidth=0, zorder=10) )
        else:
          ax1.add_patch( Rectangle((on, y1), off-on, y2-y1, fill=True, edgecolor=None, facecolor='w', alpha=0.8, linewidth=0, zorder=10) )


  if save_plot:
    if savedir is None:
      outdir = "figs" 
    else:
      outdir = "%s" % (savedir)

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    figname = "%s/triggers_stalta.%s" % (outdir, format)
    plt.savefig(figname, dpi=dpi, bbox_inches='tight', transparent=False)
    if not show_plot:
      plt.close('all')

  if show_plot:
    plt.show()
    plt.close('all')



