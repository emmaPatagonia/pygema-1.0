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
from obspy.geodetics.base import calc_vincenty_inverse
from obspy.io.xseed import Parser
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
GEMA_PATH = "%s/GEMA/PyGEMA" % (os.getenv("HOME"))
sys.path.append( GEMA_PATH )

from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.core.mysqlDB import select_event_list, select_triggers_stalta
from pygema.plot.waveforms import adjustFigAspect
from pygema.read.parameters import load_station_metadata

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

time_before = 0
time_after = 120


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

max_gap = 270
max_depth = 50

events_list = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='automatic', table="LOC", max_gap=max_gap, max_depth=max_depth)
events_list_confirmed = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='confirmed', table="LOC", max_gap=max_gap, max_depth=max_depth)
events_list_manual = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='manual', table="LOC", max_gap=max_gap, max_depth=max_depth)

if len(events_list_confirmed)>0:
  for event in events_list_confirmed:
    events_list.append(event)

if len(events_list_manual)>0:
  for event in events_list_manual:
    events_list.append(event)

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

networks, stations, stlons, stlats, stalts = load_station_metadata()

starttime = this_event[0] - time_before
endtime = this_event[0] + time_after

st, gaps = get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)
#st = st.select(network="GM")
evdists = []
for tr in st:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.005,type='hann')
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (GEMA_PATH, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  tr.filter("bandpass", freqmin=1, freqmax=10, corners=2)
  #tr.filter("highpass", freq=5)
  #tr.integrate(method='cumtrapz')
  ind = np.where(tr.stats.station == stations)[0][0]
  evdist = calc_vincenty_inverse(stlats[ind], stlons[ind], this_event[2], this_event[1])[0] / 1000.
  tr.stats.distance = evdist
  evdists.append(evdist)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

fig = plt.figure()
adjustFigAspect(fig,aspect=1.75)
date_format = DateFormatter('%H:%M:%S')

ax = plt.subplot(111)
title = "evtime = %s UTC; evmag = %1.f\nevlat = %.2f deg; evlon = %.2f deg; evdepth = %.1f km" % (this_event[0].strftime("%Y/%m/%d %H:%M:%S"), this_event[8], this_event[2], this_event[1], this_event[3] )
ax.set_title(title, fontsize=8)
ax.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
ax.tick_params(axis='both', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
ax.tick_params(axis='x', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
ax.tick_params(axis='x', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
ax.spines['right'].set_visible(True)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(False)
#ax.xaxis.set_major_formatter(date_format)
#ax.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
ax.set_xlim([ time_before, time_after ])
ax.set_ylim([ 0 , 150 ])
ax.grid(axis='x', lw=0.7, ls=':', color='0.5')
ax.set_ylabel("Distancia al hipocentro (km)", fontsize=7)
ax.set_xlabel("Tiempo relativo al origen (s)", fontsize=7)

yscale = 10
for tr in st.select():
  ypos = tr.stats.distance
  waveform = tr.data/np.nanmax( abs(tr.data) ) * yscale
  ax.plot(tr.times("utcdatetime")-this_event[0], waveform+ypos, lw=0.35, color="C0", clip_on=True )
  ax.annotate(tr.stats.station+"  ", (time_after, ypos), fontsize=4, ha='right', va='bottom', fontweight='bold', clip_on=True)

  triggers_list = select_triggers_stalta(tr.stats.station, starttime, endtime)
  for trigger in triggers_list:
    on = trigger[0] - this_event[0]
    off = trigger[1] - this_event[0]
    y1 = ypos / (ax.get_ylim()[1] - ax.get_ylim()[0]) - 0.05
    y2 = ypos / (ax.get_ylim()[1] - ax.get_ylim()[0]) + 0.05
    ax.axvline(on, ymin=y1, ymax=y2, lw=0.75, color='C1', zorder=100, clip_on=True, alpha=1) #, y1, y2

figname = "figs/record_section.jpg" 
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close("all")



"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

fig = plt.figure()
st.plot(type='section', method='full', fig=fig, title=None, color='b', orientation='horizontal', alpha=0.7)
title = "evtime = %s UTC; evmag = %1.f\nevlat = %.2f deg; evlon = %.2f deg; evdepth = %.1f km" % (this_event[0].strftime("%Y/%m/%d %H:%M:%S"), this_event[8], this_event[2], this_event[1], this_event[3] )
fig.suptitle(title, y=0.99, fontsize=10)

ax = fig.axes[0]
ax.tick_params(axis='both', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out')
ax.tick_params(axis='both', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out')
ax.tick_params(axis='x', which='major', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out', rotation=0)
ax.tick_params(axis='x', which='minor', labelsize=8, bottom='on', top='on', left='on', right='on', direction='out', rotation=0)

figname = "figs/record_section.jpg" 
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close()
"""
