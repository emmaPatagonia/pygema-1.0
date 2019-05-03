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
from obspy.io.xseed import Parser
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
GEMA_PATH = "%s/GEMA/PyGEMA" % (os.getenv("HOME"))
sys.path.append( GEMA_PATH )

from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.core.mysqlDB import select_event_list
from pygema.plot.waveforms import adjustFigAspect

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

network = "GM"
station = "COPA"

this_day = UTCDateTime(2019,2,21,12)
starttime = this_day
endtime = this_day + 86400*3

#starttime = UTCDateTime(2018,10,16,5)
#endtime   = UTCDateTime(2018,10,16,5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)
for tr in st:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.005,type='hann')
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (GEMA_PATH, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  tr.filter("bandpass", freqmin=2, freqmax=10, corners=2)
  #tr.filter("highpass", freq=5)
  #tr.integrate(method='cumtrapz')


max_gap = 360
max_depth = 250
events_list = select_event_list( starttime, endtime, status='manual', table="LOC", max_gap=max_gap, max_depth=max_depth)
events = [ ]
for event in events_list:
  dic = {"time":event[0], "text": "Ml %.1f" % (event[-1]) }
  events.append(dic)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

fig = plt.figure(dpi=300)
#fig.suptitle("%s" % (st[0].id), y=0.95, fontsize=10)

st.plot(type='dayplot', 
        interval=30, 
        one_tick_per_line=False, 
        show_y_UTC_label=False, tick_format="%d/%b  %H:%M", # 
        #events = events,
        linewidth=0.5, method='full', starttime=starttime, endtime=endtime, fig=fig, title=None)

adjustFigAspect(fig,aspect=0.75)
ax = fig.axes[0]
ax.tick_params(axis='both', which='major', labelsize=6, bottom='on', top='off', left='on', right='on', direction='out')
ax.tick_params(axis='both', which='minor', labelsize=6, bottom='on', top='off', left='on', right='on', direction='out')
ax.tick_params(axis='x', which='major', labelsize=6, bottom='on', top='off', left='on', right='on', direction='out', rotation=0)
ax.tick_params(axis='x', which='minor', labelsize=6, bottom='on', top='off', left='on', right='on', direction='out', rotation=0)
ax.set_xlabel("Tiempo (min)", fontsize=6)

figname = "figs/%s/%s_helicorder.jpg" % (station, station)
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close()

