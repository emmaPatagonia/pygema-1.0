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
import matplotlib.pyplot as plt
from obspy.core import read, UTCDateTime, Stream
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.read.parameters import load_station_metadata
from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.plot.waveforms import plot_helicorder
from pygema.core.mysqlDB import select_event_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

networks, stations, stlons, stlats, stalts = load_station_metadata()


deadtime = 0.1
while True:
  this_day = UTCDateTime().now()
  starttime = UTCDateTime(this_day.strftime("%Y-%m-%d"))
  endtime = starttime + 86400.
  for network,station in zip(networks,stations):
    try:
      if network=="GM":
        st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)
        if len(st)>0:
          for tr in st:
            #tr.detrend('demean')
            #tr.detrend('linear')
            #tr.taper(max_percentage=0.015,type='hann')
            tr.filter("bandpass", freqmin=3, freqmax=10, corners=2)

          plot_helicorder(st, dark_background=True, show_plot=False, save_plot=True, savedir="pygema/web/PyGema_Web/PyGema_Web/static/figs_html", format='jpg', dpi=300)

    except:
      continue

  time.sleep(deadtime)



