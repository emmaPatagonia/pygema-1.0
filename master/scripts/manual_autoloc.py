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
from obspy.core import UTCDateTime, read, Stream, Trace
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.signal_processing.autoloc import run_autoloc_binder, get_local_magnitude
from pygema.read.parameters import load_station_metadata
from pygema.core.mysqlDB import select_triggers_stalta, insert_event

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

networks, stations, stlons, stlats, stalts = load_station_metadata()

this_utc = UTCDateTime(2017,8,1)
last_utc = UTCDateTime(2018,1,1)

twin = 86400
overlap = 0.0

while this_utc <= last_utc:
  starttime = this_utc
  endtime = this_utc + twin 

  try:
    triggers_list = []
    for station in stations:
      triggers = select_triggers_stalta(station, starttime, endtime)
      triggers_list.append(triggers)

    count2 = 0; count1 = 0
    for line1, line2 in zip(stations, triggers_list):
      count2 += len(line2)
      if len(line2)>0:
        count1 += 1

    # run binder-nosc
    if count1>=3 and count2>=3:
      events_list = run_autoloc_binder(stations, triggers_list)
      if len(events_list)>0:
        print("+ We found %i events --->  %s"  % (len(events_list), this_utc.strftime("%d %b (%j) %Y") ) )
        for event in events_list:
          evtime = event[0]
          evlon = event[1]
          evlat = event[2]
          evdep = event[3]
          evnstats = event[4]
          evgap = event[5]
          evrms = event[6]
          evmag = get_local_magnitude(evtime, evlon, evlat, evdep, networks, stations, stlons, stlats, stalts, freqmin=1., freqmax=10.)

          # insert event
          insert_event(evtime, evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status="automatic", table="LOC")
          print("[insert] automatic  %s  %.4f %.4f %.1f    %.1f  %i %.1f %.1f " % (evtime, evlon, evlat, evdep, evmag, evnstats, evgap, evrms  ) )

  except:
    pass

  this_utc += twin*(1-overlap)





