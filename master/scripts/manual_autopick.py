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

from pygema.read.parameters import load_station_metadata, load_stalta_parameters
from pygema.core.mysqlDB import insert_triggers_stalta
from pygema.signal_processing.autopick import get_triggers_stalta


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

networks, stations, stlons, stlats, stalts = load_station_metadata()
len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, stream_length = load_stalta_parameters()
overlap_stalta = 0.2

this_utc = UTCDateTime(2016,1,1)
last_utc = UTCDateTime(2017,1,1)
#last_utc = UTCDateTime(2019,1,1)

while this_utc <= last_utc:
  starttime = this_utc
  endtime = this_utc + stream_length 

  try:
    print("\n Current time: %s" % (UTCDateTime()))
    stations_list, triggers_list = get_triggers_stalta(networks, stations, starttime, endtime, len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, remove_traces_with_gaps=True)
    if len(triggers_list)>0:
      insert_triggers_stalta(stations_list, triggers_list)

      for station,triggers in zip(stations_list, triggers_list):
        for trigger in triggers:
          print("[%s] trig_on = %s   trig_off = %s" % (station, trigger[0], trigger[1]) )

  except:
    pass

  this_utc += stream_length*(1-overlap_stalta)




