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

deadtime = 0.1
while True:
  utc_now = UTCDateTime().now()
  starttime = utc_now - stream_length 
  endtime = utc_now 

  try:
    stations_list, triggers_list = get_triggers_stalta(networks, stations, starttime, endtime, len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, remove_traces_with_gaps=True)
    if len(triggers_list)>0:
      insert_triggers_stalta(stations_list, triggers_list)

  except:
    pass

  time.sleep(deadtime)





