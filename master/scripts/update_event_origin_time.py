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
from obspy.core import UTCDateTime, read, Stream, Trace
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io, subprocess
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.read.parameters import load_station_metadata
from pygema.core.mysqlDB import select_event_list, update_origin_time
from pygema.read.parameters import find_pygema_parent_directory


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# load event lists from PyGEMA DB

max_gap = 355

events_list = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='automatic', table="LOC", max_gap=max_gap)
events_list_confirmed = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='confirmed', table="LOC", max_gap=max_gap)
events_list_manual = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='manual', table="LOC", max_gap=max_gap)


if len(events_list_confirmed)>0:
  for event in events_list_confirmed:
    events_list.append(event)

if len(events_list_manual)>0:
  for event in events_list_manual:
    events_list.append(event)

events_list = np.array(events_list)
events_list = events_list[np.argsort(events_list[:, 0])]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#ind = np.where(date2num(this_event[0].datetime)==np.array(datetimes_all))[0][0]
#for this_event in events_list[ind+1::]:
for this_event in events_list:
  evtime_old = this_event[0]
  evtime_new = UTCDateTime(this_event[0].strftime("%Y-%m-%dT%H:%M:%SZ"))
  print(" [updating]   %s     --->    %s" % (evtime_old, evtime_new) )
  update_origin_time(evtime_old, evtime_new, table='LOC')




