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
from pygema.core.mysqlDB import select_event_list, update_event_status, update_event_localization
from pygema.read.parameters import find_pygema_parent_directory


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# load event lists from PyGEMA DB

max_gap = 360

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

# print events found on the screen
count = 1; count_list = []
print("\n\n")
for event in events_list:
  count_list.append(str(count))
  if event[7] == 'manual':
    pattern = "[%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
  elif event[7] == 'confirmed':
    pattern = "[%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
  elif event[7] == 'automatic':
    pattern = "[%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;31;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )

  print(pattern)
  count += 1

#print( "[%i] TODOS!!" % (count) )
flag = input("\nType the seismic event that you want to remove: ")
while not flag in count_list:
  flag = input("Type the seismic event that you want to remove: ")
  if flag in count_list:
    break

idxev = int(flag)-1
this_event = events_list[idxev]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

print("")
idxs = range(int(idxev)-3, int(idxev)+3) 
for idx in idxs:
  try:
    event = events_list[idx]
    count = idx + 1
    if events_list[idx][0] == this_event[0]:
      if event[7] == 'manual':
        pattern = " --> [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
      elif event[7] == 'confirmed':
        pattern = " --> [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
      elif event[7] == 'automatic':
        pattern = " --> [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;31;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
    else:
      if event[7] == 'manual':
        pattern = "     [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
      elif event[7] == 'confirmed':
        pattern = "     [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
      elif event[7] == 'automatic':
        pattern = "     [%i] %s    %.4f %.4f   %.2f km Ml=%.1f  n=%i gap=%.1f rms=%.4f \x1b[0;31;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )

    print(pattern)
  except:
    continue

flag = input("\n+ Are you sure to remove this event? (yes/no): ")
while flag!="yes" and flag !="no":
  flag = input("+ Are you sure to remove this event? (yes/no): ")
  if flag=="yes" or flag =="no":
    break

if flag=="yes":
  update_event_status(origin_time=this_event[0], status="rejected", table="LOC")





