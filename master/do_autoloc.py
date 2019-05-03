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
from pygema.report.export_report import plot_map, plot_waveforms
from pygema.core.email import send_email_with_attached_files


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

networks, stations, stlons, stlats, stalts = load_station_metadata()

deadtime = 5
while True:
  utc_now = UTCDateTime().now()
  starttime = utc_now - 30*60
  endtime = utc_now 
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
        for event in events_list:
          evtime = event[0]
          evlon = event[1]
          evlat = event[2]
          evdep = event[3]
          evnstats = event[4]
          evgap = event[5]
          evrms = event[6]
          evmag = get_local_magnitude(evtime, evlon, evlat, evdep, networks, stations, stlons, stlats, stalts, freqmin=2., freqmax=10.)

          # insert event to pygema db
          insert_event(evtime, evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status="automatic", table="LOC")

          # export automatic figures
          figsdir = "%s/GEMA/PyGEMA/pygema/web/PyGema_Web/PyGema_Web/static/reports/automatic/%s" % (os.getenv("HOME"), evtime.strftime("%Y-%m-%dT%H:%M:%S+00:00") )
          if not os.path.isdir(figsdir):
            os.mkdir(figsdir)

          plot_map(evlon, evlat, evdep, dlon=0.3, dlat=0.3, add_topo=True, show_plot=False, save_plot=True, savedir=figsdir, format='jpg', dpi=200)
          plot_waveforms(evtime, evlon, evlat, freqmin=2, freqmax=10, show_plot=False, save_plot=True, savedir=figsdir, format='jpg', dpi=200)

          # send warning email
          message = "[ automatic ] \n Origin Time: %s   mag = %.1f \n evlon = %.4f deg;  evlat = %.4f deg;  evdep = %.1f km\n recorded by %i stations\n gap = %.1f deg; rms = %.1f" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evmag, evlon, evlat, evdep, evnstats, evgap, evrms)
          send_email_with_attached_files(message, figsdir=figsdir)


  except:
    pass

  time.sleep(deadtime)





