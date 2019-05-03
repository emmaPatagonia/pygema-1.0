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
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io, subprocess
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.core.kml_tools import kmlexport_stations, kmlexport_seismic_events
from pygema.core.check_myip_connection import retrieve_mysqlDB_credentials

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

max_gap = 360
max_depth = 250
outdir = "pygema/web/PyGema_Web/PyGema_Web/static" 


deadtime = 5
while True:
  utc_now = UTCDateTime().now()
  try:

    kmlexport_stations(kmlfilename='stations.kml', outdir=outdir)

    # last day events
    starttime = utc_now-86400.
    endtime = utc_now
    kmlexport_seismic_events(starttime, endtime, color='red', status='manual', 
                                                 kmlfilename='seismic_events_manual_last_day.kml', 
                                                 max_gap=max_gap, max_depth=max_depth, outdir=outdir)

    # last month events
    starttime = utc_now-86400.*30
    endtime = utc_now-86400.
    kmlexport_seismic_events(starttime, endtime, color='yellow', status='manual', 
                                                 kmlfilename='seismic_events_manual_last_month.kml', 
                                                 max_gap=max_gap, max_depth=max_depth, outdir=outdir)

    # older events
    starttime = utc_now-86400.*365*10
    endtime = utc_now-86400.*30
    kmlexport_seismic_events(starttime, endtime, color='green', status='manual', 
                                                 kmlfilename='seismic_events_manual_older.kml', 
                                                 max_gap=max_gap, max_depth=max_depth, outdir=outdir)

    # last day automatic events
    starttime = utc_now-86400
    endtime = utc_now
    kmlexport_seismic_events(starttime, endtime, color='red', status='automatic', 
                                                 kmlfilename='seismic_events_automatic_last_day.kml', 
                                                 max_gap=max_gap, max_depth=max_depth, outdir=outdir)



  except:
    pass

  for kmlfile in sorted( glob.glob(outdir+"/*kml") ):
    try:
      cmd = "cp %s %s" % (kmlfile, outdir)
      subprocess.call(cmd, shell=True)
      cmd = "scp %s user@address:/path/to/kml/files" % (kmlfile)
      subprocess.call(cmd, shell=True)
    except:
      continue

  time.sleep(deadtime)



