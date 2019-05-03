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
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.plot.map import plot_map

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

deadtime = 10
while True:
  try:
    elon = -70.8
    wlon = -71.9
    nlat = -37.65
    slat = -38.55
    zmin = 0
    zmax = 25
    plot_map(elon, wlon, slat, nlat, zmin, zmax, add_topo=True, dark_background=True, show_plot=False, save_plot=True, savedir="pygema/web/PyGema_Web/PyGema_Web/static/figs_html", format='jpg', dpi=300)

  except:
    pass

  time.sleep(deadtime)



