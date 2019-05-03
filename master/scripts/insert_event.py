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

from pygema.core.mysqlDB import insert_event


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# define hypocenter

evtime = UTCDateTime(2019,4,22,18,32,49)

evlon = -71.3148
evlat = -37.7875
evdep = 144

evmag = 2.2

evgap    = 285
evnstats = 3
evrms    = 0.8


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

print("[insert] automatic  %s  %.4f %.4f %.1f    %.1f  %i %.1f %.1f " % (evtime, evlon, evlat, evdep, evmag, evnstats, evgap, evrms  ) )
insert_event(evtime, evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status="automatic", table="LOC")


