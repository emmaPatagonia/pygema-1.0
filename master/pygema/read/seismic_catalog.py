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

from obspy.clients.fdsn import Client
from obspy.core import read, UTCDateTime, Stream
import numpy as np

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_iris_events(starttime, endtime, minlatitude=None, maxlatitude=None, minlongitude=None, maxlongitude=None, mindepth=None, maxdepth=None, minmagnitude=None, maxmagnitude=None):
  client = Client("IRIS")
  inv = client.get_events(starttime=starttime, endtime=endtime, minlatitude=minlatitude, maxlatitude=maxlatitude, minlongitude=minlongitude, maxlongitude=maxlongitude, mindepth=mindepth, maxdepth=maxdepth, minmagnitude=minmagnitude, maxmagnitude=maxmagnitude, limit=10000, orderby='time-asc')

  evtimes = []; evlons = []; evlats = []; evdepths = []; evmags = []
  for event in inv:
    evtime = event.origins[0].time
    evlon = event.origins[0].longitude
    evlat = event.origins[0].latitude
    evdep = event.origins[0].depth
    evmag = event.magnitudes[0].mag
    evtimes.append(evtime); evlons.append(evlon); evlats.append(evlat); evdepths.append(evdep); evmags.append(evmag)

  return evtimes, evlons, evlats, evdepths, evmags

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def read_iris_events(starttime, endtime, event_file='pygema/src/iris_events.txt', minlatitude=None, maxlatitude=None, minlongitude=None, maxlongitude=None, mindepth=None, maxdepth=None, minmagnitude=None, maxmagnitude=None):
  times = np.loadtxt(event_file, dtype='str', usecols=(0,) )
  origin_times = np.array([ UTCDateTime(time) for time in times ])
  lons, lats, depths, mags = np.loadtxt(event_file, usecols=(2,1,3,4) ).T

  inds = np.where( (origin_times>=starttime) & (origin_times<=endtime) &   (lats>=minlatitude) & (lats<=maxlatitude)  &   (lons>=minlongitude) & (lons<=maxlongitude)  &   (depths>=mindepth) & (depths<=maxdepth)  &   (mags<=maxmagnitude) & (mags>=minmagnitude)  )[0]

  if len(inds)>0:
    evtimes = origin_times[inds]
    evlons = lons[inds]
    evlats = lats[inds]
    evdeps = depths[inds]
    evmags = mags[inds]
  else:
    evtimes = []
    evlons = []
    evlats = []
    evdeps = []
    evmags = []

  return evtimes, evlons, evlats, evdeps, evmags

