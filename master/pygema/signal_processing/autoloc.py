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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.core import UTCDateTime, read, Stream, Trace
from obspy.signal.filter import envelope, bandpass, lowpass, highpass
from obspy.signal.trigger import coincidence_trigger, recursive_sta_lta, trigger_onset, plot_trigger
from obspy.signal.invsim import estimate_magnitude
from obspy.geodetics.base import calc_vincenty_inverse
from obspy.taup import TauPyModel
from obspy import read_inventory
from obspy.io.xseed import Parser
from obspy.imaging.cm import pqlx
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
import scipy.ndimage as ndimage
from pprint import pprint
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io, subprocess

from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.read.parameters import find_pygema_parent_directory



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def run_autoloc_binder(stations_list, triggers_list):

  all_triggers = []
  for station, triggers in zip(stations_list, triggers_list):
    for trigger in triggers:
      all_triggers.append( [station, trigger[0].timestamp] )

  all_triggers = np.array(all_triggers)
  all_triggers = all_triggers[np.argsort(all_triggers[:, 1])]

  picker_input = []
  for trigger in all_triggers:
    pattern = "%.2f %s Z P 0 _" % ( float(trigger[1]), trigger[0] )
    picker_input.append(pattern)

  localdir = os.getcwd()
  pygema_path = find_pygema_parent_directory()
  picker_input_filename = "%s/pygema/core/BINDER_NOSC/picker_input.txt" % (pygema_path)

  np.savetxt(picker_input_filename, np.array(picker_input), fmt='%s') #delimiter=','
  if os.path.isfile("%s/pygema/core/BINDER_NOSC/events_binder_output.txt" %(pygema_path) ):
    subprocess.call("rm %s/pygema/core/BINDER_NOSC/events_binder_output.txt" %(pygema_path) , shell=True) 

  if os.path.isfile("%s/pygema/core/BINDER_NOSC/events.txt" % (pygema_path) ):
    subprocess.call("rm %s/pygema/core/BINDER_NOSC/events.txt" % (pygema_path) , shell=True)

  os.chdir("%s/pygema/core/BINDER_NOSC" %(pygema_path) )

  cmd1 = "./binder_nosc_AR picker_input.txt param.txt"
  subprocess.call(cmd1, shell=True)

  if os.path.isfile("events_binder_output.txt") and sum(1 for line in open('events_binder_output.txt')) > 1:
  #if not 'No events found' in open('events_binder_output.txt').read():
  #if os.path.isfile("events_binder_output.txt") and os.path.getsize("events_binder_output.txt") > 0:
    cmd2 = "more events_binder_output.txt | awk '{if ( NF > 10 ) print $0}' > events.txt"
    subprocess.call(cmd2, shell=True)

  os.chdir(localdir)
  event_file = "%s/pygema/core/BINDER_NOSC/events.txt" % (pygema_path) 
  if os.path.isfile(event_file) and sum(1 for line in open(event_file)) >= 1:
    events = np.loadtxt(event_file)
    events_list = []
    if events.ndim == 1:
      event = events
      event_time = UTCDateTime("%04i-%02i-%02iT%02i:%02i:%02.4f" % ( int(event[1]), int(event[2]), int(event[3]), int(event[4]), int(event[5]), float(event[6]) ) ) 
      event_lon = float(event[8])
      event_lat = float(event[7])
      event_dep = -float(event[9])
      event_nstats = int(event[10])
      event_gap = float(event[11])
      event_rms = float(event[12])
      events_list.append([event_time, event_lon, event_lat, event_dep, event_nstats, event_gap, event_rms])

    elif events.ndim > 1:
      for event in events:
        event_time = UTCDateTime("%04i-%02i-%02iT%02i:%02i:%02.4f" % ( int(event[1]), int(event[2]), int(event[3]), int(event[4]), int(event[5]), float(event[6]) ) ) 
        event_lon = float(event[8])
        event_lat = float(event[7])
        event_dep = -float(event[9])
        event_nstats = int(event[10])
        event_gap = float(event[11])
        event_rms = float(event[12])
        events_list.append([event_time, event_lon, event_lat, event_dep, event_nstats, event_gap, event_rms])

  else:
    events_list = []


  return events_list


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def get_local_magnitude(evtime, evlon, evlat, evdep, networks, stations, stlons, stlats, stalts, freqmin=2., freqmax=10.):

  starttime = evtime - 60
  endtime = evtime + 2*60
  st, gaps = get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False)
  allmags = []
  for station in stations:
    st1 = st.copy().select(station=station)
    for tr in st1:
      try:
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.015,type='hann')
        tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2)

        pygema_path = find_pygema_parent_directory()
        dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (pygema_path, tr.stats.network, tr.stats.station) )[0]
        parser = Parser(dataless)
        paz = parser.get_paz(tr.id)

        amp_min = tr.data.min()
        amp_max = tr.data.max()
        ind1 = np.where(amp_min == tr.data)[0][0]
        ind2 = np.where(amp_max == tr.data)[0][0]
        amplitude = abs(amp_max) + abs(amp_min)
        timespan = tr.times("utcdatetime")[ind2] - tr.times("utcdatetime")[ind1] 

        ind = np.where(station==stations)[0][0]
        h_dist = calc_vincenty_inverse(evlat, evlon, stlats[ind], stlons[ind])[0]/1000.

        mag = estimate_magnitude(paz, amplitude, timespan, h_dist)
        allmags.append(mag)
      except:
        continue

  if len(allmags)>0:
    evmag = np.nanmean(allmags)
  else:
    evmag = 0.0

  return evmag


