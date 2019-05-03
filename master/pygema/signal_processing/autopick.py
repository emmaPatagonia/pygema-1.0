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



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_triggers_stalta(networks, stations, starttime, endtime, len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, remove_traces_with_gaps = True):

  st, gaps = get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=remove_traces_with_gaps)

  stations_list = []
  triggers_list = []
  if len(st)>0:
    #st.plot(method='full', equal_scale=False)
    for tr in st:
      try:
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.015,type='hann')
        tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=2)

        cft_rec = recursive_sta_lta(tr.data, int(len_sta*tr.stats.sampling_rate), int(len_lta*tr.stats.sampling_rate) )
        on_off = trigger_onset(cft_rec, trig_on, trig_off)
        triggers = []
        for trig in on_off:
          on = tr.times("utcdatetime")[trig[0]]
          off = tr.times("utcdatetime")[trig[1]]
          triggers.append([on, off])

        stations_list.append(tr.stats.station)
        triggers_list.append(triggers)
      except:
       continue


  return stations_list, triggers_list


