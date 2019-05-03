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
from pygema.signal_processing.deconvolve import deconvolve_trace

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

def build_event_directory_for_nonlinloc(this_event, networks, stations, freqmin=None, freqmax=None, deconvolve=True, time_before=60, time_after=300):
  starttime = this_event[0] - time_before
  endtime = this_event[0] + time_after

  st, gaps = get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=False, merge_method=1, remove_traces_with_gaps=False)
  if len(st)>0:

    if freqmin is not None and freqmax is not None:
      for tr in st:
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.015,type='hann')
        tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3)

    if deconvolve:
      for tr in st:
        tr = deconvolve_trace(tr)

    outdir = "%s" % (this_event[0].strftime("%Y-%m-%d-%H%M%S") )
    if not os.path.isdir(outdir):
      os.makedirs(outdir) 

    pygema_path = find_pygema_parent_directory()
    if os.path.isdir(  "%s/waveforms/%s" % (pygema_path, outdir)  ):
      cmd = "cp -v %s/waveforms/%s/*pick %s" % (pygema_path, outdir, outdir)
      subprocess.call(cmd, shell=True)

    for tr in st:
      outfile = "%s/%s.%s.%s.%s.%s-%s.sac" % (outdir, tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, tr.stats.starttime.strftime("%Y%m%d%H%M%S"),tr.stats.endtime.strftime("%Y%m%d%H%M%S") )
      tr.write(outfile, "SAC")

  return outdir


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def read_nonlinloc_output(locfile):
  evtime, evlon, evlat, evdepth, evrms, everrx, everry, everrz = np.loadtxt(locfile, dtype='str')

  return UTCDateTime(evtime), float(evlon), float(evlat), float(evdepth), float(evrms), float(everrx), float(everry), float(everrz)




