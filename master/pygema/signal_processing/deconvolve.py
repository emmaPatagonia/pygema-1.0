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

def deconvolve_trace(tr):
  pygema_path = find_pygema_parent_directory()
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (pygema_path, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.015,type='hann')
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  return tr


