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
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.imaging.cm import pqlx
from obspy.core import read, UTCDateTime, Stream
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io
GEMA_PATH = "%s/GEMA/PyGEMA" % (os.getenv("HOME"))
sys.path.append( GEMA_PATH )

from pygema.read.seiscomp3 import get_streams_seiscomp3


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

starttime = UTCDateTime(2018,3,24,17,41,0)
endtime   = UTCDateTime(2018,3,24,17,53,0)

network = "GM"
station = "COPA"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# read waveform vertical component
st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False, datadir=None)


# plot
fig = plt.figure()
#fig.suptitle("%s" % (st[0].id), y=0.93, fontsize=10)
#gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[2,2,3], hspace=0.15) 
gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[1,6], hspace=0.15) 
ax = fig.add_subplot(gs[0])
#ax1 = fig.add_subplot(gs[1], sharex=ax)
#ax2 = fig.add_subplot(gs[2], sharex=ax1)
date_format = DateFormatter('%Y-%m-%d\n%H:%M:%S')

for axx in [ax]:
  axx.minorticks_on()
  axx.tick_params(axis='both', which='major', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
  axx.tick_params(axis='both', which='minor', labelsize=8, bottom='off', top='off', left='off', right='off', direction='in')
  axx.tick_params(axis='x', which='major', labelsize=7, bottom='off', top='off', left='off', right='off', direction='in', rotation=0)
  axx.tick_params(axis='x', which='minor', labelsize=7, bottom='off', top='off', left='off', right='off', direction='in', rotation=0)
  axx.spines['right'].set_visible(False)
  axx.spines['top'].set_visible(False)
  axx.spines['left'].set_visible(False)
  axx.spines['bottom'].set_visible(False)
  axx.xaxis.set_major_formatter(date_format)
  axx.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  #axx.grid(axis='x', lw=0.7, ls=':', color='0.5')

#####################################################
# subplot 1

#ax.set_ylabel("Data (cuentas)", fontsize=8)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
#ax.yaxis.offsetText.set_fontsize(6)
plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)

for tr in st:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.015,type='hann')
  tr.filter("bandpass", freqmin=1.1, freqmax=10, corners=2)
  ax.plot(tr.times("matplotlib"), tr.data, lw=0.2, color="C0")



figname = "figs/%s/%s_waveforms.jpg" % (station, station)
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=True)
#plt.show()
plt.close()





