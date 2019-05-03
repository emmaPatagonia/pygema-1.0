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
from matplotlib.patches import Rectangle
from matplotlib.dates import date2num, num2date, DateFormatter
from obspy.core import UTCDateTime, read, Stream
from obspy.io.xseed import Parser
from obspy.imaging.cm import pqlx
from obspy.geodetics.base import calc_vincenty_inverse
from obspy.signal.rotate import rotate_ne_rt
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
GEMA_PATH = "%s/GEMA/PyGEMA" % (os.getenv("HOME"))
sys.path.append( GEMA_PATH )

from pygema.read.seiscomp3 import get_streams_seiscomp3
from pygema.core.mysqlDB import select_events_manual_loc
from pygema.read.parameters import load_station_metadata


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

network = "GM"
station = "COP2"

time_before = 20
time_after = 100

max_gap = 360
max_depth = 250

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


events_list = select_events_manual_loc( UTCDateTime(1970, 1,1), UTCDateTime(), table="LOC", max_gap=max_gap, max_depth=max_depth)
events_list = np.array(events_list)
events_list = events_list[np.argsort(events_list[:, 0])]

# print events found on the screen
count = 1; count_list = []
print("\n\n")
for event in events_list:
  count_list.append(str(count))
  pattern = "[%i] %s    %.4f %.4f   %.2f km Ml=%.1f  dx=%.1f km  dy=%.1f km  dz=%.1f km " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S.%f"), event[1], event[2], event[3], event[8], event[9], event[10], event[11] )
  print(pattern)
  count += 1

print("(from PyGEMA database: Table = LOC; max_gap = %.1f deg; max_depth = %.1f km)" % (max_gap, max_depth) )

#print( "[%i] TODOS!!" % (count) )
flag = input("\n+ Type the seismic event that you want to plot: ")
while not flag in count_list:
  flag = input("+ Type the seismic event that you want to plot: ")
  if flag in count_list:
    break

this_event = events_list[int(flag)-1]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

starttime = this_event[0] - time_before
endtime = this_event[0] + time_after

st, gaps = get_streams_seiscomp3([network], [station], starttime, endtime, only_vertical_channel=False, merge_method=None, remove_traces_with_gaps=False)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

fig = plt.figure(dpi=300)
title = "evtime = %s UTC; evmag = %.1f\nevlat = %.2f deg; evlon = %.2f deg; evdepth = %.1f km" % (this_event[0].strftime("%Y/%m/%d %H:%M:%S"), this_event[8], this_event[2], this_event[1], this_event[3] )
fig.suptitle(title, y=0.96, fontsize=8)
ax1 = plt.subplot(3,4,(1,2) )
ax2 = plt.subplot(3,4,(5,6) )
ax3 = plt.subplot(3,4,(9,10) )
ax11 = plt.subplot(3,4,3 )
ax22 = plt.subplot(3,4,7 )
ax33 = plt.subplot(3,4,11 )
ax111 = plt.subplot(3,4,4 )
ax222 = plt.subplot(3,4,8 )
ax333 = plt.subplot(3,4,12 )
plt.subplots_adjust(wspace=0.5, hspace=0.5)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ADD WAVEFORM

st1 = st.copy()
for tr in st1:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.005,type='hann')
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (GEMA_PATH, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  tr.filter("bandpass", freqmin=2, freqmax=10, corners=2)
  #tr.filter("highpass", freq=5)
  #tr.integrate(method='cumtrapz')

tr_z = st1.select(channel='*HZ').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_n = st1.select(channel='*HN').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_e = st1.select(channel='*HE').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]

networks, stations, stlons, stlats, stalts = load_station_metadata()
ind = np.where(tr_z.stats.station == stations)[0][0]
evdist, evaz, evbaz =   calc_vincenty_inverse(stlats[ind], stlons[ind], this_event[2], this_event[1])
RR,TT = rotate_ne_rt(tr_n.data, tr_e.data, evbaz)

ax1.plot(tr_z.times("utcdatetime")-this_event[0], tr_z.data, lw=0.5, color="C0", label=tr_z.id )
ax2.plot(tr_z.times("utcdatetime")-this_event[0], RR, lw=0.5, color="C1", label=tr_z.id[0:11]+"R")
ax3.plot(tr_z.times("utcdatetime")-this_event[0], TT, lw=0.5, color="C2", label=tr_z.id[0:11]+"T")


#date_format = DateFormatter('%H:%M:%S')
for ax in [ax1,ax2,ax3]:
  ax.legend(fontsize=6, frameon=True, fancybox=True)
  ax.minorticks_on()
  ax.tick_params(axis='both', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='x', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.tick_params(axis='x', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.spines['right'].set_visible(True)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(True)
  ax.spines['bottom'].set_visible(False)
  ax.grid(axis='x', lw=0.5, ls=':', color='0.5')
  ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax.yaxis.offsetText.set_fontsize(6)
  #ax.xaxis.set_major_formatter(date_format)
  #ax.set_xlim([ date2num(starttime.datetime), date2num(endtime.datetime) ])
  ax.set_xlim([ -time_before, time_after ])
  ax.set_xlabel("Tiempo relativo al origen (s)", fontsize=7)
  ax.set_ylabel("Velocidad (m/s)", fontsize=7)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ADD FFT

st2 = st.copy()
for tr in st2:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.005,type='hann')
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (GEMA_PATH, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  tr.filter("bandpass", freqmin=0.02, freqmax=10, corners=2)
  #tr.filter("highpass", freq=5)
  #tr.integrate(method='cumtrapz')

tr_z = st2.select(channel='*HZ').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_n = st2.select(channel='*HN').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_e = st2.select(channel='*HE').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]

networks, stations, stlons, stlats, stalts = load_station_metadata()
ind = np.where(tr_z.stats.station == stations)[0][0]
evdist, evaz, evbaz =   calc_vincenty_inverse(stlats[ind], stlons[ind], this_event[2], this_event[1])
RR,TT = rotate_ne_rt(tr_n.data, tr_e.data, evbaz)


Fs = tr_z.stats.sampling_rate       # sampling rate
Ts = 1.0/Fs                       # sampling interval
n = tr_z.stats.npts                 # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T                         # two sides frequency range
frq = frq[range(int(n/2))]        # one side frequency range
Y = np.fft.fft(tr_z.data)           # fft computing
#Y = Y/n                           # normalization
fft_amp = abs(Y[range(int(n/2))]) 
power_spectrum_db = 10*np.log10(2*(fft_amp**2)/n)
ax11.plot(frq, fft_amp**2, lw=0.5, color="C0")

Fs = tr_z.stats.sampling_rate       # sampling rate
Ts = 1.0/Fs                       # sampling interval
n = tr_z.stats.npts                 # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T                         # two sides frequency range
frq = frq[range(int(n/2))]        # one side frequency range
Y = np.fft.fft(RR)           # fft computing
#Y = Y/n                           # normalization
fft_amp = abs(Y[range(int(n/2))]) 
power_spectrum_db = 10*np.log10(2*(fft_amp**2)/n)
ax22.plot(frq, fft_amp**2, lw=0.5, color="C1")

Fs = tr_z.stats.sampling_rate       # sampling rate
Ts = 1.0/Fs                       # sampling interval
n = tr_z.stats.npts                 # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T                         # two sides frequency range
frq = frq[range(int(n/2))]        # one side frequency range
Y = np.fft.fft(TT)           # fft computing
#Y = Y/n                           # normalization
fft_amp = abs(Y[range(int(n/2))]) 
power_spectrum_db = 10*np.log10(2*(fft_amp**2)/n)
ax33.plot(frq, fft_amp**2, lw=0.5, color="C2")


for ax in [ax11,ax22,ax33]:
  ax.set_ylabel(r"Amplitud$^{2}$", fontsize=7)
  ax.set_xlabel("Frecuencia (Hz)", fontsize=7)
  ax.minorticks_on()
  ax.tick_params(axis='both', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='x', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.tick_params(axis='x', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.spines['right'].set_visible(True)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(True)
  ax.spines['bottom'].set_visible(False)
  ax.grid(axis='x', which="both", lw=0.5, ls=':', color='0.5')
  ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
  ax.yaxis.offsetText.set_fontsize(6)
  ax.set_xlim([ 0.1,10])
  ax.set_xscale('log')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ADD PARTICLE MOTION

st3 = st.copy()
for tr in st3:
  tr.detrend('demean')
  tr.detrend('linear')
  tr.taper(max_percentage=0.005,type='hann')
  dataless = glob.glob( "%s/src/dataless/%s_%s.dataless" % (GEMA_PATH, tr.stats.network, tr.stats.station) )[0]
  parser = Parser(dataless)
  paz = parser.get_paz(tr.id)
  tr.simulate(paz_remove=paz, pre_filt=(0.01, 0.02, 50, 100), paz_simulate=None, remove_sensitivity=True)
  tr.filter("bandpass", freqmin=0.1, freqmax=0.9, corners=2)
  #tr.filter("highpass", freq=5)
  #tr.integrate(method='cumtrapz')

tr_z = st3.select(channel='*HZ').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_n = st3.select(channel='*HN').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]
tr_e = st3.select(channel='*HE').merge(method=1, interpolation_samples=-1, fill_value='interpolate')[0]

networks, stations, stlons, stlats, stalts = load_station_metadata()
ind = np.where(tr_z.stats.station == stations)[0][0]
evdist, evaz, evbaz =   calc_vincenty_inverse(stlats[ind], stlons[ind], this_event[2], this_event[1])
RR,TT = rotate_ne_rt(tr_n.data, tr_e.data, evbaz)

ax111.plot(RR, tr_z.data, lw=0.5, color="C0")
ax222.plot(TT, RR, lw=0.5, color="C1")
ax333.plot(TT, tr_z.data, lw=0.5, color="C2")


ax111.set_ylabel("Vertical", fontsize=7)
ax111.set_xlabel("Radial", fontsize=7)
ax222.set_ylabel("Radial", fontsize=7)
ax222.set_xlabel("Transversal", fontsize=7)
ax333.set_ylabel("Vertical", fontsize=7)
ax333.set_xlabel("Transversal", fontsize=7)
for ax in [ax111,ax222,ax333]:
  ax.minorticks_on()
  ax.tick_params(axis='both', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='both', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in')
  ax.tick_params(axis='x', which='major', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.tick_params(axis='x', which='minor', labelsize=7, bottom='off', top='off', left='on', right='on', direction='in', rotation=0)
  ax.spines['right'].set_visible(True)
  ax.spines['top'].set_visible(False)
  ax.spines['left'].set_visible(True)
  ax.spines['bottom'].set_visible(False)
  ax.grid(axis='x', which="both", lw=0.5, ls=':', color='0.5')
  plt.setp( ax.get_xticklabels(), visible=False)
  plt.setp( ax.get_yticklabels(), visible=False)



figname = "figs/%s/%s_waveforms.jpg" % (station, station)
plt.savefig(figname, dpi=300, bbox_inches='tight', transparent=False)
#plt.show()
plt.close()



