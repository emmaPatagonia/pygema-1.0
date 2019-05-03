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
from matplotlib.transforms import blended_transform_factory
from obspy.core import UTCDateTime, read, Stream, Trace
from obspy.geodetics.base import calc_vincenty_inverse
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io, subprocess
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

from pygema.read.parameters import load_station_metadata
from pygema.core.mysqlDB import select_event_list, update_event_status, update_event_localization
from pygema.read.parameters import find_pygema_parent_directory
from pygema.signal_processing.autoloc import get_local_magnitude
from pygema.signal_processing.nonlinloc import build_event_directory_for_nonlinloc, read_nonlinloc_output
from pygema.report.export_report import plot_map, plot_waveforms
from pygema.core.email import send_email_with_attached_files


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PRINT SEISMIC EVENTS DATABASE

max_gap = 360
max_depth = 250

events_list = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='automatic', table="LOC", max_gap=max_gap, max_depth=max_depth)
events_list_confirmed = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='confirmed', table="LOC", max_gap=max_gap, max_depth=max_depth)
events_list_manual = select_event_list( UTCDateTime(1970, 1,1), UTCDateTime().now(), status='manual', table="LOC", max_gap=max_gap, max_depth=max_depth)

if len(events_list_confirmed)>0:
  for event in events_list_confirmed:
    events_list.append(event)

if len(events_list_manual)>0:
  for event in events_list_manual:
    events_list.append(event)

events_list = np.array(events_list)
events_list = events_list[np.argsort(events_list[:, 0])]

# print events found on the screen
count = 1; count_list = []
print("\n\n")
for event in events_list:
  count_list.append(str(count))
  if event[7] == 'manual':
    pattern = "[%i] %s    %.4f %.4f   %.2f km  Ml=%.1f    n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
  elif event[7] == 'confirmed':
    pattern = "[%i] %s    %.4f %.4f   %.2f km  Ml=%.1f    n=%i gap=%.1f rms=%.4f \x1b[0;32;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )
  elif event[7] == 'automatic':
    pattern = "[%i] %s    %.4f %.4f   %.2f km  Ml=%.1f    n=%i gap=%.1f rms=%.4f \x1b[0;31;40m %s \x1b[0m " % (count, event[0].strftime("%Y-%m-%d %H:%M:%S"), event[1], event[2], event[3], event[8], event[4], event[5], event[6], event[7] )

  print(pattern)
  count += 1

print("(from PyGEMA database: Table = LOC; max_gap = %.1f deg; max_depth = %.1f km)" % (max_gap, max_depth) )

#print( "[%i] TODOS!!" % (count) )
flag = input("\n+ Type the seismic event that you want to localize: ")
while not flag in count_list:
  flag = input("+ Type the seismic event that you want to localize: ")
  if flag in count_list:
    break

this_event = events_list[int(flag)-1]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# BUILD DATA DIRECTORY

pygema_path = find_pygema_parent_directory()
networks, stations, stlons, stlats, stalts = load_station_metadata()
outdir = build_event_directory_for_nonlinloc(this_event, networks, stations, 
                                              freqmin=1.5, freqmax=10, 
                                              deconvolve=True, 
                                              time_before=60, time_after=300)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# CHECK IF SEISMIC EVENT IS CONFIRMED OR REJECTED

sacfiles = glob.glob( outdir + "/*HZ*.sac" )  
st = Stream()
for sacfile in sacfiles:
  st1 = read(sacfile)
  for tr in st1:
    ind = np.where( tr.stats.station == stations )[0][0]
    tr.stats.distance = calc_vincenty_inverse( stlats[ind], stlons[ind], this_event[2], this_event[1] )[0]
    st += tr

fig = plt.figure()
st.plot(type='section', method='full', orientation='vertical', time_down=True, linewidth=.25, grid_linewidth=.25, show=False, fig=fig)
ax = fig.axes[0]
transform = blended_transform_factory(ax.transData, ax.transAxes)
for tr in st:
  ax.text(tr.stats.distance / 1e3, 1.0, tr.stats.station, rotation=270, va="bottom", ha="center", transform=transform, zorder=10)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()


flag = input("\n+ Is this a local earthquake? (yes/no): ")
while flag!="yes" and flag !="no":
  flag = input("+ Is this a local earthquake? (yes/no): ")
  if flag=="yes" or flag =="no":
    break

if flag=="yes":
  update_event_status(origin_time=this_event[0], status="confirmed", table="LOC")
elif flag=="no":
  update_event_status(origin_time=this_event[0], status="rejected", table="LOC")
  subprocess.call( "rm -r %s" % (outdir) , shell=True)
  if os.path.isdir("%s/%s" % (pygema_path, outdir) ):
    subprocess.call( "rm -r %s/%s" % (pygema_path, outdir) , shell=True)

  print("\nEOF ! ! \n\n")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# LOCALIZE SEISMIC EVENT USING NONLINLOC

if flag=="yes":
  flag = input("\n+ Do you want to localize using nonlinloc? (yes/no): ")
  while flag!="yes" and flag !="no":
    flag = input("+ Do you want to localize using nonlinloc? (yes/no): ")
    if flag=="yes" or flag =="no":
      break

  if flag=="yes":
    if len(glob.glob(outdir+'/*') ) > 0:
      cmd1 = "python %s/pygema/signal_processing/execSG2KMOD2.py" % (pygema_path)
      subprocess.call(cmd1, shell=True)
    else:
      print("+ No sacfiles found in event directory ... ")


  locfile = "%s/%s.loc" % (outdir, outdir)
  if os.path.isfile(locfile):
    try:
      flag = input("\n+ Do you want to update status to manual and recalculate event magnitude? (yes/no): ")
      while flag!="yes" and flag !="no":
        flag = input("+ Do you want to update status to manual and recalculate event magnitude? (yes/no): ")
        if flag=="yes" or flag =="no":
          break

      if flag=="yes":
        evtime_old = this_event[0]
        evtime, evlon, evlat, evdepth, evrms, everrx, everry, everrz = read_nonlinloc_output(locfile)
        evmag = get_local_magnitude(evtime, evlon, evlat, evdepth, networks, stations, stlons, stlats, stalts, freqmin=3., freqmax=10.)
        print("new evmag = %.1f" % (evmag) )

        # update database
        update_event_localization(evtime_old, evtime, evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, table="LOC")
        update_event_status(origin_time=evtime, status="manual", table="LOC")

        # send confirmation email
        figsdir = "%s/GEMA/PyGEMA/pygema/web/PyGema_Web/PyGema_Web/static/reports/manual/%s" % (os.getenv("HOME"), evtime.strftime("%Y-%m-%dT%H:%M:%S+00:00") )
        if not os.path.isdir(figsdir):
          os.mkdir(figsdir)

        plot_map(evlon, evlat, evdep, dlon=0.3, dlat=0.3, add_topo=True, show_plot=False, save_plot=True, savedir=figsdir, format='jpg', dpi=200)
        plot_waveforms(evtime, evlon, evlat, freqmin=2, freqmax=10, show_plot=False, save_plot=True, savedir=figsdir, format='jpg', dpi=200)
        message = "[ localized event ] \n Origin Time: %s   mag = %.1f \n evlon = %.4f deg;  evlat = %.4f deg;  evdep = %.1f km \n errX = %.1f km; errY = %.1f km; errZ = %.1f km" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evmag, evlon, evlat, evdepth, everrx, everry, everrz)
        send_email_with_attached_files(message, figsdir=figsdir)


        outdir_new = "%s" % (evtime.strftime("%Y-%m-%d-%H%M%S") )
        os.rename(outdir+"/"+outdir+".pick", outdir+"/"+outdir_new+".pick")
        os.rename(outdir+"/"+outdir+".loc", outdir+"/"+outdir_new+".loc")
        os.rename(outdir, outdir_new)
        subprocess.call( "rm -r %s/*sac" % (outdir_new) , shell=True)
        subprocess.call( "cp -r %s %s/waveforms/" % (outdir_new, pygema_path) , shell=True)
        subprocess.call( "rm -r %s" % (outdir_new) , shell=True)
        if os.path.isdir(pygema_path+"/waveforms/"+outdir):
          subprocess.call( "rm -r %s/waveforms/%s" % (pygema_path, outdir) , shell=True)

      elif flag=="no":
        subprocess.call( "rm -r %s/*sac" % (outdir) , shell=True)
        subprocess.call( "cp -r %s %s/waveforms/" % (outdir, pygema_path) , shell=True)
        subprocess.call( "rm -r %s" % (outdir) , shell=True)

        flag = input("\n+ Do you want to reject event? (yes/no): ")
        while flag!="yes" and flag !="no":
          flag = input("+ Do you want to reject event? (yes/no): ")
          if flag=="yes" or flag =="no":
            break

        if flag=="yes":
          update_event_status(origin_time=this_event[0], status="rejected", table="LOC")


    except:
      print("+ It seems that the event file is empty (%s) ..." % (locfile) )

  else:
    if os.path.isdir( outdir ):
      subprocess.call( "rm -r %s/*sac" % (outdir) , shell=True)
      subprocess.call( "cp -r %s %s/waveforms/" % (outdir, pygema_path) , shell=True)
      subprocess.call( "rm -r %s" % (outdir) , shell=True)

  print("\nEOF ! ! \n\n")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




