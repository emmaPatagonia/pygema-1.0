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
from obspy.core import UTCDateTime, read, Stream, Trace
import sys, os, glob, datetime, MySQLdb, imp, time, socket, io

from pygema.read.parameters import load_station_metadata
from pygema.core.mysqlDB import select_event_list, select_events_manual_loc
from pygema.read.parameters import find_pygema_parent_directory

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def kmlexport_stations(kmlfilename='stations.kml', outdir=None):
  networks, stations, stlons, stlats, stalts = load_station_metadata()
  pygema_path = find_pygema_parent_directory()
  if outdir is None:
    outdir = pygema_path + "/src/kml"

  if not os.path.isdir(outdir):
    os.makedirs(outdir)

  f = open(outdir+'/'+kmlfilename, 'w')
  f.write("<?xml version='1.0' encoding='UTF-8'?>\n")
  f.write("<kml xmlns='http://earth.google.com/kml/2.1'>\n")
  f.write("<Document>\n")
  f.write("   <name>" + kmlfilename +"</name>\n")

  f.write("   <Style id='normalPlacemark'>\n")
  f.write("     <IconStyle>\n")
  f.write("       <scale>1.2</scale>\n")
  f.write("       <Icon>\n")
  f.write("         <href>http://maps.google.com/mapfiles/kml/pal4/icon56.png</href>\n")
  f.write("       </Icon>\n")
  f.write("     </IconStyle>\n")
  f.write("   </Style>\n")

  f.write("   <StyleMap id='exampleStyleMap'>\n")
  f.write("     <Pair>\n")
  f.write("       <key>normal</key>")
  f.write("       <styleUrl>#normalPlacemark</styleUrl>\n")
  f.write("     </Pair>\n")
  f.write("   </StyleMap>\n")

  for network, station, stlon, stlat, stalt in zip(networks, stations, stlons, stlats, stalts):
    f.write("   <Placemark>\n")

    f.write("       <name>" + network + '.' + station + "</name>\n")
    f.write("       <styleUrl>#exampleStyleMap</styleUrl>\n")
    #f.write("       <description> ... write any description here ... </description>\n")

    f.write("       <ExtendedData>\n")
    f.write("           <Data name='Network'>\n")
    f.write("               <value>\n")
    f.write("                   " + network + "\n")
    f.write("               </value>\n")
    f.write("           </Data>\n")
    f.write("           <Data name='Station'>\n")
    f.write("               <value>\n")
    f.write("                   " + station + "\n")
    f.write("               </value>\n")
    f.write("           </Data>\n")
    f.write("           <Data name='Longitude'>\n")
    f.write("               <value>\n")
    f.write("                   " + str(stlon) + "\n")
    f.write("               </value>\n")
    f.write("           </Data>\n")
    f.write("           <Data name='Latitude'>\n")
    f.write("               <value>\n")
    f.write("                   " +  str(stlat) + "\n")
    f.write("               </value>\n")
    f.write("           </Data>\n")
    f.write("           <Data name='Altitude'>\n")
    f.write("               <value>\n")
    f.write("                   " + str(stalt) + "\n")
    f.write("               </value>\n")
    f.write("           </Data>\n")
    f.write("       </ExtendedData>\n")

    f.write("       <Point>\n")
    f.write("           <coordinates>" + str(stlon) + "," + str(stlat) + "," + str(stalt) + "</coordinates>\n")
    f.write("       </Point>\n")

    f.write("   </Placemark>\n")

  f.write("</Document>\n")
  f.write("</kml>\n")
  f.close()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def kmlexport_seismic_events(starttime, endtime, color, status, kmlfilename='seismic_events_automatic.kml', max_gap=270, max_depth=50, outdir=None):
  pygema_path = find_pygema_parent_directory()
  if outdir is None:
    outdir = pygema_path + "/src/kml"

  events = select_event_list(starttime, endtime, status=status, max_gap=max_gap, max_depth=max_depth)
  if status=='manual':
    events = select_events_manual_loc(starttime, endtime, table="LOC", max_gap=max_gap, max_depth=max_depth)
  elif status=='automatic':
    events = select_event_list( starttime, endtime, status='automatic', table="LOC", max_gap=max_gap, max_depth=max_depth)

  if len(events)>0:
    events = np.array(events)[np.argsort(np.array(events)[:, 0])]

    if not os.path.isdir(outdir):
      os.makedirs(outdir)

    f = open(outdir+'/'+kmlfilename, 'w')
    f.write("<?xml version='1.0' encoding='UTF-8'?>\n")
    f.write("<kml xmlns='http://earth.google.com/kml/2.1'>\n")
    f.write("<Document>\n")
    f.write("   <name>" + kmlfilename +"</name>\n")

    f.write("   <Style id='normalPlacemark'>\n")
    f.write("     <IconStyle>\n")
    f.write("       <scale>1.4</scale>\n")
    f.write("       <Icon>\n")
    if color=='red' and status=='manual':
      f.write("         <href>http://maps.google.com/mapfiles/kml/paddle/red-stars.png</href>\n")
    elif color=='yellow' and status=='manual':
      f.write("         <href>http://maps.google.com/mapfiles/kml/paddle/ylw-stars.png</href>\n")
    elif color=='green' and status=='manual':
      f.write("         <href>http://maps.google.com/mapfiles/kml/paddle/grn-stars.png</href>\n")
    elif status=='automatic':
      f.write("         <href>http://maps.google.com/mapfiles/kml/paddle/red-stars-lv.png</href>\n")

    f.write("       </Icon>\n")
    f.write("     </IconStyle>\n")
    f.write("   </Style>\n")

    f.write("   <StyleMap id='exampleStyleMap'>\n")
    f.write("     <Pair>\n")
    f.write("       <key>normal</key>")
    f.write("       <styleUrl>#normalPlacemark</styleUrl>\n")
    f.write("     </Pair>\n")
    f.write("   </StyleMap>\n")

    for event in events:
      evtime = event[0]
      evlon = event[1]
      evlat = event[2]
      evdep = event[3]
      evmag = event[8]
      status = event[7]

      f.write("   <Placemark>\n")

      #f.write("       <name>" + evtime.strftime("Origin Time: %Y-%m-%d %H:%M:%S.%f") + "</name>\n")
      f.write("       <styleUrl>#exampleStyleMap</styleUrl>\n")
      #f.write("       <description> ... write any description here ... </description>\n")

      f.write("       <ExtendedData>\n")
      f.write("           <Data name='Origin time'>\n")
      f.write("               <value>\n")
      f.write("                   " + evtime.strftime("%Y-%m-%d %H:%M:%S.%f") + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      f.write("           <Data name='Magnitude'>\n")
      f.write("               <value>\n")
      f.write("                   " + str(evmag) + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      f.write("           <Data name='Longitude'>\n")
      f.write("               <value>\n")
      f.write("                   " + str(evlon) + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      f.write("           <Data name='Latitude'>\n")
      f.write("               <value>\n")
      f.write("                   " +  str(evlat) + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      f.write("           <Data name='Depth (km)'>\n")
      f.write("               <value>\n")
      f.write("                   " + str(evdep) + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      f.write("           <Data name='Status'>\n")
      f.write("               <value>\n")
      f.write("                   " + status + "\n")
      f.write("               </value>\n")
      f.write("           </Data>\n")
      if status=='manual':
        errx = event[9]
        erry = event[10]
        errz = event[11]
        f.write("           <Data name='Error X (km)'>\n")
        f.write("               <value>\n")
        f.write("                   " + str(errx) + "\n")
        f.write("               </value>\n")
        f.write("           </Data>\n")
        f.write("           <Data name='Error Y (km)'>\n")
        f.write("               <value>\n")
        f.write("                   " + str(erry) + "\n")
        f.write("               </value>\n")
        f.write("           </Data>\n")
        f.write("           <Data name='Error Z (km)'>\n")
        f.write("               <value>\n")
        f.write("                   " + str(errz) + "\n")
        f.write("               </value>\n")
        f.write("           </Data>\n")

      f.write("       </ExtendedData>\n")

      f.write("       <Point>\n")
      f.write("           <coordinates>" + str(evlon) + "," + str(evlat) + "," + str(evdep) + "</coordinates>\n")
      f.write("       </Point>\n")

      f.write("   </Placemark>\n")

    f.write("</Document>\n")
    f.write("</kml>\n")
    f.close()

  else:
    f = open(outdir+'/'+kmlfilename, 'w')
    f.write("<?xml version='1.0' encoding='UTF-8'?>\n")
    f.write("<kml xmlns='http://earth.google.com/kml/2.1'>\n")
    f.write("<Document>\n")
    f.write("   <name>" + kmlfilename +"</name>\n")

    f.write("   <Style id='normalPlacemark'>\n")
    f.write("     <IconStyle>\n")
    f.write("       <scale>1.0</scale>\n")
    if color=='red':
      f.write("       <color>501400FF</color>\n")
    elif color=='yellow':
      f.write("       <color>5014F0FF</color>\n")
    elif color=='green':
      f.write("       <color>5014B400</color>\n")

    f.write("       <Icon>\n")
    f.write("         <href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n")
    f.write("       </Icon>\n")
    f.write("     </IconStyle>\n")
    f.write("   </Style>\n")

    f.write("   <StyleMap id='exampleStyleMap'>\n")
    f.write("     <Pair>\n")
    f.write("       <key>normal</key>")
    f.write("       <styleUrl>#normalPlacemark</styleUrl>\n")
    f.write("     </Pair>\n")
    f.write("   </StyleMap>\n")

    f.write("</Document>\n")
    f.write("</kml>\n")
    f.close()
