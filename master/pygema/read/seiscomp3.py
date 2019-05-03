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

import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
from obspy.core import UTCDateTime, read, Stream

from pygema.core.check_myip_connection import retrieve_seiscomp3_datadirs

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_streams_buffer(station, channels, starttime, endtime, num_msfiles=2, datadir=None):
  """ Script para leer la base de datos de Seiscomp3 (buffer)
  + station (str):        nombre de la estacion
  + channels (str):       canal(es) de la estacion (e.g. HHZ, BH*, *H[ZNE])
  + starttime (datetime): tiempo inicial
  + endtime (datetime):   tiempo final
  + num_msfiles (int):    numero de miniseed files para leer datos (notar que cada archivo tiene una duracion aleatoria). Se recomienda 1 (default) o 2.
  """

  buffer_dir, archive_dir = retrieve_seiscomp3_datadirs(datadir=datadir)

  path = '%s/%s/segments' % (buffer_dir, station)
  name_list = os.listdir(path)
  full_list = [os.path.join(path,i) for i in name_list]
  time_sorted_list = sorted(full_list, key=os.path.getmtime)
  if num_msfiles==0:
    msfiles = time_sorted_list
  else:
    msfiles = time_sorted_list[-num_msfiles::]

  st = Stream()
  for msfile in msfiles:
    st += read(msfile)

  st.trim(starttime, endtime)
  st = st.select(channel=channels)

  return st

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_streams_archive(network, station, channels, starttime, endtime, datadir=None):
  """ Script para leer la base de datos de Seiscomp3 (archive)
  + network (str):        nombre de la red
  + station (str):        nombre de la estacion
  + channels (list):      lista de canales (e.g. [HHZ, HHN, HHE]). No se aceptan wildcards
  + starttime (datetime): tiempo inicial
  + endtime (datetime):   tiempo final
  """

  buffer_dir, archive_dir = retrieve_seiscomp3_datadirs(datadir=datadir)

  msfiles = []
  this_day = starttime
  while this_day <= endtime:
    for channel in channels:
      path = '%s/%s/%s/%s/%s.D/*%s' % (archive_dir, this_day.strftime("%Y"), network, station, channel, this_day.strftime("%Y.%03j") )
      msfile = glob.glob(path)
      if len(msfile)>0:
        msfiles.append(msfile[0])

    this_day += 86400

  st = Stream()
  for msfile in msfiles:
    st += read(msfile)

  st.trim(starttime, endtime)

  return st


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_streams_seiscomp3(networks, stations, starttime, endtime, only_vertical_channel=True, merge_method=None, remove_traces_with_gaps=False, datadir=None):
  """ Script para leer la base de datos de Seiscomp3.
  + network (list):        lista de la redes
  + station (list):        lista de la estaciones
  + starttime (datetime): tiempo inicial
  + endtime (datetime):   tiempo final
  + datadir (str):        directorio donde se encuentra la carpeta archive y buffer (check mount_maniedba_seiscomp3_database.py)
  """

  st = Stream()
  gaps = []
  utc_now = UTCDateTime().now()

  for network, station in zip(networks, stations):

    # read archive data
    if utc_now - starttime > 86400:
      if only_vertical_channel:
        st1 = get_streams_archive(network, station, ["BHZ"], starttime, endtime, datadir=datadir)
      else:
        st1 = get_streams_archive(network, station, ["BHZ", "BHN", "BHE"], starttime, endtime, datadir=datadir)

      gaps1 = st1.get_gaps()
      if len(gaps1) > 0 and remove_traces_with_gaps:
        gaps.append(gaps1)
      else:
        st += st1

    else:
      if only_vertical_channel:
        st1 = get_streams_archive(network, station, ["BHZ"], starttime, endtime)
      else:
        st1 = get_streams_archive(network, station, ["BHZ", "BHN", "BHE"], starttime, endtime)

      gaps1 = st1.get_gaps()
      if len(gaps1) > 0 and remove_traces_with_gaps:
        gaps.append(gaps1)
      else:
        st += st1

   # read data from buffer
    if only_vertical_channel:
      st2 = get_streams_buffer(station, "BHZ", starttime, endtime, num_msfiles=3)
    else:
      st2 = get_streams_buffer(station, "BH[ZNE]", starttime, endtime, num_msfiles=3)

      gaps2 = st2.get_gaps()
      if len(gaps2) > 0 and remove_traces_with_gaps:
        gaps.append(gaps2)
      else:
        st += st2

  # combine all data readed
  if len(st)>0:
    st.sort()

    if merge_method==0:
      st.merge(method=0)
    elif merge_method==1:
      st.merge(method=1, interpolation_samples=-1, fill_value='interpolate')
    
  return st, gaps



