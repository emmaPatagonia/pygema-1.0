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
import sys, os

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def find_pygema_parent_directory(version="PyGEMA"):
  parent_directory = []
  for path in sorted(sys.path):
    if (version in path):
      parent_directory.append(path)

  return parent_directory[0]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def load_station_metadata(stations_file=None):
  if stations_file is None:
    pygema_path = find_pygema_parent_directory()
    stations_file = pygema_path + "/src/stationALL.lst"

  networks, stations = np.loadtxt(stations_file, usecols=(0,1), dtype='str').T
  stlons, stlats, stalts = np.loadtxt(stations_file, usecols=(3,2,4)).T

  return networks, stations, stlons, stlats, stalts



 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def load_rsam_parameters(rsam_params_file=None):
  if rsam_params_file is None:
    pygema_path = find_pygema_parent_directory()
    rsam_params_file = pygema_path + "/src/params_rsam.txt"

  twin, freqmin, freqmax = np.loadtxt(rsam_params_file).T

  return twin, freqmin, freqmax




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def load_ssam_parameters(ssam_params_file=None):
  if ssam_params_file is None:
    pygema_path = find_pygema_parent_directory()
    ssam_params_file = pygema_path + "/src/params_ssam.txt"

  twin, freqmin, freqmax, nfreqs = np.loadtxt(ssam_params_file).T
  nfreqs = int(nfreqs)

  return twin, freqmin, freqmax, nfreqs



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def load_stalta_parameters(stalta_file=None):
  if stalta_file is None:
    pygema_path = find_pygema_parent_directory()
    stalta_file = pygema_path + "/src/params_stalta.txt"

  len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, stream_length = np.loadtxt(stalta_file).T

  return len_sta, len_lta, trig_on, trig_off, freqmin, freqmax, stream_length




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def load_warnings_rsam(station, warning_file_rsam=None):
  if warning_file_rsam is None:
    pygema_path = find_pygema_parent_directory()
    warning_file_rsam = pygema_path + "/src/warning_file_rsam.txt"

  stations = np.loadtxt(warning_file_rsam, usecols=(0,), dtype='str' ).T
  lower_corners, upper_corners, yellow_warnings, orange_warnings = np.loadtxt(warning_file_rsam, usecols=(1,2,3,4)).T
  ind = np.where(stations==station)[0][0]
  lower_corner, upper_corner, yellow_warning, orange_warning = lower_corners[ind], upper_corners[ind], yellow_warnings[ind], orange_warnings[ind]

  return lower_corner, upper_corner, yellow_warning, orange_warning





