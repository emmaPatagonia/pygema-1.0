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
from pygema.read.parameters import find_pygema_parent_directory

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def retrieve_seiscomp3_datadirs(datadir=None):

  s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  s.connect(('8.8.8.8', 1))
  local_ip = s.getsockname()[0]

  if local_ip=="192.168.1.201" and socket.gethostname()=='maniedba':
    buffer_dir = "/home/gema/seiscomp3/var/lib/seedlink/buffer"
    archive_dir = "/home/gema/seiscomp3/var/lib/archive"

  elif socket.gethostname()=='sirius' or socket.gethostname()=='tremor':
    if not datadir is None:
      buffer_dir = "%s/seiscomp_data_buffer" % (datadir)
      archive_dir = "%s/seiscomp_data_archive" % (datadir)
    else:
      buffer_dir = "/home/gema/SHARED/seiscomp_data_buffer"
      archive_dir = "/home/gema/SHARED/seiscomp_data_archive"

  elif socket.gethostbyname(socket.gethostname())=="152.74.5.244" and socket.gethostname()=='patagonia':
    if not datadir is None:
      buffer_dir = "%s/seiscomp_data_buffer" % (datadir)
      archive_dir = "%s/seiscomp_data_archive" % (datadir)
    else:
      buffer_dir = "/store2/SHARED/seiscomp_data_buffer"
      archive_dir = "/store2/SHARED/seiscomp_data_archive"

  else:
    if not datadir is None:
      buffer_dir = "%s/seiscomp_data_buffer" % (datadir)
      archive_dir = "%s/seiscomp_data_archive" % (datadir)
    else:
      buffer_dir = "%s/seiscomp3_msdata/seiscomp_data_buffer" % (os.getenv("HOME"))
      archive_dir = "%s/seiscomp3_msdata/seiscomp_data_archive" % (os.getenv("HOME"))
      #raise Exception('Check your internet connection or mount seiscomp3 database !!')

  return buffer_dir, archive_dir


def retrieve_mysqlDB_credentials():

  s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
  s.connect(('8.8.8.8', 1))
  local_ip = s.getsockname()[0]

  if local_ip=="192.168.1.201" and socket.gethostname()=='maniedba':
    pygema_path = find_pygema_parent_directory()
    credentials_file = pygema_path + "/src/credentials_db_mysql_maniedba"

  elif socket.gethostname()=='sirius':
    pygema_path = find_pygema_parent_directory()
    credentials_file = pygema_path + "/src/credentials_db_mysql_sirius"

  elif socket.gethostname()=='tremor':
    pygema_path = find_pygema_parent_directory()
    credentials_file = pygema_path + "/src/credentials_db_mysql_tremor"

  elif socket.gethostbyname(socket.gethostname())=="152.74.5.244" and socket.gethostname()=='patagonia':
    pygema_path = find_pygema_parent_directory()
    credentials_file = pygema_path + "/src/credentials_db_mysql_patagonia"

  else:
    pygema_path = find_pygema_parent_directory()
    credentials_file = pygema_path + "/src/credentials_db_mysql"

  return credentials_file

