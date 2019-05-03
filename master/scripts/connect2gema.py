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
sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# montar discos con datos
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.connect(('8.8.8.8', 1))
local_ip = s.getsockname()[0]

if socket.gethostbyname(socket.gethostname())=="152.74.5.244" and socket.gethostname()=='patagonia':
  cmd1 = "sshfs gema@gemaudec.dyndns.org:/home/gema/seiscomp3/var/lib/seedlink/buffer /store2/SHARED/seiscomp_data_buffer -p 2221"
  cmd2 = "sshfs gema@gemaudec.dyndns.org:/home/gema/seiscomp3/var/lib/archive /store2/SHARED/seiscomp_data_archive -p 2221"
elif socket.gethostname()=='sirius' or socket.gethostname()=='tremor':
  cmd1 = "sshfs gema@192.168.1.201:/home/gema/seiscomp3/var/lib/seedlink/buffer /home/gema/SHARED/seiscomp_data_buffer"# -p 2221"
  cmd2 = "sshfs gema@192.168.1.201:/home/gema/seiscomp3/var/lib/archive /home/gema/SHARED/seiscomp_data_archive"# -p 2221"
elif local_ip=="192.168.1.201" and socket.gethostname()=='maniedba':
  raise Exception("You are at MANIEDBA SERVER ... not necesary to mount a virtual disk ... ")
else:
  outdir_buffer = "%s/seiscomp3_msdata/seiscomp_data_buffer" % (os.getenv("HOME"))
  outdir_archive = "%s/seiscomp3_msdata/seiscomp_data_archive" % (os.getenv("HOME"))
  cmd1 = "sshfs gema@gemaudec.dyndns.org:/home/gema/seiscomp3/var/lib/seedlink/buffer %s -p 2221" % (outdir_buffer)
  cmd2 = "sshfs gema@gemaudec.dyndns.org:/home/gema/seiscomp3/var/lib/archive %s -p 2221" % (outdir_archive)

subprocess.call(cmd1, shell=True)
subprocess.call(cmd2, shell=True)
