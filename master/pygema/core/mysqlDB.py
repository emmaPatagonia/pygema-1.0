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

import sys, os, glob, datetime, MySQLdb, imp, time, socket
from obspy.core import UTCDateTime, read, Stream, Trace
from matplotlib.dates import date2num, num2date, DateFormatter
import numpy as np

from pygema.core.check_myip_connection import retrieve_mysqlDB_credentials
from pygema.core.email import send_email

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def select_event_list(starttime, endtime, status='automatic', table="LOC", max_gap=270, max_depth=250):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()

  if table=="AUTOLOC":
    #flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude FROM AUTOLOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f"), status )
    flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude FROM AUTOLOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S"), status )
  elif table=="LOC":
    #flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude FROM LOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f"), status )
    flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude FROM LOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S"), status )

  cursor.execute(flag)
  res = cursor.fetchall()
  db.close()

  events_list = []
  for row in res:
    time, lon, lat, depth, nstats, gap, rms, status, magnitude = row
    if gap<=max_gap and depth<=max_depth:
      events_list.append( [UTCDateTime(time), lon, lat, depth, int(nstats), gap, rms, status, magnitude] )


  return events_list


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def select_events_manual_loc(starttime, endtime, table="LOC", max_gap=270, max_depth=250):
  status = "manual"

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()

  if table=="AUTOLOC":
    #flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude, dx, dy, dz FROM AUTOLOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f"), status )
    flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude, dx, dy, dz FROM AUTOLOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S"), status )
  elif table=="LOC":
    #flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude, dx, dy, dz FROM LOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f"), status )
    flag = "SELECT time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude, dx, dy, dz FROM LOC where time BETWEEN '%s' AND '%s' && status='%s'" % ( starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S"), status )


  cursor.execute(flag)
  res = cursor.fetchall()
  db.close()

  events_list = []
  for row in res:
    time, lon, lat, depth, nstats, gap, rms, status, magnitude, dx, dy, dz = row
    if gap<=max_gap and depth<=max_depth:
      events_list.append( [UTCDateTime(time), lon, lat, depth, int(nstats), gap, rms, status, magnitude, dx, dy, dz] )



  return events_list



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def insert_rsam_warnings(station, rsam_values, rsam_times, rsam_orange_warning, rsam_yellow_warning):

  credentials_file = retrieve_mysqlDB_credentials()


  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  inds_yellow = np.where( (np.array(rsam_values)>=rsam_yellow_warning) & (np.array(rsam_values)<rsam_orange_warning)  )[0]
  inds_orange = np.where(  (np.array(rsam_values)>=rsam_orange_warning)  )[0]
  flags = []

  if len(inds_orange)>0:
    for ind in inds_orange:
      color = "ORANGE"
      comment = "RSAM value of %.1f (RSAM > %.1f)" % ( rsam_values[ind], rsam_orange_warning )
      flag = "INSERT INTO  WARNING (time, color, commentary, station) VALUES ('%s', '%s', '%s', '%s')" % ( (rsam_times[ind]).strftime("%Y-%m-%d %H:%M:%S UTC"), color, comment, station )
      flags.append(flag)
      #message = "[%s warning] %s" % (color, comment)
      #send_email(message)

  elif len(inds_yellow)>0:
    for ind in inds_yellow:
      color = "YELLOW"
      comment = "RSAM value of %.1f (%.1f > RSAM > %.1f)" % ( rsam_values[ind], rsam_orange_warning, rsam_yellow_warning )
      flag = "INSERT INTO  WARNING (time, color, commentary, station) VALUES ('%s', '%s', '%s', '%s')" % ( (rsam_times[ind]).strftime("%Y-%m-%d %H:%M:%S UTC"), color, comment, station )
      flags.append(flag)
      #message = "[%s warning] %s" % (color, comment)
      #send_email(message)

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  for flag in flags:
    try:
      cursor.execute(flag)
      db.commit()
    except:
      continue

  db.close()




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def select_rsam_warnings(station, starttime, endtime):

  credentials_file = retrieve_mysqlDB_credentials()


  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  #flag = "SELECT * FROM WARNING WHERE station= '%s' && time BETWEEN '%s' and '%s'" % (station, starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f") )
  flag = "SELECT * FROM WARNING WHERE station= '%s' && time BETWEEN '%s' and '%s'" % (station, starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S") )
  cursor.execute(flag)
  res = cursor.fetchall()
  db.close()

  colors = []; times = []; times_matplotlib = []
  for row in res:
    time, color, comment, station = row
    colors.append( color )
    times.append( UTCDateTime(time) )
    times_matplotlib.append( date2num(time) )

  return colors, times, times_matplotlib




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def insert_triggers_stalta(stations_list, triggers_list):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  for station, triggers in zip(stations_list, triggers_list):
    for trigger in triggers:
      try:
        #flag = "INSERT INTO AUTOTRIGGERS (station, trig_on, trig_off) VALUES ('%s', '%s', '%s')" % ( station, trigger[0].strftime("%Y-%m-%d %H:%M:%S.%f"), trigger[1].strftime("%Y-%m-%d %H:%M:%S.%f") )
        flag = "INSERT INTO AUTOTRIGGERS (station, trig_on, trig_off) VALUES ('%s', '%s', '%s')" % ( station, trigger[0].strftime("%Y-%m-%d %H:%M:%S"), trigger[1].strftime("%Y-%m-%d %H:%M:%S") )
        cursor.execute(flag)
        db.commit()
      except:
        continue

  db.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def select_triggers_stalta(station, starttime, endtime):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  #flag = "SELECT station, trig_on, trig_on FROM AUTOTRIGGERS  WHERE station='%s' && trig_on BETWEEN '%s' AND '%s'" % ( station, starttime.strftime("%Y-%m-%d %H:%M:%S.%f"), endtime.strftime("%Y-%m-%d %H:%M:%S.%f") )
  flag = "SELECT station, trig_on, trig_on FROM AUTOTRIGGERS  WHERE station='%s' && trig_on BETWEEN '%s' AND '%s'" % ( station, starttime.strftime("%Y-%m-%d %H:%M:%S"), endtime.strftime("%Y-%m-%d %H:%M:%S") )
  cursor.execute(flag)
  res = cursor.fetchall()
  db.close()

  triggers_list = []
  for row in res:
    station, on, off = row
    triggers_list.append( [UTCDateTime(on), UTCDateTime(off)] )

  return triggers_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def ask_last_triggers_stalta_inserted(N = 10):
  credentials_file = retrieve_mysqlDB_credentials()
  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T
  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  flag = "SELECT * FROM AUTOTRIGGERS ORDER BY trig_on DESC LIMIT %i" % ( N )
  cursor.execute(flag)
  res = cursor.fetchall()
  for line in res:
    print( "[%s]  %s    %s " % (line[0],line[1].strftime("%Y-%m-%d %H:%M:%S"),line[2].strftime("%Y-%m-%d %H:%M:%S")) )

  db.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def insert_event(evtime, evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status="automatic", table="LOC"):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  if table=="AUTOLOC":
    #flag = "INSERT INTO AUTOLOC (time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude) VALUES ('%s', %f, %f, %f, %i, %f, %f,'%s',%f)" % (evtime.strftime("%Y-%m-%d %H:%M:%S.%f"), evlon, evlat, evdep, evnstats, evgap, evrms, status, evmag)
    flag = "INSERT INTO AUTOLOC (time, longitude, latitude, depth, number_of_stations, gap, rms, status, magnitude) VALUES ('%s', %f, %f, %f, %i, %f, %f,'%s',%f)" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evlon, evlat, evdep, evnstats, evgap, evrms, status, evmag)
  elif table=="LOC":
    #flag = "INSERT INTO LOC (time, longitude, latitude, depth, number_of_stations, gap, rms, magnitude, status) VALUES ('%s', %f, %f, %f, %i, %f, %f, %.1f,'%s')" % (evtime.strftime("%Y-%m-%d %H:%M:%S.%f"), evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status)
    flag = "INSERT INTO LOC (time, longitude, latitude, depth, number_of_stations, gap, rms, magnitude, status) VALUES ('%s', %f, %f, %f, %i, %f, %f, %.1f,'%s')" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evlon, evlat, evdep, evnstats, evgap, evrms, evmag, status)

  cursor.execute(flag)
  db.commit()

  #message = "[ %s event ] \n Origin Time: %s   mag = %.1f \n evlon = %.4f deg;  evlat = %.4f deg;  evdep = %.1f km\n recorded by %i stations\n gap = %.1f deg; rms = %.1f" % (status, evtime.strftime("%Y-%m-%d %H:%M:%S"), evmag, evlon, evlat, evdep, evnstats, evgap, evrms)
  #send_email(message)

  db.close()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def update_event_status(origin_time, status, table="LOC"):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  if table=="AUTOLOC":
    #flag = "UPDATE AUTOLOC SET status='%s' WHERE  time='%s'" % (status, origin_time.strftime("%Y-%m-%d %H:%M:%S.%f") )
    flag = "UPDATE AUTOLOC SET status='%s' WHERE  time='%s'" % (status, origin_time.strftime("%Y-%m-%d %H:%M:%S") )
  elif table=="LOC":
    flag = "UPDATE LOC SET status='%s' WHERE  time='%s'" % (status, origin_time.strftime("%Y-%m-%d %H:%M:%S") )
    #flag = "UPDATE LOC SET status='%s' WHERE  time='%s'" % (status, origin_time.strftime("%Y-%m-%d %H:%M:%S") )

  cursor.execute(flag)
  db.commit()
  db.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def update_event_localization(evtime_old, evtime, evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, table="LOC"):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  if table=="AUTOLOC":
    #flag = "UPDATE  AUTOLOC  SET  time  = '%s',  longitude  = %f,  latitude  = %f,  depth  = %f,  rms  = %.4f,  dx  = %f,  dy  = %f,  dz  = %f,  magnitude  = %.1f WHERE  LOC . time  = '%s'" % (evtime.strftime("%Y-%m-%d %H:%M:%S.%f"), evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, evtime_old.strftime("%Y-%m-%d %H:%M:%S.%f") )
    flag = "UPDATE  AUTOLOC  SET  time  = '%s',  longitude  = %f,  latitude  = %f,  depth  = %f,  rms  = %.4f,  dx  = %f,  dy  = %f,  dz  = %f,  magnitude  = %.1f WHERE  LOC . time  = '%s'" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, evtime_old.strftime("%Y-%m-%d %H:%M:%S") )
  elif table=="LOC":
    #flag = "UPDATE  LOC  SET  time  = '%s',  longitude  = %f,  latitude  = %f,  depth  = %f,  rms  = %.4f,  dx  = %f,  dy  = %f,  dz  = %f,  magnitude  = %.1f WHERE  LOC . time  = '%s'" % (evtime.strftime("%Y-%m-%d %H:%M:%S.%f"), evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, evtime_old.strftime("%Y-%m-%d %H:%M:%S.%f") )
    flag = "UPDATE  LOC  SET  time  = '%s',  longitude  = %f,  latitude  = %f,  depth  = %f,  rms  = %.4f,  dx  = %f,  dy  = %f,  dz  = %f,  magnitude  = %.1f WHERE  LOC . time  = '%s'" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evlon, evlat, evdepth, evrms, everrx, everry, everrz, evmag, evtime_old.strftime("%Y-%m-%d %H:%M:%S") )

  cursor.execute(flag)
  db.commit()
  db.close()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def update_origin_time(evtime_old, evtime_new, table="LOC"):

  credentials_file = retrieve_mysqlDB_credentials()

  host, port, user, password, database = np.loadtxt(credentials_file, dtype='str').T

  db = MySQLdb.connect(host=host, port=int(port), user=user, passwd=password, db=database)
  cursor = db.cursor()
  if table=="AUTOLOC":
    flag = "UPDATE  AUTOLOC  SET  time  = '%s' WHERE  LOC . time  = '%s'" % (evtime_new.strftime("%Y-%m-%d %H:%M:%S"), evtime_old.strftime("%Y-%m-%d %H:%M:%S") )
  elif table=="LOC":
    flag = "UPDATE  LOC  SET  time  = '%s'  WHERE  LOC . time  = '%s'" % (evtime_new.strftime("%Y-%m-%d %H:%M:%S"), evtime_old.strftime("%Y-%m-%d %H:%M:%S.%f") )

  cursor.execute(flag)
  db.commit()
  db.close()

