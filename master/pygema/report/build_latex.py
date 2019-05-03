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

from obspy.core import UTCDateTime
import sys, os, glob, datetime, MySQLdb, imp, time, socket, subprocess, logging
from pylatex import Document, Section, Subsection, Command, Figure, Package
from pylatex.utils import italic, bold, NoEscape

sys.path.append("%s/GEMA/PyGEMA" % (os.getenv("HOME")) )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def build_report_automatic_event(evtime, evmag, evlon, evlat, evdep, evnstats, evgap, evrms, outpdf):

  doc = Document(outpdf, document_options=['a4paper','14pt'], geometry_options={"margin": "0.7in"}, lmodern=False, textcomp=False, fontenc=None, page_numbers=False)
  doc.packages.add(Package('babel', options='spanish'))
  doc.preamble.append(Command('title', 'REPORTE AUTOMÁTICO DE ACTIVIDAD VOLCÁNICA EN ALTO BIO-BIO'))
  doc.preamble.append(Command('author', 'ENEL - Universidad de Concepción'))
  doc.preamble.append(Command('date', NoEscape('%s' % (UTCDateTime().strftime("%d de %B de %Y (%H:%M:%S UTC)") ) )))
  doc.append(NoEscape(r'\maketitle'))

  doc.append("La CENTRAL DE MONITOREO DE ACTIVIDAD VOLCÁNICA da a conocer la siguiente información obtenida a través de los equipos de monitoreo en tiempo real de ENEL - UNIVERSIDAD DE CONCEPCIÓN:\n\n")
  doc.append("Hoy %s a las %s hora UTC (hora local = UTC - 03:00), las estaciones sismológicas instaladas en Alto Bio-Bio detectaron automáticamente un evento sísmico de magnitud local %.1f, el cual está localizado a %.1f km de profundidad. Los datos de localización son los siguientes: \n\n" % (evtime.strftime("%d de %B de %Y"), evtime.strftime("%H:%M:%S"), evmag, evdep))

  with doc.create(Figure(position = 'ht!')) as fig:
    fig.add_image("figs/map_automatic.jpg", width=NoEscape(r"0.8\textwidth") )

  with doc.create(Figure(position = 'ht!')) as fig:
    fig.add_image("figs/waveforms_automatic.jpg", width=NoEscape(r"0.7\textwidth") )

  doc.append("TIEMPO DE ORIGEN: %s UTC\nMAGNITUD: %.1f\nLONGITUD: %.4f\nLATITUD: %.4f\nPROFUNDIDAD: %.1f km\nESTACIONES: %i\nGAP: %.1f\nRMS: %.1f\n" % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evmag, evlon, evlat, evdep, evnstats, evgap, evrms) )
  
  doc.append("\n\nIMPORTANTE: Los datos mostrados en este reporte son generados automáticamente y deben ser revisados por un analista. Cualquier consulta dirigirla directamente a:\n" )
  doc.append("Email: gemaudec.concepcion@gmail.com")

  doc.generate_pdf(clean_tex = True, compiler='pdflatex')
  doc.generate_tex()










 