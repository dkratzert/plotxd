#!/usr/bin/python
#-*- coding: utf-8 -*-

import sys, os
import re
import optparse
import csv
from optparse import OptionParser
progress = True
try:
   from progressbar import ProgressBar
except:
   progress = False

headers = [
   [('Pval',1), ('Kappa',2), ("Kappa'",4), ('Net charge',5)],
   [('D11+',1), ('D11-',2), ('D10',3)],
   [('Q20',1), ('Q21+',2), ('Q21-',3), ('Q22+',4), ('Q22-',5)],
   [('O30',1), ('O31+',2), ('O31-',3), ('O32+',4), ('O32-',5), ('O33+',6), ('O33-',7)],
   [('H40',1),('H41+',2),('H41-',3),('H42+',4),('H42-',5),('H43+',6),('H43-',7),('H44+',8),('H44-',9)],
   [('X', 1), ('Y', 2), ('Z', 3), ('OZ', 4), ('ISO', 5)],
   [('U11',1),('U22',2),('U33',3),('U12',4),('U13',5),('U23',6)],
]
  
def output_csv(output, ofile):
   csv_writer = csv.writer(ofile, dialect='excel',delimiter=';',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
   head = []
   for column in output[0]:
      head.append(column.keys()[0])
   csv_writer.writerow(head)   

   for cycle in output:
      line = []
      for column in cycle:
         line.append(column.values()[0])
      csv_writer.writerow(line)
   

def output_standard(output, ofile):
      for column in output[0]:  # Column names in the first line
         print >> ofile, "%s\t"%column.keys()[0],
      print >> ofile

      for cycle in output:
         for column in cycle:
            print >> ofile, "%s\t"%column.values()[0], 
         print >> ofile

def find_line(lines, token):
   for i, line in enumerate(lines):
      if line.find(token) >= 0:
         return i


def find_table(lines, num):
   j = 0
   tofind = 'Table %d.'%num
   for i in range(len(lines)):
      if lines[i].find(tofind) > 0:  # find the Table head first
         j=i
         break
   while lines[j].find('---') != 0: # find the next line with ---
      j+=1
   return j+1

# Returns a list of all atoms in the geofile
def find_atoms(geofile):
   gfile = None
   try:
      gfile = file(geofile, 'r')
   except IOError, e:
      print "Error: No such file or directory: '%s'"%geofile
      sys.exit(1)

   start = 0
   lines = gfile.readlines()
   start = find_table(lines, 1)

   i = 0
   atoms = []
   regexp = re.compile('\s+([^\s]+)\s')
   while 1:
      line = lines[start+i]
      if line[:3] == '---':  #atoms start at the third line after ---
         break
      atoms.append( re.match(regexp, line).group(1) )
      i += 1
   
   gfile.close()
   return atoms

# Searches for the given atom in a table and
# returns the line where the atom was found or "None" otherwise.
def find_atomline(lines, atom, start):
   i = start
   while 1:
      line = lines[i]
      pos = line.find(atom)
      if pos >= 0 and pos <= 2:
         return line
      i += 1
      if line.find('---') == 0 or line.strip() == '': #Hier ist die Tabelle zu Ende
         return None

def split_atomline(line):
   # Sometimes columns with negative value can be too 
   # long so they lack a white space in between.
   # When this happens, we split them manually
   columns = line.split()
   ret = columns
   if line.find(')-') > 0:
      for i,col in enumerate(columns):
         pos = col.find(')-')
         if pos > 0:
            a,c = col.split(')-')
            ret.insert(i, '-'+c)
            ret.insert(i, a+')')
            ret.pop(i+2)
   return ret
      

def extract_features(lines, atom, cyclecount):
      cycle = []
      #Loops through all tables containing a header
      #In each table the features of the current atom are extracted
      for tnum in range(len(headers)-2): #loop over table number
         start = find_table(lines, tnum+1)
         atom_line = find_atomline(lines, atom, start)
         if atom_line == None: continue #Atom not in this table

         cols = atom_line.split()
         for c in headers[tnum]: #column
            value = cols[c[1]]
            if not options.deviation:
               pos = value.find('(')
               if pos >=0: 
                  value = value[:pos]
            cycle.append( {c[0]: value } )

      #These tables need extra handling
      for tnum, token in enumerate(('Coordinates', 'Uij values')):
         start = find_line(lines, token)
         atom_line = find_atomline(lines, atom, start)
         if atom_line == None: continue

         cols = split_atomline( atom_line )
         for c in headers[tnum+5]:
            tmp = c[1]
            value = cols[tmp]
            if not options.deviation:
               pos = value.find('(')
               if pos >=0: 
                  value = value[:pos].rstrip()
            cycle.append( {c[0]: value } )
      
      cycle.insert(0, {'Nummer': "%03d"%cyclecount})
      return cycle

# Assembles the features of an atom
# and returns a list of column dicts, 
# where the key of each dict is the name of the feature.
def handle_geofiles(atom, geofile_pattern):
   cyclecount = 1 
   output = []
   cycle = []

   while 1: #Loops as long as geo.out files are found
      fname = geofile_pattern%cyclecount
      if not os.access(fname, os.R_OK): # check if geo file exists
         break

      geofile = file(fname, 'r')
      lines = tuple(geofile.readlines())

      cycle = extract_features(lines, atom, cyclecount)
      output.append( cycle )
      cyclecount += 1
   return output
   

def handle_atom(atom, geofile_pattern):
   output = handle_geofiles(atom, geofile_pattern)

   ofile=None
   if options.csvout:
      ofile = file(options.outdir+os.sep+'%s.csv'%atom, 'w')
      output_csv(output, ofile)
   else:
      ofile = file(options.outdir+os.sep+'%s.txt'%atom, 'w')
      output_standard(output, ofile)

   ofile.close()

def extract_geofile_pattern(geofile):
   match = re.match("(.*?)([0]+)(\d)(\w*.\w+)", geofile)
   out = ""
   if match:
      groups = match.groups()
      start = ""
      if groups[0]:
         start = groups[0]
      zeros = groups[1]

   out = start + "%0" + str(len(zeros)+1) + "d" + groups[3]
   return out

if __name__ == '__main__':
   parser = OptionParser(version="%prog-03")
   parser.add_option("-s", "--std-dev",  dest="deviation", action="store_true",
      help="Include standard deviation in output. Default: False", default=False)
   parser.add_option("-c", "--csv",  dest="csvout", action="store_true",
      help="Enables output format CSV. Default: False", default=False)
   parser.add_option("-i", "--input", default="xd01_geo.out", dest="geofile", 
      help="Name of the first input file. Default: %s"%"xd01_geo.out")
   parser.add_option("-o", "--outdir", default="pyout/", dest="outdir",
      help="Output directory. Default: pyout/")

   (options, args) = parser.parse_args()

   try:
      os.mkdir(options.outdir)
   except:
      pass

   atoms = find_atoms(options.geofile)

   if progress:
      bar = ProgressBar(width=70)
   num_atom = len(atoms)
   i = 1.0
   for atom in atoms: #Each atom will be written to an extra file
      if progress:
         bar( (i/num_atom) * 100 )
      else:
         print atom,
         sys.stdout.flush()
      pattern = extract_geofile_pattern( options.geofile )
      handle_atom(atom, pattern)
      i += 1
   else:
      print





# vim:sts=3:ts=3:sw=3
