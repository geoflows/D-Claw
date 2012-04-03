#!/usr/bin/env python
#
# Update netcdf Makefile and setplot.py to implement geoclaw netcdf routines
#
# Works in 2d geoclaw applications.
#
# Writes output to Makefile, and renames original Makefile.original.

import sys,string,os

os.rename('Makefile', 'Makefile.original')
ifile = open('Makefile.original','r')
ofile = open('Makefile','w')

nclean = 0
for line in ifile:
    if string.find(line,'$(GEOLIB)/valout_geo.f') > -1:
        print 'replacing geoclaw/2d/lib/valout_geo.f'
        ofile.write('  $(GEOLIB)/valout_nc_geo.f \\\n')
    elif string.find(line,'FFLAGS ?=') > -1:
        print 'replacing compiler flag "FFLAGS ?="'
        ofile.write('FFLAGS ?=-I/usr/include  -lnetcdf -lnetcdff\n')
    else:
        ofile.write('%s' % line)

print 'new version is now in Makefile (old version in Makefile.original)'

ifile.close()
ofile.close()


os.rename('setplot.py', 'setplot.py.original')
ifile = open('setplot.py.original','r')
ofile = open('setplot.py','w')

for line in ifile:
    if string.find(line,'return plotdata') > -1:
        print 'replacing setplot.py'
        ofile.write('    plotdata.format = "netcdf"\t\t\t# use netcdf output\n    return plotdata\n')
    else:
        ofile.write('%s' % line)

print 'new version is now in setplot.py (old version in setplot.py.original)'

ifile.close()
ofile.close()
