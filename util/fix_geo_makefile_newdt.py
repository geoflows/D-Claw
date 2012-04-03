#!/usr/bin/env python

"""
Fix geoclaw Makefile for changes introduced for automatic choice of
refinement ratios on each level. 

Usage:
   $ python $CLAW/util/fix_geo_makefile_newdt.py Makefile_name
  
"""
import sys,string,os

def fixmake(file_list):

    for file in file_list:
        oldfile = open(file).read()
        if oldfile.find('tick_geo.f') > -1:
            print "*** %s already in proper form?  Found tick_geo.f" % file
            print "*** Not changing ",file
        else:
            os.rename(file, file+'.original')
            lines = open(file+'.original','r').readlines()
            ofile = open(file,'w')
    
            #for i,line in enumerate(lines):
            for line in lines:
    
                if string.find(line,'upbnd_geo.f') > -1:
                    print 'Adding new _geo routines'
                    ofile.write('%s' % line)
                    ofile.write(' $(GEOLIB)/tick_geo.f \\\n')
                    ofile.write(' $(GEOLIB)/setgrd_geo.f \\\n')
                    ofile.write(' $(GEOLIB)/gfixup_geo.f \\\n')
                    ofile.write(' $(GEOLIB)/ginit_geo.f \\\n')
                    ofile.write(' $(GEOLIB)/getmaxspeed_geo.f \\\n')
                elif string.find(line,'tick.f') > -1:
                    pass
                elif string.find(line,'setgrd.f') > -1:
                    pass
                elif string.find(line,'gfixup.f') > -1:
                    pass
                elif string.find(line,'ginit.f') > -1:
                    pass
                else:
                    ofile.write('%s' % line)
    
            print 'new version is now in %s ' % file
            print '  (old version in %s.original)' % file
    
            ofile.close()

if __name__=='__main__':
    import sys
    print "Calling fixmake on ",sys.argv[1:]
    fixmake(sys.argv[1:])

