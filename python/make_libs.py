"""
Performs 'make .objs' in each library.
This should be run before make_all.py so module files *.mod will exist.

Sends output and errors to separate files to simplify looking for errors.
"""

import os,sys,glob
try:
    import subprocess
except:
    print '*** Error: require subprocess module from Python 2.4 or greater'
    raise ImportError()


def make_libs(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$CLAW')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)

    print "Will 'make .objs' in library subdirectories of"
    print "    ", rootdir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()
    
    fname_output = 'make_libs_output.txt'
    fout = open(fname_output, 'w')
    fout.write("ALL OUTPUT FROM MAKING LIBRARIES\n\n")

    fname_errors = 'make_libs_errors.txt'
    ferr = open(fname_errors, 'w')
    ferr.write("ALL ERRORS FROM MAKING LIBRARIES\n\n")


    #os.chdir(rootdir)
    goodlist = []
    badlist = []

    liblist = """clawpack/1d/lib clawpack/2d/lib amrclaw/2d/lib
                 geoclaw/2d/lib""".split()

    for lib in liblist:
        libdir = os.path.join(rootdir,lib)
        if os.path.isdir(libdir):
            os.chdir(libdir)
            fout.write("\n=============================================\n")
            fout.write(libdir)
            fout.write("\n=============================================\n")
            ferr.write("\n=============================================\n")
            ferr.write(libdir)
            ferr.write("\n=============================================\n")

            # flush I/O buffers:
            fout.flush()
            ferr.flush()

            print "Running 'make .objs' in ",libdir
            job = subprocess.Popen(['make','.objs'], \
                             stdout=fout, stderr=ferr)
            return_code = job.wait()
            if return_code == 0:
                print "   Successful completion"
                goodlist.append(lib)
            else:
                print "   *** Errors encountered: see ", fname_errors
                badlist.append(lib)
        else:
            print "*** Library not found:",libdir

    
    print ' '
    print 'Libraries created:'
    for d in goodlist:
        print '   ',d
    print ' '
    
    if len(badlist) > 0:
        print 'Errors encountered in the following libraries:'
        for d in badlist:
            print '   ',d
        print ' '
    
    fout.close()
    ferr.close()
    print 'For all output see ', fname_output
    print 'For all errors see ', fname_errors

if __name__=='__main__':
    make_libs(sys.argv[1:])
