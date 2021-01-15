#!/usr/bin/python
"""
gaugedata
=========
   Provides routines for analyzing gauge data
   as output from Geoclaw files: fort.gauge.

   see:
      gaugedata.fortgaugeread
      gaugedata.fortgaugeread1d
      gaugedata.plotgauge
      gaugedata.plotfortgauge
      gaugedata.writegdata
      gaugedata.readgdata
      gaugedata.fortgauge2gdata
      gaugedata.samplesinglegauge
      gaugedata.Lagrangian_Xoft

"""

from pylab import * #I regret this import...need to clean-up at some point
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot
import os
import string

#import pdb


#========================================================================
def fortgaugeread (datafile="fort.gauge",setgaugefile="setgauges.data"):
    """
    Read data from fort.gauge files output by GeoClaw.

    Reads the gauge data and returns a list with mgauge elements.
    Each element of the list is a dictionary containing the data for a particular gauge.

    example:

    for N gauges: allgaugedata=[dictionary_1,...,dictionary_N] where
    dictionary_n.keys() = gauge#,x,y,level,t,q1,...qmq+1.
    where level,t,q1...q mq+1 are numpy arrays. 'q mq+1' is the surface elevation.
    """

    fid=open(setgaugefile)
    inp='#'
    while inp == '#':
        inpl=fid.readline()
        inp=inpl[0]

    inp = fid.readline()
    mgauges=int(inp.split()[0])
    gaugelocs=[]
    linesread=0
    while linesread < mgauges :
        row=string.split(fid.readline())
        if row!=[]:
           gaugelocs.append(row)
           linesread=linesread+1

    fid.close()

    data=loadtxt(datafile)

    allgaugedata=[]

    for n in range(mgauges) :
        onegaugedata=data[mlab.find(data[:,0]==int(gaugelocs[n][0]))]
        dict={}
        dict['gauge']=int(gaugelocs[n][0])
        dict['x']=float(gaugelocs[n][1])
        dict['y']=float(gaugelocs[n][2])
        onegaugedata = data[mlab.find(data[:,0]==dict['gauge'])]
        dict['level']=onegaugedata[:,1]
        dict['t']=onegaugedata[:,2]
        dict['mq'] = len(onegaugedata[0])-3
        for m in range(dict['mq']) :
            dict['q'+str(m+1)]=onegaugedata[:,3 + m]

        allgaugedata.append(dict)

    return allgaugedata
    # end fortgaugeread ======================================================

#=============================================================================
def selectgauge (gaugenumber, allgaugedata=[], datafile="fort.gauge", \
    setgaugefile="setgauges.data"):
    """
    select a single gauge in fort.gauge files output by GeoClaw
    or from an allgaugedata list of dictionaries produced by fortgaugeread
    returns a single  dictionary for a single gauge
    """

    if allgaugedata==[]:
        allgaugedata=fortgaugeread(datafile,setgaugefile)
    for g in range(len(allgaugedata)):
        gnumber = allgaugedata[g]['gauge']
        if gnumber==gaugenumber:
            gg = g
    try:
        gaugedata = allgaugedata[gg]
    except:
        print('Gauge number %i does not exist in %s' % (gaugenumber,datafile))

    return gaugedata



#====================================================================================
def writegdata (gaugedata,fname='',level=1,var='q1'):

    """
    given dictionary gaugedata for data at a single gauge write a file in .gdata format
        ie. with header:
         'gauge#' int
         'level' int
         'x' float
         'y' float
         'rows' int
         'column definitions: 't, variable'
    """

    dict=gaugedata

    if fname=='':
        fname= '%s %03i %s %03i %s' % ('gauge_',dict['gauge'],'_level',level,'.gdata')
        #fname='gauge'+str(dict['gauge'])+'.gdata'

    fid=open(fname,'w')
    keystrlist=['gauge','x','y']
    varstrlist=['t','h','hu','hv','eta','b']

    for stri in keystrlist:
        fid.write("%s %s \n" % (dict[stri],stri))
    sep='   '
    fid.write(sep.join(varstrlist))
    fid.write('\n')

    for tp in range(len(dict['t'])):
        fid.write(" %012.6e %012.6e %012.6e %012.6e %012.6e %012.6e \n" % (dict['t'][tp],dict['h'][tp],dict['hu'][tp],dict['hv'][tp],dict['eta'][tp],dict['b'][tp]))

    fid.close()
    return

##
##
## functions below need to be updated

#=============================================================================
def plotfortgauge (gaugenumber, allgaugedata=[],gaugevar1='t',gaugevar2='q1',\
                   datafile="fort.gauge",setgaugefile="setgauges.data"):
    """
    plot the data output in fort.gauge files output by GeoClaw
    """

    if allgaugedata==[]:
        allgaugedata=fortgaugeread(datafile,setgaugefile)
    for g in range(len(allgaugedata)):
        gnumber == allgaugedata[g]['gauge']
        if gnumber==gaugenumber:
            gg = g
    try:
        plotdata = allgaugedata[gg]
    except:
        print('Gauge number %i does not exist in %s' % (gaugenumber,datafile))

    plotdata1=plotdata[gaugevar1]

    if gaugevar2=='b':
        plotdata2=plotdata['q4']-plotdata['q1']
    elif gaugevar2=='eta':
        plotdata2=plotdata['q4']
    else:
        plotdata2=plotdata[gaugevar2]

    lines=pyplot.plot(plotdata1,plotdata2)
    return lines



#=====================================================================================
def plotgauge (ingauge,var1='t',var2='q1'):

    """
    plot the gauge data in ingauge.
    indata may be a file name for a .gdata file or a gaugedata dictionary
    """

    if type(ingauge)==str:
        gaugedata=readgdata(ingauge)
    else:
        gaugedata=ingauge
    handle=pyplot.plot(gaugedata[var1],gaugedata[var2])
    return handle



#=====================================================================================
def readgdata (fname):

    """
    read the data from a file containing gauge data in the .gdata format
        ie. with header:
         int 'gauge#'
         float 'x'
         float 'y'
         'column definitions: eg., t,h,hu,hv,eta,b'
    """

    gaugedata={}
    fid=open(fname,'r')

    row = string.split(fid.readline())
    gaugedata['gauge'] = int(row[0])
    row = string.split(fid.readline())
    gaugedata['x'] = float(row[0])
    row = string.split(fid.readline())
    gaugedata['y'] = float(row[0])
    row = string.split(fid.readline())

    fid.close()

    data=np.loadtxt(fname,skiprows=4)
    keylist=['t','h','hu','hv','eta','b']
    for i in range(len(keylist)) :
        gaugedata[keylist[i]]=data[:,i]

    fid.close()
    return gaugedata


#================================================================================
def fortgauge2gdata (indatafile="fort.gauge",setgaugefile="setgauges.data"):

    """
    write out data in fort.gauge into separate files in the .gdata format
        ie. with headers:
         int 'gauge#'
         float 'x'
         float 'y'
         't,h,hu,hv,eta,b'
    """

    allgaugedata = fortgaugeread (indatafile,setgaugefile)

    N=len(allgaugedata)

    for ig in range(N):
        writegdata(allgaugedata[ig],fname=allgaugedata[ig]['gauge'])
    return

#================================================================================
def samplesinglegauge (ingauge,ntimes,t0=-inf,tend=+inf,output=True,outname=''):

    """
    subsample the data in ingauge to return a new dictionary
    with output only at ntimes equally spaced points between t0 and tend
    data is linearly interpolated to times.
    ingauge may be a dictionary of gauge data, or an input file name for a .gdata file
    if output is True, a new file is written.
    """

    if type(ingauge)==str:
        gaugedata=readgdata(ingauge)
    else:
        gaugedata=ingauge

    sampledat={}
    sampledat['gauge']=gaugedata['gauge']
    sampledat['x']=gaugedata['x']
    sampledat['y']=gaugedata['y']

    sampledat['t']=[]
    #keylist=['h','hu','hv','eta','b']
    keylist=['q1','q2','q3','q4','q5','q6','q7','q8']
    for key in keylist:
        sampledat[key]=[]

    if t0==-inf:
        t0=gaugedata['t'][0]
    if tend==inf:
        tend=gaugedata['t'][-1]

    t=linspace(t0,tend,ntimes)
    for tn in t:
        sampledat['t'].append(tn)
        if (tn<gaugedata['t'][0])|(tn>gaugedata['t'][-1]):
            for key in keylist :
                sampledat[key].append(nan)
        elif (tn==gaugedata['t'][0]):
            for key in keylist :
                sampledat[key].append(gaugedata[key][0])
        elif (tn==gaugedata['t'][-1]):
            for key in keylist :
                sampledat[key].append(gaugedata[key][-1])
        else:
            itm=mlab.find(gaugedata['t']<tn)[-1]
            itp = itm + 1
            tm=gaugedata['t'][itm]
            tp=gaugedata['t'][itp]
            delt=tp-tm
            for key in keylist:
                dely=gaugedata[key][itp]-gaugedata[key][itm]
                slope=dely/delt
                sampledat[key].append(gaugedata[key][itm] + slope*(tn-tm))

    if output:
        if outname=='':
            outname= '%s%03i_%i_.%s' % ('gauge',sampledat['gauge'],ntimes,'gdata')
        writegdata(sampledat,outname)

    if type(ingauge)==str:
        return
    else:
        return sampledat



def getgaugedata(fortfile='fort.gauges',setgaugefile='setgauges.data'):

    #-----get all gauge data--------------------------
    allgaugedata = fortgaugeread(fortfile,setgaugefile)

    #------find gauge locations and numbers-----------

    fid=open(setgaugefile)
    inp='#'
    while inp == '#':
        inpl=fid.readline()
        inp=inpl[0]

    inp = fid.readline()
    mgauges=int(inp.split()[0])
    gaugelocs=[]
    linesread=0
    while linesread < mgauges :
        row=string.split(fid.readline())
        if row!=[]:
            gaugelocs.append(row)
            linesread=linesread+1
    fid.close()

    xgauges=[]
    ygauges=[]
    gauge_nums=[]
    for gauge in range(mgauges):
        xgauges.append(allgaugedata[gauge]['x'])
        ygauges.append(allgaugedata[gauge]['y'])
        gauge_nums.append(allgaugedata[gauge]['gauge'])

    return allgaugedata,xgauges,ygauges,gauge_nums
#---------------------------------------------------

#===============================================================================

def Lagrangian_Xoft(allgaugedata,xgauges,gauge_nums,x0,t):
    """
    return an array, Xoft that shows position vs time starting at x0
    """
    import pdb

    Xoft=np.array([])
    dt = t[1]-t[0]
    Xoft=np.hstack((Xoft,x0))

    for nt in range(len(t)-1):
        tn = t[nt]
        tnp= t[nt+1]
        x = Xoft[nt]
        #find appropriate gauges
        if xgauges[-1]>x:
            ind = np.where(xgauges>x)[0][0]
        else:
            ind = -1
        indm1 = max(ind-1,0)
        gaugeupper=gauge_nums[ind]
        gaugelower=gauge_nums[indm1]
        gdataupper = selectgauge(gaugeupper,allgaugedata)
        gdatalower = selectgauge(gaugelower,allgaugedata)
        #find position of gauges relative to x(tn)
        xlower=gdatalower['x']
        xupper=gdataupper['x']
        dx = xupper - xlower
        if dx > 0:
            chi = (x-xlower)/dx
        else:
            chi = 1.0
        tlower = gdatalower['t']
        tupper = gdataupper['t']
        #spatially intepolated velocity at tn and xn
        #pdb.set_trace()
        indlower = np.where(tlower>tn)[0][0]
        indupper = np.where(tupper<=tn)[-1][-1]
        uupper = gdataupper['q2'][indupper]/gdataupper['q1'][indupper]
        ulower= gdatalower['q2'][indlower]/gdatalower['q1'][indlower]
        if gdataupper['q1'][indupper]==0.0:
            uupper=0.0
            chi = 0.0
        if gdatalower['q1'][indlower]==0.0:
            ulower=0.0
            chi = 1.0
        u = chi*uupper + (1.0-chi)*ulower
        x = x + u*dt
        #spatially intepolated velocity at tn and xn
        #find appropriate gauges at tn + 0.5dt
        if xgauges[-1]>x:
            ind = np.where(xgauges>x)[0][0]
        else:
            ind = -1
        indm1 = max(ind-1,0)
        gaugeupper=gauge_nums[ind]
        gaugelower=gauge_nums[indm1]
        gdataupper = selectgauge(gaugeupper,allgaugedata)
        gdatalower = selectgauge(gaugelower,allgaugedata)
        #find position of gauges relative to x(tn + 0.5dt)
        xlower=gdatalower['x']
        xupper=gdataupper['x']
        dx = xupper - xlower
        if dx > 0:
            chi = (x-xlower)/dx
        else:
            chi = 1.0
        tlower = gdatalower['t']
        tupper = gdataupper['t']
        #spatially intepolated velocity at tn and xn
        indlower = np.where(tlower>tn+0.5*dt)[0][0]
        indupper = np.where(tupper<=tn+0.5*dt)[-1][-1]
        uupper = gdataupper['q2'][indupper]/gdataupper['q1'][indupper]
        ulower= gdatalower['q2'][indlower]/gdatalower['q1'][indlower]
        if gdataupper['q1'][indupper]==0.0:
            uupper=0.0
            chi = 0.0
        if gdatalower['q1'][indlower]==0.0:
            ulower=0.0
            chi = 1.0
        u2 = chi*uupper + (1.0-chi)*ulower
        x = Xoft[nt] + 0.5*(u + u2)*dt

        #pdb.set_trace()
        Xoft=np.hstack((Xoft,x))

    return Xoft









