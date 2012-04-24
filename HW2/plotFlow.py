#This function reads in the output of projection_method.f90 and plots a contour map of the different variables

#Jonathan Varkovitzky
#April 23, 2012

from numpy import *
from matplotlib import *
from matplotlib import rc
rc("text",usetex=True)
from pylab import *
from scipy.special import *
import time

################################
## Read in simulation results ##
################################

def readIn(i,j,fname):
    uvp = loadtxt(fname)
    return(uvp)

##########################
## Make Frames and Plot ##
##########################

def makeFrames(uvp,i,j):
    framenos = size(uvp)/(i*j)
    print "The number of frames to be plotted is %s..." %framenos

    for frameno in range(0,2):
        plotFrame(uvp,frameno,framenos)

################
## Plot Frame ##
################

def plotFrame(uvp,k,framenos):
    print "Extracting data for frame number %s" %k
    nPts = i*j
    xArr = uvp[k*nPts:(k+1)*nPts,0]
    yArr = uvp[k*nPts:(k+1)*nPts,1]
    uArr = uvp[k*nPts:(k+1)*nPts,2]
    vArr = uvp[k*nPts:(k+1)*nPts,3]
    pArr = uvp[k*nPts:(k+1)*nPts,4]
    
    #Reshape vector into properly sized array
    x = xArr.reshape(j,i)
    y = yArr.reshape(j,i)
    u = uArr.reshape(j,i)
    v = vArr.reshape(j,i)
    p = pArr.reshape(j,i)
   

    print v
    if k == 1:
        figure(1)
        clf()
        pcolormesh(x,y,v)
        colorbar()
    return(x)
##################
## Main Program ##
##################
close('all')

fname = '_output/UVP.dat'
i = 200
j = 80
uvp = readIn(200,80,fname)
x = makeFrames(uvp,i,j)
show()
