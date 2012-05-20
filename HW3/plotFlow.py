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
    framenos = shape(uvp)[0]/(i*j)
    print "The number of frames to be plotted is %s..." %framenos

    for frameno in range(0,nframes):
        if frameno == nframes-1:
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
   
    (xb,ub,vb) = blasiusVals()
    nu = 1e-3
    eta = y[:,0]/(nu*x[:,0])**(1/2)
    print shape(ub)
    print "Josh said to look at  ",shape(eta)


    if k == nframes-1:
                
        figure(1)
        clf()
        title('X Velocity Component')
        xlabel('x')
        ylabel('y')
        contourf(x,y,u)
        colorbar()
        axis([4,10,0,0.4])
        savefig('u_full.png')
        axis([4.8,8,0,0.1])
        savefig('u_zoom.png')
        """
        figure(2)
        clf()
        title('Y Velocity Component')
        xlabel('x')
        ylabel('y')
        contourf(x,y,v)
        plot(xb,vb)
        colorbar()
        axis([4,10,0,0.4])
        savefig('v_full.png')
        axis([4.8,5.4,0.0,0.1])
        savefig('v_zoom')

        figure(3)
        clf()
        title('Pressure')
        xlabel('x')
        ylabel('y')
        contourf(x,y,p)
        colorbar()
        axis([4,10,0,0.4])
        savefig('p_full.png')
        axis([4.8,5.4,0.0,0.1])
        savefig('p_zoom')


        figure(4)
        clf()
        title('Blasius Comparison, u Velocity')
        xlabel('eta')
        ylabel('u')
        plot(eta*10,u[:,-1],'b')
        plot(xb,ub,'r')
        axis([0,8,0,1.1])
        legend(("Computational Solution","Blasius Solution"),loc='lower right')
        savefig('blasius_u.png')

        figure(5)
        clf()
        title('Blasius Comparison, v Velocity')
        xlabel('eta')
        ylabel('v')
        plot(eta*10,v[:,-1],'b')
        plot(xb,vb/4,'r')
        axis([0,8,-0.05,0.1])
        legend(("Computational Solution","Blasius Solution"),loc='lower right')
        savefig('blasius_v.png')
        """
        figure(10)
        clf()
        omega = 10
        U = 1
        nu = 1#10**(-3)
        t = 5
        yE = arange(0,81)
        uExact = U*exp(-1*yE*sqrt(omega/(2*nu)))#*cos(omega*t-yE*sqrt(omega/(2*nu)))
        title('u(y)')
        xlabel('u(y)')
        ylabel('y')
        print "the shapes of y and UExact are: %r %s"%(shape(yE),shape(uExact))
        print uExact
        plot(u[:,10],'b')
        plot(zeros(80),'k--')
        plot(uExact,yE,'g')
        axis([0,80,-0.001,0.001])
    return(x)


##################
## Get  Blasius ##
##################
def blasiusVals():
    blasius = loadtxt('blasius2.txt')
#    blasius = blasius/10.
    U = 1.
    nu = 1e-3
    x = blasius[:,0]
    u = U*blasius[:,1]
    v = 1./sqrt(2)*(U*nu/x)**(1/2)*blasius[:,2]/10


    return(x,u,v)
######################
## readIn and Parse ##
######################
def readInUVP(i,j,fname):
    uvp = loadtxt(fname)

    print "Reading in frame %s" %fname
    k = shape(uvp)[0]/(i*j)-1
    print k
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

    return(x,y,u,v,p)
################
## Skin Coeff ##
################
def skinFric():
    
    ic = 200
    jc = 80
    
    im = 300
    jm = 60
    
    ifi = 400
    jfi = 80

    (xc,yc,uc,vc,pc) = readInUVP(ic,jc,'_output/UVP_coarse.dat')
    (xm,ym,um,vm,pm) = readInUVP(im,jm,'_output/UVP_med.dat')
    (xf,yf,uf,vf,pf) = readInUVP(ifi,jfi,'_output/UVP_fine.dat')

    cf = 2*nu/1.*(uf[0,:]-uf[-1,:])/(yf[0,:]-yf[1,:])
    cm = 2*nu/1.*(um[0,:]-um[-1,:])/(ym[0,:]-ym[1,:])
    cc = 2*nu/1.*(uc[0,:]-uc[-1,:])/(yc[0,:]-yc[1,:])
    
    (xb,ub,vb) = blasiusVals()


    figure(6)
    title('Skin-Friction Coefficient')
    plot(xc[0,:],cc,'b')
    plot(xm[0,:],cm,'r')
    plot(xf[0,:],cf,'g')
    plot(xf[0,:],0.664*(xf[0,:]/nu)**(-1./2),'k')
    legend(['Coarse','Medium','Fine','Exact'],loc = 'upper left')
#    savefig('skin_fric.png')

##################
## Main Program ##
##################
close('all')

fname = '_output/UVP.dat'
i = 200
j = 80
nu = 1e-3
print "Reading in %s" %fname
uvp = readIn(i,j,fname)
nframes = shape(uvp)[0]/(i*j)
print "Making Images"
x = makeFrames(uvp,i,j)
#skinFric()
#show()


