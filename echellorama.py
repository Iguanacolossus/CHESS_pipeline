

#!/usr/bin/env python

import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
import matplotlib.figure as fig

# syntax: $ python squashmatrix2.py 'filename.fits'
import sys

import astropy.io.fits as fits

#import pylab
#import aplpy

import numpy as np
from scipy import signal


def main():

    if (len(sys.argv) < 2):
        raise Exception("Need to pass file name")
    
    # reads in the file name
    filename = sys.argv[1]

    #this opens up the fits file
    hdulist = fits.open(filename)
    
    targname = hdulist[0].header['targname']
    print('Target Name', targname)
    
    
    scidata = hdulist['sci',2].data #this picks out the actual data from fits file, and turns it into numpy array
    
    hdulist.close()
    
    
    plt.subplot(211)
    plt.xlabel('pixels(x)')
    plt.ylabel('pixels(y)')
    img = plt.imshow(scidata, vmin = 0, vmax = 255)
    
    # making an empty array to put the new values in
    y=[]
    
    # filling in y with the sums of the rows of scidata
    for i in range(0,len(scidata[0])):
        t = np.sum(scidata[i,:])
        y.append(t)
        

    # making an x axis with same dementions as y
    x = np.linspace(0,len(scidata[0]), num = len(scidata[0]))

    #reversing x for the sake of the visualization
    xrev = x[::-1]

    # the orders seem to blend around lyman alpha so i have to select
    # the biggest chunk I could use to get a delta y between the
    # orders. i wil have to make this selectable on the GUI some how
    chunk = [0,500]
    
    # cuting off back ground here to get zeros and peaks. i will ahve
    # to make this selectable some how on the GUI
    bgcutoff = 475
    minv = chunk[0]
    maxv = chunk[1]
    plt.subplot(212)
    plt.plot(y,xrev)
    plt.ylabel('pixels(y)')
    plt.xlabel('Intensity')
    #plt.vlines(bgcutoff,[0],1000,label='background cutoff')
    plt.hlines(chunk,[-1000],[5000],color='r')
    plt.vlines([-1000,5000],minv,maxv,color='r')
    plt.show()
    
    #cutting out chunk of data
    xchunk = x[max(x)-chunk[1]:max(x)]
    index1 = int(max(x)-chunk[1])
    index2 = int(max(x))
    ychunk = y[index1:index2]
    #reversing x for the sake of the plot
    xrevchunk = xchunk[::-1]
    #plt.subplot(123)
####  this is my manual code for making a zero backgrouns then finding the value at each peak############################
#    #creating a new baseline
#   for i in range(0,xchunk.size):
#        if ychunk[i] < bgcutoff:
#            ychunk[i] = 0
#
#
#    #plt.plot(ychunk,xrevchunk)
#    #plt.show()
#
#    bounds = []
#    a = ychunk
#    # find the bounds of a order by indice in the array
#    for i in range(1, len(a)):
#        if a[i-1] == 0 and a[i] > 0:
#           # print(i , "begin hump index")
#            bounds.append(i)
#        if a[i-1] > 0 and a[i] == 0:
#           # print(i-1, "end hump index")
#            bounds.append(i - 1)
#
#    # now 'bounds' is a list of the indices of the array were a #order begins and ends
#  
#    # 'bounds' but now grouped
#    grpbounds = list(group(bounds, 2))
#
#    centers = []
#    # finding max with in the bouds of the orders and reading them #in
#    # to 'centers' list
#    for i in range(len(grpbounds)):
#        centers.append(max(a[grpbounds[i][0]:grpbounds[i][1]+1]))
#   
# 
#    #print a.index(centers)
    
# defining a function to group the indecies of 'bounds' to # sets of
# (begin order, end order)
#def group(lst,n):
#    for i in range(0, len(lst), n):
#       val=lst[i:i+n]
#       if len(val) == n:
#          yield tuple(val)
############################################################################################################################3
    plt.figure(figsize=(7.5,8.4))
    #using scipy.signal.find_peaks_cwt() to find centers of orders.  this required scipy version .11.0 or greater
    peakind=signal.find_peaks_cwt(ychunk,np.arange(3,15))
    
    plt.subplot(211)   
    plt.plot(xchunk,ychunk)
    plt.vlines(xchunk[peakind],0,3000,color='purple')
    plt.title('chunk of data with centers found by find_peaks_cwt()')
    
    #find w
    w=[]
    for i in range(1,len(peakind)):
    	t=peakind[i]-peakind[i-1]
    	w.append(t)
    
    avew=sum(w)/len(w)
    #w.append(avew)
    maxw=max(w)
    w.append(maxw)
    
    wvlines=[]
    b=xchunk[peakind[0]]-avew/2
    wvlines.append(b)
    for i in np.arange(len(peakind)):
    	f=xchunk[peakind[i]]+w[i]/2
    	wvlines.append(f)
    plt.vlines(wvlines,0,2800,color='r',linestyle='--')
    
    ####extraction
    plt.figure()
    zeros=np.zeros((len(x),1))
    zeros[1024-max(w):1024]=1
    order1=scidata*zeros
    
    
    o1=[]
    for i in range(0,len(scidata[0])):
        t = np.sum(scidata[:,i])
        o1.append(t)
    x1 = np.linspace(0,len(order1[0]), num = len(order1[0]))
    plt.subplot(212)
    plt.title('one order')
    plt.plot(x1,o1)
    	
    
####### fitting centers to line#####
    #plt.subplot(324)
    #plt.plot(np.zeros(len(peakind)),xchunk[peakind],marker='o',linestyle='none',color='purple',label="order centers")
#    plt.xlabel('yaxis on detector')
#    plt.ylabel('order center index')
#    #plt.title('center possitions with repect to verticle postion on detector')
    
#    print len(range(int(min(xchunk)),int(max(xchunk)),len(peakind)))
#    print len(xchunk[peakind])
    
    #fit=np.polyfit(np.zeros(len(peakind)),xchunk[peakind],2)
     
    
    #xfit=np.zeros(len(peakind))
    #print xchunk[peakind]
    #linefit=fit[2] + fit[1]*xfit + fit[0]*(xfit**2)
#    deg0 = str.format('{0:.2f}', fit[0])
#    deg1 = str.format('{0:.2f}', fit[1])
#    deg2 = str.format('{0:.2f}', fit[2])
    #plt.plot(np.zeros(len(peakind)),linefit,color='r',linestyle='--')#,label='$x^2$'+deg0+'$+x$'+deg1+'$+$'+deg2+'')
#    #plt.legend(loc=2)
    #
#   #arange= int(min(xchunk[peakind]))
#   # brange= int(max(xchunk[peakind]))
#   # crange= int(len(xchunk))+1
#   # print len(range(arange,brange,crange))
    
#    #print len(xchunk[peakind])
#    #plt.plot(range(arange,brange,crange),xchunk[peakind],marker='o',linestyle='none',color='purple',label="centers of orders")
    
#    # filling in for the missing centers with the new centers funtion(linefit)
#    plt.subplot(323)
#    allcen=fit[2] + fit[1]*x + fit[0]*(x**2)
    
#    allcen=[i for i in allcen if i < max(x)]
#    plt.vlines(x[allcen],0,3000,color='purple')
#    plt.plot(x,y)
    
#   # print allcen
#  #  print peakind
    
#    plt.subplot(322)

###### finding and fitting spacing between orders####    
#    orderspace=[]
#    for i in range(1,len(peakind)):
#    	t=peakind[i]-peakind[i-1]
#    	orderspace.append(t)
    	
#    peakind.remove(peakind[0])
#    #print len(peakind)
    plt.figure()
    #for i in range(1,len(w)):
    #	w[i]=w[0]+w[i]
    #w[0]=w[0]*2
    plt.plot(xchunk[peakind],w,marker='o',linestyle='none',color='purple',label="space between orders")
    plt.ylabel('w')
    plt.xlabel('pixels along yaxis of detector')
#    plt.ylabel('order spaceing by pixel')
#    plt.title('center possitions with repect to verticle postion on detector')
    
    ##fitting orderspaceing to line
    fit=np.polyfit(xchunk[peakind],w,2)
    linefit=fit[2] + fit[1]*xchunk + fit[0]*(xchunk**2)
    deg0 = str.format('{0:.2f}', fit[0])
    deg1 = str.format('{0:.2f}', fit[1])
    deg2 = str.format('{0:.2f}', fit[2])
    plt.plot(xchunk, linefit, color='r', linestyle='--', label='$x^2$'+deg0+'$+x$'+deg1+'$+$'+deg2+'')
#    plt.legend(loc=4)
   

    
    
    
    plt.show()
    


if __name__ == '__main__':
    main()
