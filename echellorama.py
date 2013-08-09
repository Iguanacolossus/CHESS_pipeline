
# coding: utf-8
#!/usr/bin/env python

## example:
##           $ python squashmatrix2.py 'filename.fits'

import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.text as text
import sys
import astropy.io.fits as fits
import numpy as np
from scipy import signal
from gi.repository import Gtk, GObject
import pygtk
from matplotlib.figure import Figure
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo as FigureCanvas
#from matplotlib.backends.backend_GTKAgg import FigureCanvasGTK as FigureCanvas
#import pylab
#import aplpy


class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self,title="Echelle Reduction GUI")
  
   ## setting up GUI canvase ###     
        self.set_default_size(700,600)
    	self.f = Figure(figsize=(5,7), dpi=100)
        self.a = self.f.add_subplot(221)
        self.b = self.f.add_subplot(212)
        self.c = self.f.add_subplot(222)
        self.a.grid(False)
        self.a.set_xlabel("x pixels", fontsize=10)
        self.a.tick_params(axis='both', labelsize=7)
        self.a.set_ylabel("y pixels")
        self.a.set_title("2D raw data")
        self.b.set_title("1D extracted data")
        self.b.set_xlabel('pixels')
        self.b.set_ylabel('intensity')
        self.b.tick_params(axis='both', labelsize=7)
        self.c.set_title('PHD')
        self.c.tick_params(axis='both', labelsize=7)
        
        
        self.canvas = FigureCanvas(self.f)
        self.add(self.canvas)
    # passing in filename or sends error
    	if (len(sys.argv) < 2):
        	raise Exception("Need to pass file name")
    
    # reads in the file name
    	filename = sys.argv[1]

    #this opens up the fits file
    	hdulist = fits.open(filename)
    
    	targname = hdulist[0].header['targname']    
    
    #this picks out the actual data from fits file, and turns it into numpy array
    	scidata = hdulist['sci',2].data 
    	hdulist.close()
    
   # old plot stuff  
    	#plt.subplot(211)
   	#plt.xlabel('pixels(x)')
    	#plt.ylabel('pixels(y)')
    	#plt.title(targname)
    	#img = plt.imshow(scidata, vmin = 0, vmax = 255)
    
    # sends 2d data to my gui plotting funtion	
    	self.update_plot(scidata)
    
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
	minv = chunk[0]
	maxv = chunk[1]
 
  # plotting for my own sake to check myself while making code
	#plt.subplot(212)
        #plt.plot(y,xrev)
	#plt.ylabel('pixels(y)')
	#plt.xlabel('Intensity')
    #drawing box around my chunk to check makesure I have a good portion
	#plt.hlines(chunk,[-1000],[5000],color='r')
        #plt.vlines([-1000,5000],minv,maxv,color='r')
	#plt.show()
    
    #cutting out the chunk of data that i selected
	xchunk = x[max(x)-chunk[1]:max(x)]
	index1 = int(max(x)-chunk[1])
	index2 = int(max(x))
	ychunk = y[index1:index2]
    #reversing x for the sake of the plot
    	xrevchunk = xchunk[::-1]
    
####  this is my manual code for making a zero backgrounds then finding the value at each peak############################
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
    
    #plotting chunk of data with lines through the centers of the orders to double check how the peak finder did
    	#plt.subplot(211)   
    	#plt.plot(xchunk,ychunk)
    	#plt.vlines(xchunk[peakind],0,2000,color='purple',label='centers')
    	#plt.title('chunk of data with centers found by find_peaks_cwt()')
    
    #find w, the widths of the orders
    # first i make an array of the difference between peaks
    	w=[]
    	for i in range(1,len(peakind)):
    		t=peakind[i]-peakind[i-1]
    		w.append(t)
    
    
    	avew=sum(w)/len(w)
    # i have to add an extra w at the end of the array to make it the right size i (hopefully this is kosher)
    	maxw=max(w)-4
    	w.append(maxw)
    
    # placeing verticle lines in the visualization where the widths would stand
    	#wvlines=[]
    	#b=xchunk[peakind[0]]-avew/2
    	#wvlines.append(b)
    	#for i in np.arange(len(peakind)):
    		#f=xchunk[peakind[i]]+w[i]/2
    		#wvlines.append(f)
    	#plt.vlines(wvlines,0,2500,color='r',linestyle='--',label='w boundries')
    	#plt.legend(loc=1)
    
    ####extraction of orders
    
    	plt.figure()
    
    ### making arrays of 1s ans 0s and extracting the 1d orders into a dictionary called  oneDorders
   	zeros=np.zeros((len(x),1))
    	index = range(0,len(w))
    	reindex = index[::-1]
    
    	oneDorders = {}
    	for i in reindex:
    		zeros1=np.copy(zeros)
    		zeros1[ 1024 - (np.sum(w[(i):18])) : 1024 - np.sum(w[(i+1):18]) ] = 1
    		twoD = scidata*zeros1
    		if i == 16:
    			plt.subplot(211)
    			plt.imshow(twoD)
    			plt.title('sample of one of the orders in 1d and 2d') 	
    		Y=[]
     # making 2d orders in to 1d orders
    		for j in range(0,len(scidata[0])):
        		t = np.sum(twoD[:,j])
        		Y.append(t)
        	oneDorders['order'+str(i)]=Y
     # sending plotting info to update_1dplot for gui
    	x = np.linspace(0,len(scidata[0]), num = len(scidata[0]))
    	#plt.subplot(212)
    	#plt.plot(x,oneDorders['order16'])
    	odo=oneDorders['order16']
    	self.update_1dplot(odo,x)

         
    ### turning 2D orders into 1D orders
    #o1=[]
    #1Dorders={}
   # for i in range(0,len(scidata[0])):
    #    t = np.sum([:,i])
     #   o1.append(t)
#    x1 = np.linspace(0,len(order1[0]), num = len(order1[0]))
#    plt.subplot(212)
#    plt.title('one order')
#    plt.plot(x1,o1)
    	
###this next part was me fitting with air glow#######
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
    #plt.figure()
    #for i in range(1,len(w)):
    #	w[i]=w[0]+w[i]
    #w[0]=w[0]*2
    #plt.plot(xchunk[peakind],w,marker='o',linestyle='none',color='purple',label="space between orders")
    #plt.ylabel('w')
    #plt.xlabel('pixels along yaxis of detector')
#    plt.ylabel('order spaceing by pixel')
#    plt.title('center possitions with repect to verticle postion on detector')
    
    ##fitting orderspaceing to line
    #fit=np.polyfit(xchunk[peakind],w,2)
    #linefit=fit[2] + fit[1]*xchunk + fit[0]*(xchunk**2)
    #deg0 = str.format('{0:.2f}', fit[0])
    #deg1 = str.format('{0:.2f}', fit[1])
    #deg2 = str.format('{0:.2f}', fit[2])
    #plt.plot(xchunk, linefit, color='r', linestyle='--', label='$x^2$'+deg0+'$+x$'+deg1+'$+$'+deg2+'')
#    plt.legend(loc=4)
   

### fake PDH stuff####\
	x=np.arange(10)
	y=(.2*x)**3
	self.update_PHDplot(x,y)
    
    
    	#plt.show()
    	
  ## plotting in gui
    def update_plot(self, scidata):
        self.plt= self.a.imshow(scidata, vmin = 0, vmax = 255)
    	
        self.canvas.draw()
    def update_1dplot(self,odo,x):
    	
    	self.plt=self.b.plot(x,odo)
        
        self.canvas.draw()
    def update_PHDplot(self,x,y):
        self.plt=self.c.plot(x,y)
        self.canvas.draw
    



if __name__ == '__main__':
    win = MyWindow()
    win.connect("delete-event", Gtk.main_quit)
    win.show_all()
    Gtk.main()

