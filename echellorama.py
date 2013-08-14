
# coding: utf-8
#!/usr/bin/env python


# syntax:           >>> python squashmatrix2.py 'filename.fits' 'Background_filename.fits'
# 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.text as text
import sys
import astropy.io.fits as fits
import numpy as np
from scipy import signal
from gi.repository import Gtk, GObject

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar





class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self,title="Echelle Reduction GUI")
  
   ## setting up canvase for plotting ###     
        self.set_default_size(1000,700)
    	self.f = Figure(figsize=(5,7), dpi=100)
        self.b = self.f.add_subplot(212) # 1D
        self.c = self.f.add_subplot(231) # PHD
        self.e = self.f.add_subplot(233) # orders 
        self.a = self.f.add_subplot(232) # raw
        self.e.tick_params(axis='both',labelsize=6)
        self.e.set_title("orders")
        self.a.tick_params(axis='both', labelsize=7)
        self.a.set_title("2D raw data")
        self.b.set_title("1D extracted data")
        self.b.set_xlabel('pixels')
        self.b.set_ylabel('intensity')
        self.b.tick_params(axis='both', labelsize=7)
        self.c.set_title('PHD')
        self.c.tick_params(axis='both', labelsize=7)
        
        self.canvas = FigureCanvas(self.f) 
             
        
    # Navigtion toolbar stuff     
        
        toolbar = NavigationToolbar(self.canvas, self)
        main_box = Gtk.Box(spacing = 2, orientation = Gtk.Orientation.VERTICAL)  
        self.add(main_box)
        
    # button box
        vbutton_box = Gtk.HButtonBox()
        button1 = Gtk.Button('Raw count rate')
        button2 = Gtk.Button('Filter PHD')
        button3 = Gtk.Button('Fit 1D Gauss')
        
        vbutton_box.pack_start(button1,True,True, 0)
        vbutton_box.pack_start(button2,True, True, 0)
        vbutton_box.pack_start(button3,True, True, 0)
        
        
        
        
    # packing in main_box
        main_box.pack_start(self.canvas, True, True, 0)
        main_box.pack_start(vbutton_box,False,False,0)
        main_box.pack_start(toolbar, False, False, 0)
            
        
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
    
    # sends 2d data to my gui plotting funtion	
    	self.update_plot(scidata)
    
   
    ##### smashing 2D data into 1D to view orders as peaks ####
    	
    
    # filling in y with the sums of the rows of scidata
        y=[]
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
    	
    	#plt.figure(figsize=(7.5,8.4))
    #
    # using scipy.signal.find_peaks_cwt() to find centers of orders.  this required scipy version .11.0 or greater
    	peakind=signal.find_peaks_cwt(ychunk,np.arange(3,15))
    
    # plotting chunk of 1D data with lines through the centers of the orders to double check how the peak finder did

    	lines=xchunk[peakind]
    	revlines=lines[::-1]
    	self.update_ordersplot(ychunk,xchunk,lines)
    
    
    ### extraction of orders ###
    
    
    # find w, the widths of the orders (difference between peaks)

    	w=[]
    	for i in range(1,len(peakind)):
    		t=peakind[i]-peakind[i-1]
    		w.append(t)
    	
    # i have to add an extra w at the end of the array to make it the right size i (hopefully this is kosher)
    	maxw=max(w)-4
    	w.append(maxw)
    
    # placeing verticle lines in the visualization where the widths would stand
        #avew=sum(w)/len(w)
    	#wvlines=[]
    	#b=xchunk[peakind[0]]-avew/2
    	#wvlines.append(b)
    	#for i in np.arange(len(peakind)):
    		#f=xchunk[peakind[i]]+w[i]/2
    		#wvlines.append(f)
    	#plt.vlines(wvlines,0,2500,color='r',linestyle='--',label='w boundries')
    	#plt.legend(loc=1)
    
    	#plt.figure()
    
    ### making arrays of 1s ans 0s and extracting the 1d orders by matrix multiplication 
   	zeros=np.zeros((len(x),1))
    	index = range(0,len(w))
    	reindex = index[::-1]
    
    	oneDorders = {}
    	for i in reindex:
    		zeros1=np.copy(zeros)
    		zeros1[ 1024 - (np.sum(w[(i):18])) : 1024 - np.sum(w[(i+1):18]) ] = 1
    		twoD = scidata*zeros1
    			   		
            # making 2d orders in to 1d orders
                Y=[]
    		for j in range(0,len(scidata[0])):
        		t = np.sum(twoD[:,j])
        		Y.append(t)
            # placing 1d orders in dictionary called oneDorders
        	oneDorders['order'+str(i)]=Y
        	
     # sending plotting info to update_1dplot for gui (for now using just on order until cross coralation is added to script
    	x = np.linspace(0,len(scidata[0]), num = len(scidata[0]))
        odo=oneDorders['order16']
    	self.update_1dplot(odo,x)

    	
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
   

### fake PDH stuff ### (fake data for now)
	xphd=np.arange(10)
	yphd=(.2*xphd)**3
	self.update_PHDplot(xphd,yphd)
	
 ### preperation for guass fitting ###
 	#data = the 1d data between two mouse clicks or a box drawn
 	#mu = average inside that data area (center of peak) 	
 	#FWHM = sqrt(abs(sum((x-mu)**2*data)/sum(data)))
 	#max = data.max()
 	#gaussfit = lambda t : max*exp(-(t-mu)**2/(2*FWHM**2))
 	#plot(gaussfit(x))
 	
 	

    	
### plotting in gui ###
    def update_plot(self, scidata):
        self.plt= self.a.imshow(scidata, vmin = 0, vmax = 255)
        cbar=self.f.colorbar(self.plt,shrink=.84,pad=0.01)
        #cbar.self.a.tick_params(labelsize=5) 
        self.canvas.draw()
        
    def update_1dplot(self,odo,x):
    	self.plt=self.b.plot(x,odo)
        self.canvas.draw()
        
    def update_PHDplot(self,xphd,yphd):
        self.plt=self.c.plot(xphd,yphd)
        self.canvas.draw()
        
    def update_ordersplot(self,ychunk,xchunk,lines):
        self.plt=self.e.plot(ychunk,xchunk)
        self.e.hlines(lines,0,2000,color='purple',label='centers')
        self.canvas.draw()
        
 ##  preping for mouse interaction
    #def onclick(event):
    	#print event.button, event.x, event.y, event.xdata, event.ydata

        cid = aelf.f.canvas.mpl_connect('button_press_event', onclick)


        

    



if __name__ == '__main__':
    win = MyWindow()
    win.connect("delete-event", Gtk.main_quit)
    win.show_all()
    Gtk.main()

