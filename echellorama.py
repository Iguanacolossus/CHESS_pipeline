
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
import math
from scipy.optimize import curve_fit

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo as FigureCanvas
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar







class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self,title="Echelle Reduction GUI")
  
   ## setting up canvase for plotting ### 
   
        #self.box = Gtk.EventBox()
        
       
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
        main_box = Gtk.Box( orientation = Gtk.Orientation.VERTICAL)  
        self.add(main_box)
     
    # status bar
        self.statusbar = Gtk.Statusbar()
        context_id=self.statusbar.get_context_id("stat bar example")    
        
    # button box
        vbutton_box = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL)
        self.button1 = Gtk.Button('Raw count rate')
        self.button1.connect("clicked", self.on_button1_clicked, context_id)
        button2 = Gtk.Button('Filter PHD')
        
        self.button3 = Gtk.Button(label='Fit 1D Gauss')
        self.button3.connect("clicked", self.on_button3_clicked,  context_id)
        
        
        vbutton_box.pack_start(self.button1,True,True, 0)
        vbutton_box.pack_start(button2,True, True, 0)
        vbutton_box.pack_start(self.button3,True, True, 0)
        
     
        
    # packing in main_box
        main_box.pack_start(self.statusbar, False,False,0)
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
    	self.x = np.linspace(0,len(scidata[0]), num = len(scidata[0]))
        self.odo=oneDorders['order16']
        odo = self.odo[:]
    	self.update_1dplot(odo,x)

    	

   

### fake PDH stuff ### (fake data for now)
	PHDfake = '/home/rachel/codes/chesstest.fits'
	hdu = fits.open(PHDfake)
	PHD = hdu[1].data['PHD']
	self.update_PHDplot(PHD)
	hdu.close()
	
	
	
 ### guass fitting ###
        self.xdata=[]
    def gauss_fit(self,xdata):
        
        x = list(self.x)
        xg=[]     
        xg.append(int(xdata[0]))
        xg.append(int(xdata[1]))        
        xgauss = x[xg[0]:xg[1]]
        ygauss = self.odo[xg[0]:xg[1]]
        
        def gauss(xgauss,MAX,mu,sigma):
        	return MAX*np.exp(-(xgauss-mu)**2/(2.*sigma**2))
        avex = sum(xgauss)/len(xgauss) 
        guess = [1.,avex,1.]
        coeff, var_matrix = curve_fit(gauss,xgauss,ygauss,p0=guess) 
        fit = gauss(xgauss, *coeff)
      
 
 	plt.figure()
 	fitplot = plt.plot(xgauss,ygauss,xgauss,fit)
 	#fitplot.text()
 	print coeff
        plt.show(fitplot)
       
        
        
    	
    	
 	
 	
 	

    	
### plotting in gui ###
    def update_plot(self, scidata):
        self.plt= self.a.imshow(scidata, vmin = 0, vmax = 255)
        cbar=self.f.colorbar(self.plt,shrink=.84,pad=0.01)
        #cbar.self.a.tick_params(labelsize=5) 
        self.canvas.draw()
        
    def update_1dplot(self,odo,x):
    	self.plt=self.b.plot(x,self.odo)
    	
        self.canvas.draw()
        
    def update_PHDplot(self,PHD):
        self.plt=self.c.hist(PHD,bins=80,histtype='stepfilled')
        self.canvas.draw()
        
    def update_ordersplot(self,ychunk,xchunk,lines):
        self.plt=self.e.plot(ychunk,xchunk)
        self.e.hlines(lines,0,2000,color='purple',label='centers')
        self.canvas.draw()
        
      
  ## gauss fitting button 
    
          
    def on_button3_clicked(self, widget, data):
         self.statusbar.push(data,'Ready to fit.  Click on both sides of the emission feature you wish to fit')
         
         
    	 def onclick(event):
    	         
    	         self.xdata.append(event.xdata)
        	 #xdata = event.xdata
        	 #self.xdata.append(xdata)
        	 if len(self.xdata) == 2:
        	 	self.statusbar.push(data,'fitting...')
        	 	xdata=self.xdata
    	                self.gauss_fit(xdata)
        	 
      
         
               
   ###  mouse click event 	
    	 cid = self.canvas.mpl_connect('button_press_event', onclick)

    	 #xdata=self.xdata
    	 #self.gauss_fit(xdata)
    	
    	 
    	 
    def on_button1_clicked(self, widget,data):
    	self.statusbar.push(data,'Use mouse to draw box in 2d data area to see count rate')
    	
    	def onclick2(event):
                  print event.xdata, event.ydata
                  
        def offclick2(event):
                  print event.xdata, event.ydata
                  self.statusbar.push(data,'count rate = ')
         
        cid2 = self.canvas.mpl_connect('button_press_event',onclick2 )
        cid3 = self.canvas.mpl_connect('button_release_event',offclick2 )
        
    
    



if __name__ == '__main__':
    win = MyWindow()
    win.connect("delete-event", Gtk.main_quit)
    win.show_all()
    Gtk.main()

