
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# syntax:  $ python squashmatrix2.py 'filename.fits'

import sys
filename=sys.argv[1] # reads in the file name

import astropy.io.fits
hdulist = astropy.io.fits.open(filename) #this opens up the fits  file

targname = hdulist[0].header['targname']
print 'Target Name', targname

import pylab
import aplpy

scidata = hdulist['sci',2].data #this picks out the actual data from fits file, and turns it into numpy array

plt.subplot(221)
img=plt.imshow(scidata, vmin = 0, vmax = 255)


y=[]          # making an empty array to put the new values in


import numpy as np
for i in range(0,len(scidata[0])):  # filling in y with the sums of the rows of scidata 
    t=np.sum(scidata[i,:])
    y.append(t)


x=np.linspace(0,len(scidata[0]),num=len(scidata[0]))# making an x axis with same dementions as y

xrev=x[::-1] #reversing x for the sake of the visualization

plt.subplot(212)
plt.plot(y,xrev)
plt.show()



  
  


