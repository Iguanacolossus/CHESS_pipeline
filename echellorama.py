
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
plt.xlabel('pixels(x)')
plt.ylabel('pixels(y)')
img=plt.imshow(scidata, vmin = 0, vmax = 255)


y=[]          # making an empty array to put the new values in


import numpy as np
for i in range(0,len(scidata[0])):  # filling in y with the sums of the rows of scidata 
    t=np.sum(scidata[i,:])
    y.append(t)


x=np.linspace(0,len(scidata[0]),num=len(scidata[0]))# making an x axis with same dementions as y

xrev=x[::-1] #reversing x for the sake of the visualization

chunk=[0,500]# the orders seem to blend around lyman alpha so i have to select the biggest chunk I could  use to get a delta y between the orders. i wil have to make this selectable on the GUI some how
bgcutoff=475 # cuting off back ground here to get zeros and peaks. i will ahve to make this selectable some how on  the GUI 
minv=chunk[0]
maxv=chunk[1]
plt.subplot(212)
plt.plot(y,xrev)
plt.ylabel('pixels(y)')
plt.xlabel('Intensity')
plt.vlines(bgcutoff,[0],1000,label='background cutoff')
plt.hlines(chunk,[-1000],[5000],color='r')
plt.vlines([-1000,5000],minv,maxv,color='r')
plt.show()

#cutting out chunk of data
xchunk= x[max(x)-chunk[1]:max(x)]
index1=int(max(x)-chunk[1])
index2=int(max(x))
ychunk=y[index1:index2]
#reversing x for the sake of the plot
xrevchunk=xchunk[::-1]

 #creating a new baseline
for i in range(0,xchunk.size):
  if ychunk[i] < bgcutoff:
     ychunk[i] = 0


#plt.plot(ychunk,xrevchunk)
#plt.show()

bounds=[]
a=ychunk
# find the bounds of a order by indice in the array
for i in range(1,len(a)):
     if a[i-1] == 0 and a[i] > 0:
          print i , "begin hump index"
          bounds.append(i)
     if a[i-1] > 0 and a[i] == 0:
          print i-1, "end hump index"
          bounds.append(i-1)

# now 'bounds' is a list of the indices of the array were a #order begins and ends
  

# defining a function to group the indecies of 'bounds' to # sets of (begin order, end order)
def group(lst,n):
   for i in range(0, len(lst), n):
       val=lst[i:i+n]
       if len(val) ==n:
          yield tuple(val)

grpbounds= list(group(bounds,2)) # 'bounds' but now grouped

centers=[]
# finding max with in the bouds of the orders and reading them #in to 'centers' list
for i in range(len(grpbounds)):
   centers.append(max(a[grpbounds[i][0]:grpbounds[i][1]+1]))
   
 
#print a.index(centers)



  
  


