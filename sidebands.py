# -*- coding: utf-8 -*-
"""
Code to calculate the yield of Y(1S) in the data set using the sideband method. 

Authors: Isaac Watson and Emma Carrie
"""

import pylab
import numpy as np

f  =  open("ups-15-small.bin","r")
datalist  =  np.fromfile(f,dtype=np.float32)
#  number  of  events
n_event  =  len(datalist)/6
x_data  =  np.split(datalist,n_event)
#print(x_data[0])
#  make  list  of  invariant  mass  of  events
xmass  =  []
peak_1 = []
peak_range = []

for  i  in  range(0, int(n_event)):
        xmass.append(x_data[i][0])

#print(mass)

# Make a histogram and store the results as arrays
entries, binedges, patches = pylab.hist(xmass, bins=int((np.max(xmass)-np.min(xmass))/0.005), range=(np.min(xmass), np.max(xmass)))


#Define peak value as max count. Define lower and upper edges of the peak
peak1 = binedges[np.argmax(entries)]
peak1_lower = peak1 - 0.15
peak1_upper = peak1 + 0.15

#Find counts in bands
counts1 = []
lower_side_counts = []
upper_side_counts = []

# Append data to band lists
for m in xmass:
     if m>(float(peak1_lower)) and m<(float(peak1_upper)):
          counts1.append(m)
     elif m>(float(peak1_lower - 0.3)) and m<(float(peak1_lower - 0.15)):
          lower_side_counts.append(m)
     elif m>(float(peak1_upper + 0.15)) and m<(float(peak1_upper + 0.3)):
          upper_side_counts.append(m)
  
# Show central band and sidebands on histogram           
entries4, binedges4, patches4 = pylab.hist(lower_side_counts, bins=int((np.max(xmass)-np.min(xmass))/0.005), range=(np.min(xmass), np.max(xmass)))
entries5, binedges5, patches5 = pylab.hist(upper_side_counts, bins=int((np.max(xmass)-np.min(xmass))/0.005), range=(np.min(xmass), np.max(xmass)))
entries6, binedges6, patches6 = pylab.hist(counts1, bins=int((np.max(xmass)-np.min(xmass))/0.005), range=(np.min(xmass), np.max(xmass)))
          
pylab.title("Number of Particles of Different Masses")
pylab.xlabel("Mass (GeV/C^2)")
pylab.ylabel("Number of Particles")

# Show graph
pylab.show()