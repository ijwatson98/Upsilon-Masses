# -*- coding: utf-8 -*-
""" 
Authors: Emma Carrie and Isaac Watson

Code to plot 2d histograms of different observables against invariant mass for the data set

"""

import  numpy  as  np
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt

f  =  open("ups-15-small.bin","r")
datalist  =  np.fromfile(f,dtype=np.float32)
#  number  of  events
n_event  =  len(datalist)/6
x_data  =  np.split(datalist,n_event)
print(x_data[0])
#  make  list  of  invariant  mass  of  events

xmass  =  []
for  i  in  range(0, int(n_event)):
        xmass.append(x_data[i][0])
        if  i  <  10:
                print(xmass[i])
                
# Create obseravable data arrays 
trans_mom_pair = []
for  i  in  range(0, int(n_event)):
        trans_mom_pair.append(x_data[i][1])
      #  if  i  <  10:
     #           print(trans_mom_pair[i])              
rapidity = []
for  i  in  range(0, int(n_event)):
        rapidity.append(x_data[i][2])
       # if  i  <  10:
                #print(rapidity[i])
mom_pair = []
for  i  in  range(0, int(n_event)):
        mom_pair.append(x_data[i][3])
       # if  i  <  10:
               # print(mom_pair[i])
trans_mom_1 = []
for  i  in  range(0, int(n_event)):
        trans_mom_1.append(x_data[i][4])
        if  i  <  10:
                print(trans_mom_1[i])
trans_mom_2 = []
for  i  in  range(0, int(n_event)):
        trans_mom_2.append(x_data[i][5])
    #    if  i  <  10:
           #     print(trans_mom_2[i])

# Plot 2D histograms adjusted inputs where neccesary              
entries1, xbinedges1, ybinedges1, patches1 = pylab.hist2d(xmass, rapidity, bins=100, density=True)
pylab.title("2D histogram of rapidity of muon pair and invariant mass of pair")
plt.xlabel("Invariant mass of muon pair (GeV/c^2)")
plt.ylabel("Rapidity")
plt.ylim(2, 10)
pylab.show

#entries2, binedges2, patches2 = pylab.hist(trans_mom_pair, bins=int(30/0.01), range=(0, 30))
#pylab.title("Number of Pairs of Particles of different transverse momentum")
#pylab.xlabel("Momentum (GeV/C)")
#pylab.ylabel("Number of Pairs")
#pylab.show
