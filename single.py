# -*- coding: utf-8 -*-
"""
Code to fit a gaussian plus exponential to the Y(1S) peak.
Includes the curve_fit and minimize method.
Optimal parameter values and their errors are found and used to produce plots. 

Authors: Isaac Watson and Emma Carrie
"""

import  numpy  as  np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize

#  Import  data
#  xmass  =  np.loadtxt(sys.argv[1])
f  =  open("ups-15-small.bin","r")
datalist  =  np.fromfile(f,dtype=np.float32)
#  Number  of  events
n_event  =  len(datalist)/6
x_data  =  np.split(datalist,n_event)
#print(x_data[0])

# Make  lists  of  invariant  mass  of  events
xmass  =  []
peak_1 = []
peak_range = []

# Append all invariant mass items to xmass list
for  i  in  range(0, int(n_event)):
        xmass.append(x_data[i][0])

# Append invariant mass items around the Y(1S) peak to peak_1
#
for i in xmass:
    if i >9.3 and i<9.6:
        peak_1.append(i)
    if i >9.35 and i<9.56:
        peak_range.append(i)

#print(peak_1)

#print(np.max(xmass))
#print(np.min(xmass))                

# Create histogram of all Y(1S) mass data and return the bin entries and bin edges 
entries1, binedges1 = np.histogram(peak_1, bins=int((np.max(xmass)-np.min(xmass))/0.005), range=(9.3, 9.6), density=True)

# Define the x data as bin centres
xdata = (binedges1[:-1]+binedges1[1:])/2
#print(xdata)
# Define the y data as bin entries
ydata = entries1
#print(ydata)
# Calculate binwdiths 
binwidth = binedges1[1]-binedges1[0]

# Calculate the mean and standard deviation to use as initial parameter guesses of the fit
mu = np.mean(peak_range)
sigma = np.std(peak_range)

# Define the gaussian part of the fit 
def gaussian(x, mu, sigma):
    return (1/(sigma*((2*np.pi)**(1/2))))*np.exp(-0.5*((x - mu)/sigma)**2)

# Define the exponetial part of the fit
def exponential(x, b):
    return np.exp(b*x)

# Combine the gaussian and expoential to produce the complete fit function
# Iterate through each value in peak_1 list and calculate the fit value (f(x))
# Append each fit value to flist and return
def fit(x, mu, sigma, b, N):
    flist = []
    a = 1/((1/b)*(np.exp(b*np.max(x))-np.exp(b*np.min(x)))) #+c*(np.max(x)-np.min(x)))
    for i in x:
        f = N*gaussian(i, mu, sigma) + (1-N)*(a*exponential(i, b))
        flist.append(f)
    return flist

# Define the fucntion to calculate the sum of the negative logs 
# This value is minimised using the minimze method
def neg_log(params, x, func):
    ll = []
    data = func(x, *params)
    for i in data:
        ll.append(np.log(i))    
    return -np.sum(ll)

# Return the optimal parameters covariance matrix using curve_fit 
popt, pcov = curve_fit(fit, xdata, ydata, p0=[mu, sigma, -0.3, 0.2])
print(popt)
# Caluculate parameter errors by finding the square roots of the diagonals
print(np.sqrt(np.diag((pcov))))

# Return the optimal data by minimisng the NLL
# Use the parameters found in curve_fit as initial guesses
res = minimize(neg_log, x0=[*popt], args=(peak_1, fit), method='L-BFGS-B')
print(res)
# Obtain the optimal parameters from the data
popt2 = res.x
# Caluclate the errors on the parameters through the square roots of the diagonals of the inverse hessian
print(np.sqrt(np.diag(res.hess_inv.todense())))

# Create a list containig just the background contribution to the fit
# Can be used to observe the gaussian in isolation 
exp_data = []
a = 1/((1/popt[2])*(np.exp(popt[2]*np.max(xdata))-np.exp(popt[2]*np.min(xdata))))
for i in xdata:
    exp_data.append((1-popt[3])*a*exponential(i, popt[2]))

# Produce all plots
plt.title("curve_fit vs minimize fitting of Y(1S) invariant mass data")
plt.xlabel("Invariant mass of muon pairs (GeV/c^2)")
plt.ylabel("PDF")
plt.legend(loc="upper left")
plt.bar(xdata, ydata, width=binwidth, color='blue', label='Histogram')
plt.plot(xdata, fit(xdata, *popt), color='green', label='curve_fit')
plt.plot(xdata, fit(xdata, *popt2), color='red', label='minimize')
# To view the gaussian in isolation remove # below:
#plt.plot(xdata, np.array(fit(xdata, *popt2))-np.array(exp_data), color='yellow')
plt.legend(loc="upper left")

# Show plots
plt.show()
