# -*- coding: utf-8 -*-
"""
Code to fit a triple gaussian plus exponential to all mass data.
Includes the curve_fit method only.
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

# Create histogram of all mass data and return the bin entries and bin edges 
entries1, binedges1 = np.histogram(xmass, bins=int((np.max(xmass)-np.min(xmass))/0.005), density=True)

# Define the x data as bin centres
xdata = (binedges1[:-1]+binedges1[1:])/2
#print(xdata)
# Define the y data as bin entries
ydata = entries1
#print(ydata)
# Calculate binwdiths 
binwidth = binedges1[1]-binedges1[0]

# Calculate the mean and standard deviation to use as initial parameter guesses of the fit
mu_1 = np.mean(peak_range)
sigma_1 = np.std(peak_range)

# Define the gaussian part of the fit 
def gaussian(x, mu, sigma):
    return (1/(sigma*((2*np.pi)**(1/2))))*np.exp(-0.5*((x - mu)/sigma)**2)

# Define the exponetial part of the fit
def exponential(x, b):
    return np.exp(b*x)

# Combine three gaussians and the expoential to produce the complete fit function
# Iterate through each value in peak_1 list and calculate the fit value (f(x))
# Append each fit value to flist and return
def fit(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, sigma_3, b, N, M, P):
    flist = []
    a = 1/((1/b)*(np.exp(b*np.max(x))-np.exp(b*np.min(x)))) #+c*(np.max(x)-np.min(x)))
    for i in x:
        f = N*gaussian(i, mu_1, sigma_1) + M*gaussian(i, mu_2, sigma_2) + P*gaussian(i, mu_3, sigma_3) + (1-N-M-P)*(a*exponential(i, b))
        flist.append(f)
    return flist

# Define the function to calculate the sum of the negative logs 
# This value is minimised using the minimze method (not used here)
def neg_log(params, x, func):
    ll = []
    data = func(x, *params)
    for i in data:
        ll.append(np.log(i))   
    return -np.sum(ll)

# Return the optimal parameters covariance matrix using curve_fit 
popt, pcov = curve_fit(fit, xdata, ydata, p0=[mu_1, sigma_1, 10, 0.05, 10.3, 0.05, -0.3, 0.2, 0.2, 0.2])
print(popt)
# Caluculate parameter errors by finding the square roots of the diagonals
print(np.sqrt(np.diag((pcov))))


#res = minimize(neg_log, x0=[*popt], args=(xmass, fit), method='L-BFGS-B')
#print(res)
#popt2 = res.x

# Produce all plots
plt.title("Triple gaussian fitting of all invariant mass data")
plt.xlabel("Invariant mass of muon pairs (GeV/c^2)")
plt.ylabel("PDF")
plt.bar(xdata, ydata, width=binwidth, color='blue', label='Histogram')
plt.plot(xdata, fit(xdata, *popt), color='green', label='curve_fit')
#plt.plot(xdata, fit(xdata, *popt2), color='red', label='fit')
plt.legend(loc="upper left")

# Show plots
plt.show()
