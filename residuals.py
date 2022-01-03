"""
Code to compare the single and double gaussian fit residuals.
Includes the curve_fit method only.

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

# Combine 2 gaussians and exponential to produce the complete fit function
# Iterate through each value in peak_1 list and calculate the fit value (f(x))
# Append each fit value to flist and return
def fit(x, mu, sigma, mu2, sigma2, b, N, M):
    flist = []
    a = 1/((1/b)*(np.exp(b*np.max(x))-np.exp(b*np.min(x)))) 
    for i in x:
        f = N*gaussian(i, mu, sigma) + M*gaussian(i, mu2, sigma2) + (1-N-M)*(a*exponential(i, b))
        flist.append(f)
    return flist

def fit2(x, mu, sigma, b, N):
    f2list = []
    a = 1/((1/b)*(np.exp(b*np.max(x))-np.exp(b*np.min(x)))) #+c*(np.max(x)-np.min(x)))
    for i in x:
        f2 = N*gaussian(i, mu, sigma) + (1-N)*(a*exponential(i, b))
        f2list.append(f2)
    return f2list

# Define the fucntion to calculate the sum of the negative logs 
# This value is minimised using the minimze method (not used here)
def neg_log(params, x, func):
    ll = []
    data = func(x, *params)
    for i in data:
        ll.append(np.log(i))    
    return -np.sum(ll)

# Return the optimal parameters covariance matrix for doubel gaussian 
popt1, pcov1 = curve_fit(fit, xdata, ydata, p0=[mu, sigma, mu, sigma, -0.3, 0.2, 0.2])
print(popt1)
# Caluculate parameter errors by finding the square roots of the diagonals
print(np.sqrt(np.diag((pcov1))))

# Return the optimal parameters covariance matrix for single gaussian 
popt2, pcov2 = curve_fit(fit2, xdata, ydata, p0=[mu, sigma, -0.3, 0.2])
print(popt2)
# Caluculate parameter errors by finding the square roots of the diagonals
print(np.sqrt(np.diag((pcov2))))

# Obtain residuals for double gaussian
residuals1 = ydata - fit(xdata, *popt1)
# Obtain residuals for single gaussian
residuals2 = ydata - fit2(xdata, *popt2)

# Produce all plots
#plt.title("Residual data for double gaussian")
plt.title("Residual data for single gaussian")
plt.xlabel("Invariant mass of muon pairs (GeV/c^2)")
plt.ylabel("Residual")
# Plot resiudals 
#plt.plot(xdata, residuals1, color='red', label='fit')
plt.plot(xdata, residuals2, color='green', label='fit')

plt.show()