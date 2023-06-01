#conversion from decibels referenced to milliwatts. Values defined in lab manual.
#Power (mW) per RBW = 10^(dBM/10)
#convert to brightness distribution
#21cm line should match 1000k black body
#20.9cm matches with 100k 

from astropy.io import fits as pyfits
from astropy import units as u
import csv
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from glob import glob
from scipy.stats import norm
from scipy.optimize import curve_fit
import scipy.stats as stats

# with open('TRACE051.CSV', newline='') as f:
#     reader = csv.reader(f)
#     for row in reader:
#         print(row)

def csv_plot(csv_file):
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        next(reader, None)  # Skip header row

        x = []
        y = []

        for i, row in enumerate(reader):
            if i >= 11 and i <= 471:
                xaxis = float(row[1])
                yaxis = float(row[2])
                x.append(xaxis)
                y.append(yaxis)
        return x, y


csv_files = ['TRACE051.CSV', 'TRACE052.CSV', 'TRACE053.CSV', 'TRACE054.CSV', 'TRACE055.CSV',
             'TRACE056.CSV', 'TRACE057.CSV', 'TRACE058.CSV', 'TRACE059.CSV',]

x, y = csv_plot('TRACE051.CSV')
x = np.array(x)
y = np.array(y)
#print(x)
def func(x, a, b,c,d,e):
    return a*np.exp(-((x-b)**2)/(c))+d*x+e
#next line from curve_fit documentation. p0=(guess values)
popt, pcov = curve_fit(func, x, y, p0=(-103,1.42025*1e9,0.0050*1e9,0,7.085*1e-8))# bounds=(0, [.5, ,.5, 1., 0.5]))
plt.plot(x, func(x, *popt), 'r-',
         label='fit: a=%f, b=%f, c=%f, d=%f, e=%f' % tuple(popt))
#feeding my x values into my function with guess values
guess = func(x, a=-103, b=1.42025*1e9, c=4.0*1e9, d=0, e=7.085*1e-8)
plt.plot(x, y, 'b-', label='data')
#plt.plot(x, guess, 'g-', label='guess')

plt.figure(figsize = (10,5))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dBm)')

for csv_file in csv_files:
    x, y = csv_plot(csv_file)
    plt.plot(x, y)

plt.legend(['051', '052', '053', '054', '055', '056', '057', '058', '059'])
#plt.savefig('Plot.png')
plt.show()
