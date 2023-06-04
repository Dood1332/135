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
for csv_file in csv_files:
    x, y = csv_plot(csv_file)
    x = np.array(x)
    y = np.array(y)
    a = np.max(y)
    print(a)
    max_index = np.argmax(y)
    b = x[max_index]

    def func(x, a, b,c,d,e):
        return a*np.exp(-((x-b)**2)/(c))+d*x+e
    #next line from curve_fit documentation. p0=(guess values)
    popt, pcov = curve_fit(func, x, y, p0=(a, b, 0.0050*1e9, 0 , 7.085*1e-8))# bounds=(0, [.5, ,.5, 1., 0.5]))
    plt.figure(figsize = (10,5))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (dBm)')
    plt.plot(x, func(x, *popt), 'r-',
            label='fit: a=%f, b=%f, c=%f, d=%f, e=%f' % tuple(popt))
    #feeding my x values into my function with guess values
    plt.plot(x, y, 'b-', label='data')
    plt.legend(['Fitted Curve', f'{csv_file[0:8]}'])
    plt.savefig(f'{csv_file}.png')
    plt.show()

# plt.figure(figsize = (10,5))
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Amplitude (dBm)')

# for csv_file in csv_files:
#     x, y = csv_plot(csv_file)
#     plt.plot(x, y)

# plt.legend(['051', '052', '053', '054', '055', '056', '057', '058', '059'])
# #plt.savefig('Plot.png')
# plt.show()
