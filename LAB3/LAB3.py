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
from scipy.interpolate import UnivariateSpline

#report for extra lab due next thursday

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

    #Unit conversion
    y += 3  #Add 3dB for cable attenuation 
    y -= 37  #Subtract by 37dB. Power to LNA
    y = (10**(y/10))/10  #dBm to W
    y = y*(0.001/30)  #W/RBW to W/Hz
    y *= 2  #Factor of 2


    a = np.max(y)
    #print(a)
    max_index = np.argmax(y)
    b = x[max_index]

    def make_norm_dist(x, mean, sd):
        return 1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x - mean)**2/(2*sd**2))
    ytemp = make_norm_dist(x, b, 1)
    # create a spline of x and blue-np.max(blue)/2 
    spline = UnivariateSpline(x, ytemp-np.max(ytemp)/2, s=0)
    r1, r2 = spline.roots()
    c = r2 - r1

    def func(x, a, b,c,d,e):
        return a*np.exp(-((x-b)**2)/(2*(c**2)))+d*x+e
    #next line from curve_fit documentation. p0=(guess values)
    popt, pcov = curve_fit(func, x, y, p0=(a, b, c, 0 , 7.085*1e-8))# bounds=(0, [.5, ,.5, 1., 0.5]))
    
    def sub(x, d, e):
        sub = []
        s = d*x+e
        sub.append(s)
        return sub

    # plt.figure(figsize = (10,5))
    # plt.xlabel('Frequency (Hz)')
    # plt.ylabel('Amplitude (dBm)')
    # plt.plot(x, func(x, *popt), 'r-',
    #         label='fit: a=%f, b=%f, c=%f, d=%f, e=%f' % tuple(popt))
    # #feeding my x values into my function with guess values
    # plt.plot(x, y, 'b-', label='data')
    # plt.legend(['Fitted Curve', f'{csv_file[0:8]}'])
    # #plt.savefig(f'{csv_file}.png')
    # plt.show()

    y = y - sub(x, popt[-2], popt[-1])
    a = np.max(y[0,:])
    #print(a)
    max_index = np.argmax(y[0,:])
    b = x[max_index]

    def make_norm_dist(x, mean, sd):
        return 1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x - mean)**2/(2*sd**2))
    ytemp = make_norm_dist(x, b, 1)
    # create a spline of x and blue-np.max(blue)/2 
    spline = UnivariateSpline(x, ytemp-np.max(ytemp)/2, s=0)
    r1, r2 = spline.roots()
    c = r2 - r1

    #v = ((1.4204*10**9) - x)/ (1.4204*10**9)
    v = (3.0*10**8)*((-(1.4204*10**9)/x)+1)/1000

    # def func(x, a, b,c,d,e):
    #     return a*np.exp(-((x-b)**2)/(2*(c**2)))+d*x+e
    # popt, pcov = curve_fit(func, x, y[0,:], p0=(a, b, c, 0 , 3.085*1e-16))
    # #print(len(y[0,:]))
    # plt.figure(figsize = (10,5))
    # plt.xlabel('Frequency (Hz)')
    # plt.ylabel('Watts per Hertz (w/Hz)')
    # plt.plot(x, func(x, *popt), 'r-',
    #         label='fit: a=%f, b=%f, c=%f, d=%f, e=%f' % tuple(popt))
    # plt.plot(x, y[0,:], 'b-', label='data')
    # plt.legend(['Fitted Curve', f'{csv_file[0:8]}'])
    # #plt.savefig(f'{csv_file}.png')
    # plt.show()

    plt.figure(figsize = (10,5))
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Watts per Hertz (w/Hz)')
    plt.grid()
    plt.vlines(x=0,ymin=np.min(y[0,:]),ymax=np.max(y[0,:]), color='r')
    plt.plot(v, y[0,:])
    plt.legend(['21 cm line', f'{csv_file[0:8]}'])
    plt.savefig(f'VELOCITY{csv_file[0:8]}.png')
    plt.show()

    #calculating the hydrogen column density   a
    dv = (x[2]-x[1])
    beta = 0.5
    k = 1.38*(10**(-23))
    Tb = (y[0,:])/(2*beta*k)

    N = ((1.8224*(10**18))*dv*np.sum(Tb))
    print(N)


# plt.figure(figsize = (10,5))
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Amplitude (dBm)')

# for csv_file in csv_files:
#     x, y = csv_plot(csv_file)
#     plt.plot(x, y)

# plt.legend(['051', '052', '053', '054', '055', '056', '057', '058', '059'])
# #plt.savefig('Plot.png')
# plt.show()
