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
import sep
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from scipy.signal import argrelextrema
from astropy import units as u

fits = ['fc_150.336238+55.898786_sdss(dr7)_g.fits',
        'fc_150.336238+55.898786_sdss(dr7)_i.fits',
        'fc_150.336238+55.898786_sdss(dr7)_r.fits',
        'fc_150.336238+55.898786_sdss(dr7)_u.fits',
        'fc_150.336238+55.898786_sdss(dr7)_z.fits']

#Open data and subtracts the background
for fit in fits:
    hdu = pyfits.open(fit)

    data = hdu[0].data

    def bg_subtraction(data):
        #Cropping the image so we only get the twin quasar. Comment the line below to get the full image.
        data = data[700:740 , 1942:1987]
        data = data.astype(float)
        bkg = sep.Background(data)
        
        bkg_image = bkg.back()
        bkg_rms = bkg.rms()

        data_sub = data - bkg
        objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)

        
        #print('Objects detected:',len(objects))

        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray', vmin=0, vmax=75)

        return data_sub

    #Plotting the image
    data = bg_subtraction(data)
    plt.imshow(data, interpolation='nearest', cmap='gray', vmin=0, vmax=75)
    plt.show()

    #Taking the average of the cropped image into plottable 1D graph
    avg = np.average(data, axis=0)
    plt.plot(avg)
    plt.show()

    #Guessing a gaussian for the graph
    def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
        return (A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) +
                A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2)))

    x = np.arange(len(avg))
    initial_guess = [375, 16, 3, 303, 30.5, 3]#a1, b1, c1, a2, b2, c2
    params, _ = curve_fit(double_gaussian, x, avg, p0=initial_guess)

    A1_fit, mu1_fit, sigma1_fit, A2_fit, mu2_fit, sigma2_fit = params

    x_fit = np.linspace(0, len(avg) - 1, 100)
    y_fit = double_gaussian(x_fit, A1_fit, mu1_fit, sigma1_fit, A2_fit, mu2_fit, sigma2_fit)

    maxima_indices = argrelextrema(y_fit, np.greater)
    x_maxima = x_fit[maxima_indices]

    #Finding values for the true guess/fit
    a1 = ((y_fit[maxima_indices])[::1][0])
    a2 = ((y_fit[maxima_indices])[::1][1])

    sorted_indices = np.argsort(y_fit[maxima_indices])[::-1]
    x_peaks = x_maxima[sorted_indices][:2]

    b1 = x_peaks[0]
    b2 = x_peaks[1]

    c = 3

    #Using true values for Gaussian fit
    new_guess = [a1, b1, c, a2, b2, c]
    params, _ = curve_fit(double_gaussian, x, avg, p0=initial_guess)

    A1_fit, mu1_fit, sigma1_fit, A2_fit, mu2_fit, sigma2_fit = params

    x_fit = np.linspace(0, len(avg) - 1, 100)
    y_fit = double_gaussian(x_fit, A1_fit, mu1_fit, sigma1_fit, A2_fit, mu2_fit, sigma2_fit)

    maxima_indices = argrelextrema(y_fit, np.greater)
    x_maxima = x_fit[maxima_indices]
    sorted_indices = np.argsort(y_fit[maxima_indices])[::-1]
    x_peaks = x_maxima[sorted_indices][:2]

    #Difference between peaks
    diff = x_peaks[1] - x_peaks[0]

    #Plotting the data and fitted curve
    print('Difference between peaks:', diff, 'pixels')
    plt.plot(avg)
    plt.plot(x_fit, y_fit, 'r')
    plt.legend(['Data', 'Fited Curve'])
    plt.xlabel('Pixel Number')
    plt.ylabel('Counts')
    plt.title(f'Twin Quasar {fit[34:35]}')
    plt.savefig(f'Twin_Quasar_{fit[34:35]}.png')
    plt.show()

    #Converting pixel to arcsec
    sdss_pixelscale = u.pixel_scale(0.4*u.arcsec/u.pixel)
    scale = (diff*u.pixel).to(u.arcsec, sdss_pixelscale)
    print(scale)

