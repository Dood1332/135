import math
import numpy as np
from matplotlib import pyplot as plt
import astropy.io.fits as fits

import matplotlib.pyplot as plt

from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename

#HW 2.1
#Plot a 2D histogram of one of the bias fits files ->save the image
biats_fits_image_data_2D = fits.getdata("bias04.fits")

#code to display the image 2D histogram
print("displaying bias04.fits")
plt.style.use(astropy_mpl_style)
plt.figure("bias04.fits")
plt.imshow(biats_fits_image_data_2D)
plt.colorbar()
plt.savefig('hw21-bias-2D-hist.pdf')
plt.show()

"""
#2.2 Plot the count distribution of the bias image.
Plotting 1D taking 5+ minutes, hence commented out for now
print("plotting histogram for bias 1D data")
plt.figure("bias04 1D histogram")
plt.hist(bias_fits_image_data_1D, bins='auto')
plt.savefig('hw21-bias-1D-hist.pdf')
plt.show()
"""

#2.3 Create three new 2D histograms: median of the biases fits, median of the dark fits, median of the flat fits ->
#plot and save 2D histograms of these 3 images (master_bias, master_dark, master_flat)
# **each pixel of your new 2D histogram is the median value of the same pixel of your previous 2D histograms**

#create 3 list variables for file names
image_list_bias = ["bias01.fits","bias02.fits","bias03.fits","bias04.fits","bias05.fits"]
image_list_dark = ["dark01.fits","dark02.fits","dark03.fits","dark04.fits","dark05.fits"]
image_list_flat = ["flat01.fits","flat02.fits","flat03.fits","flat04.fits","flat05.fits"]

print("compute medians master bias, master data, master flat ")
#read image data from bias files into array of bias images
bias_image_data_array = [fits.getdata(image) for image in image_list_bias]

#create 2D histogram for the median of the biases fit
master_bias = np.median(bias_image_data_array, axis=0)

#read image data from dark fit files into array of dark images
dark_image_data_array = [fits.getdata(image) for image in image_list_dark]

#create 2D histogram for the median of the dark fits
master_dark= np.median(dark_image_data_array, axis=0)

#read image data from fits data files into array of fits images
flat_image_data_array = [fits.getdata(image) for image in image_list_flat]

#create 2D histogram for the median of the flat fits
master_flat= np.median(flat_image_data_array, axis=0)

#plot 2D histogram for master bias
print("displaying master bias")
plt.style.use(astropy_mpl_style)
plt.figure("Master Bias 2D")
plt.imshow(master_bias)
plt.colorbar()
plt.savefig('hw21-bias-median-2D.pdf')
plt.show()

#plot 2D histogram for master dark
print("displaying master dark")
plt.figure("Master Dark 2D")
plt.imshow(master_dark)
plt.colorbar()
plt.savefig('hw21-dark-median-2D.pdf')
plt.show()

#plot 2D histogram for master flat
print("displaying master flat")
plt.figure("Master Flat 2D")
plt.imshow(master_flat)
plt.colorbar()
plt.savefig('hw21-flat-median-2D.pdf')
plt.show()

#HW2 3b)
#Subtract the master_bias from the master_flat (clean_flat)
clean_flat = master_flat - master_bias
#HW2 4)
#Normalize the clean_flat to its mean (which means dividing each pixel by the pixel array mean value) ->plot and save
import statistics
print('normalize clean fit')
clean_flat_1D = clean_flat.flatten() #first make 1D out of 2D clean-flat data
clean_flat_mean = statistics.mean(clean_flat_1D) #compute mean
print('clean_flat_mean=',clean_flat_mean)
clean_flat_norm= clean_flat/clean_flat_mean #normalize the clean flat

#Display normalized clean flat
plt.figure("Norm Clean Flat")
plt.imshow(clean_flat_norm)
plt.colorbar()
plt.savefig('Hw2-clean-flat-2D.pdf')
plt.show()

#create array to hold the science fit file names
image_list_science = ["science01.fits","science02.fits","science03.fits","science04.fits","science05.fits"]

#read science fit data files into an array of 2D images
science_image_data_array = [fits.getdata(image) for image in image_list_science]

#use for loop to substract master bias from each off science 2D histogram data
print("reading science data files")
for i in range(5):
   science_image_data_array[i] = science_image_data_array[i] - master_bias
   #print(fits.getheader(image_list_science[i]).persec_science)
   #divide each pixel by exposure time for counts/sec
   #I couldn't use the exptime dynamically but above print statement showed 30 sec for all science images
   science_image_data_array[i] = science_image_data_array[i]/30.0

print("calibrate science data")
#create 2D histogram using median of all science fit data files
master_science = np.median(science_image_data_array, axis=0)
norm_science = master_science/clean_flat_norm #normalize by dividing with clean flat norm
plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(norm_science,vmin=0,vmax=np.amax(norm_science))
plt.colorbar()
plt.savefig('Hw2-norm_science.pdf')
plt.show()