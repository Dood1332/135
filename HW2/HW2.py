from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import statistics

#NOTE: I had all the fits files in a folder called Fits to keep files organized

##### 1 #####
bias1 = fits.open('Fits/bias01.fits')
bias1data = bias1['PRIMARY'].data

#Graph of one of the bias.fits
plt.hist2d(bias1data[:, 0], bias1data[:, 1], bins = 35)
plt.savefig('bias2dhist.png')
plt.show()
bias1.close()
##############

##### 2a #####
counts, bins = np.histogram(bias1data[:, 0], bins = 35)
plt.bar(bins[:-1], counts, width=np.diff(bins), align='edge')
plt.xlabel('# of counts per pixel')
plt.ylabel('# of pixels associated with count')
plt.savefig('bias1dhist.png')
plt.show()
#############

##### 3 #####
bias_data = []
#For loop creates list of bias.fits data
for i in range(1, 6):
    filename = f"Fits/bias0{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        bias_data.append(data)

#Median of all bias data lists
master_bias = np.median([bias_data[0], bias_data[1], bias_data[2], bias_data[3], bias_data[4]], axis=0)

#Graph of master bias


plt.hist2d(master_bias[:, 0], master_bias[:, 1], bins = 20)
plt.savefig('master_bias.png')
plt.show()

#Same process as bias
dark_data = []
for i in range(1, 6):
    filename = f"Fits/dark0{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        dark_data.append(data)

master_dark = np.median([dark_data[0], dark_data[1], dark_data[2], dark_data[3], dark_data[4]], axis=0)

plt.hist2d(master_dark[:, 0], master_dark[:, 1], bins = 19)
plt.savefig('master_dark.png')
plt.show()

#Same process as dark and bias
flat_data = []
for i in range(1, 6):
    filename = f"Fits/flat0{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        flat_data.append(data)

master_flat = np.median([flat_data[0], flat_data[1], flat_data[2], flat_data[3], flat_data[4]], axis=0)

plt.hist2d(master_flat[:, 0], master_flat[:, 1], bins = 25)
plt.savefig('master_flat.png')
plt.show()
##################

##### 3b & 4 #####
clean_flat = master_flat - master_bias
clean_flat_1D = clean_flat.flatten()
clean_mean = statistics.mean(clean_flat_1D)
normalized_clean = clean_flat / clean_mean

plt.imshow(normalized_clean)
plt.show()
##################

##### 5 #####

#Pulls exposure times from each science.fits into a list. Used in science.fits loop.
exptimes = []
for i in range(1, 6):
    filename = f"Fits/science0{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

#Creates data list for science.fits, subtracts master_bias from it and into a list, and then divides the clean list by exposure times
science_data = []
clean_science_data = []
persec_science = []

for i in range(1, 6):
    filename = f"Fits/science0{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        science_data.append(data)
        cleandata = data - master_bias
        clean_science_data.append(cleandata)
        persec = cleandata / exptimes[i-1]
        persec_science.append(persec)

#Median of all science data lists
master_science = np.median([persec_science[0], persec_science[1], persec_science[2], persec_science[3], persec_science[4]], axis=0)

#Cleans up zeroes for when master_science is divided by normalized_clean
nonzeroclean = np.where(normalized_clean == 0, 1e-10, normalized_clean)
nonzero_master = np.where(master_science == 0, 1e-10, master_science)
calibrated_science = master_science / normalized_clean

#Graph of calibrated_science data
plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(calibrated_science,interpolation='nearest', cmap='gray', vmin=0, vmax=50)
plt.colorbar()
plt.savefig('calibrated_science.png')
plt.show()