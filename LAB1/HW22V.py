from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import statistics
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
plt.style.use(astropy_mpl_style)

bias_data = []
for i in range(100, 110):
    filename = f"bias/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        bias_data.append(data)

#Median of all bias data lists
master_bias = np.median([bias_data[0], bias_data[1], bias_data[2], bias_data[3], bias_data[4]], axis=0)

plt.hist2d(master_bias[:, 0], master_bias[:, 1], bins = 25)
plt.savefig('master_bias.pdf')
plt.show()

flat_data = []
flat_exp = []
persec_flat = []
for i in range(110, 115):
    filename = f"flat/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        flat_data.append(data)
        exptime = hdulist[0].header['EXPTIME']
        flat_exp.append(exptime)

for i in range(0, 5):
    persec = flat_data[i] / flat_exp[i]
    persec_flat.append(persec)

master_flat = np.median(flat_data, axis=0)
plt.hist2d(master_flat[:, 0], master_flat[:, 1], bins = 25)
plt.savefig('master_bias.pdf')
plt.show()
# plt.imshow(master_flat)
# plt.savefig('master_flat.pdf')
# plt.show()
##################

##### 3b & 4 #####
clean_flat = master_bias - master_flat
clean_flat_1D = clean_flat.flatten()
clean_mean = statistics.mean(clean_flat_1D)
normalized_clean = clean_flat / clean_mean

plt.imshow(normalized_clean)
plt.show()
##################

##### 5 #####

#Pulls exposure times from each science.fits into a list. Used in science.fits loop.
exptimes = []
for i in range(141, 146):
    filename = f"M67/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

#Creates data list for science.fits, subtracts master_bias from it and into a list, and then divides the clean list by exposure times
M67B_data = []
clean_M67B_data = []
persec_M67B = []

for i in range(141, 146):
    filename = f"M67/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        M67B_data.append(data)
        cleandata = data - master_bias
        clean_M67B_data.append(cleandata)
        persec = cleandata / (exptimes[i - 141] * 900)
        persec_M67B.append(persec)

# for i in range(0, 5):
#     cleandata = M53_data[i] - master_bias
#     clean_M53_data.append(cleandata)
#     persec = cleandata / exptimes[i]
#     persec_M53B.append(persec)

#Median of all science data lists
master_M67B = np.median(persec_M67B, axis=0)

#Cleans up zeroes for when master_science is divided by normalized_clean
nonzeroclean = np.where(normalized_clean == 0, 1e-10, normalized_clean)
# master_M53B = np.where(master_M53B == 0, 1e-10, master_M53B)
calibrated_science = master_M67B / normalized_clean

#Graph of calibrated_science data
# plt.imshow(calibrated_science, cmap='gray', vmin=0, vmax=np.max(calibrated_science))
plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(calibrated_science,cmap='gray', vmin=0,vmax=np.max(calibrated_science))
plt.colorbar()
plt.savefig('calibrated_science.pdf')
plt.show()
