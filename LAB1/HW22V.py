from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import statistics
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
#plt.style.use(astropy_mpl_style)

bias_data = []
for i in range(100, 110):
    filename = f"bias/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        bias_data.append(data)

#Median of all bias data lists
master_bias = np.median(bias_data, axis=0)

plt.hist2d(master_bias[:, 0], master_bias[:, 1], bins = 19)
plt.savefig('master_bias.pdf')
plt.show()

flat_data = []

for i in range(110, 115):
    filename = f"flat/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        flat_data.append(data)

master_flat = np.median(flat_data, axis=0)
plt.hist2d(master_flat[:, 0], master_flat[:, 1], bins = 25)
plt.savefig('master_bias.pdf')
plt.show()

##################
clean_flat = master_bias - master_flat
clean_mean = np.mean(clean_flat)
normalized_clean = clean_flat / clean_mean
##################

exptimes = []
for i in range(141, 146):
    filename = f"M67/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

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
        persec = cleandata / (exptimes[i - 141])
        persec_M67B.append(persec)

#Median of all science data lists
master_M67B = np.median(persec_M67B, axis=0)

calibrated_science = master_M67B / normalized_clean

#Graph of calibrated_science data
plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(calibrated_science, vmin=0, vmax=50)
plt.colorbar()
plt.savefig('calibrated_science.pdf')
plt.show()
