import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import statistics
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
plt.style.use(astropy_mpl_style)

bias_data = []
for i in range(100,110):
    filename = f"bias/d{i}.fits"
    data = fits.getdata(filename, ext=0)
    bias_data.append(data)

master_bias = np.median(bias_data, axis=0)
plt.hist2d(master_bias[:, 0], master_bias[:, 1], bins = 25)
plt.show()

flat_data = []
for i in range(110,115):
    filename = f"flat/B/d{i}.fits"
    data = fits.getdata(filename, ext=0)
    flat_data.append(data)

master_flat = np.median(flat_data, axis=0)
plt.hist2d(master_flat[:, 0], master_flat[:, 1], bins = 25)
plt.show()

clean_flat = master_flat - master_bias
clean_flat_mean = np.mean(clean_flat)
norm_clean_flat = clean_flat / clean_flat_mean

plt.hist2d(norm_clean_flat[:, 0], norm_clean_flat[:, 1], bins = 25)
plt.show()

M53B_data = []
for i in range(171,176):
    filename = f"M53/B/d{i}.fits"
    data = fits.getdata(filename, ext=0)
    M53B_data.append(data)

exptimes = []
for i in range(171,176):
    filename = f"M53/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

persec_B = []
for i in range(0, 5):
    clean_data = M53B_data[i] - master_bias
    persec = clean_data[i] / exptimes[i]
    persec_B.append(persec)

master_M53B = np.median(persec_B, axis=0)

calibrated_M53B = master_M53B / norm_clean_flat
plt.figure("calibrated science")
plt.imshow(calibrated_M53B,cmap='gray', vmin=0,vmax=np.amax(calibrated_M53B))
plt.show()