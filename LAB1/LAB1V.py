from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

bias = fits.open('bias/d100.fits')
biasdata = bias['PRIMARY'].data

B = fits.open('M53/B/d171.fits')
Bdata = B['PRIMARY'].data

flat = fits.open('flat/B/d111.fits')
flatB = flat['PRIMARY'].data

plt.hist2d(flatB[:, 0], flatB[:, 1], bins = 25)
plt.show()

clean = biasdata - flatB
mean = np.mean(flatB)
clean_mean = clean / mean
normalized_clean = np.where(clean_mean == 0, 1e-10, clean_mean)
clean_B = Bdata - biasdata

exptime = B[0].header['EXPTIME']

persecB = clean_B / exptime

calibrated = persecB / normalized_clean

plt.imshow(calibrated, cmap = 'gray', vmin=0, vmax=np.max(calibrated))
#plt.savefig('calibrated_M53B.pdf')
plt.show()