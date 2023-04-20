from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import statistics
#hi
bias_data = []
#For loop creates list of bias.fits data
for i in range(100, 110):
    filename = f"bias/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        bias_data.append(data)
master_bias = np.median(bias_data, axis=0)

flatB_data = []
for i in range(110, 115):
    filename = f"flat/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        flatB_data.append(data)
master_flatB = np.median(flatB_data, axis=0)

clean_flat = master_flatB - master_bias
clean_flat_1D = clean_flat.flatten()
clean_mean = statistics.mean(clean_flat_1D)
normalized_clean = clean_flat / clean_mean

raw_M53B = []
for i in range(171, 176):
    filename = f"M53/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        raw_M53B.append(data)

normalized_clean_non = np.where(normalized_clean == 0, 1e-10, normalized_clean)

exptimes = []
for i in range(171, 176):
    filename = f"M53/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

clean_M53B = []
persec_M53B = []
for i in range(0, 5):
    cleandata = raw_M53B[i] - master_bias[i]
    clean_M53B.append(cleandata)
    persec = cleandata[i] / exptimes[i]
    persec_M53B.append(persec)

master_M53B = np.median(persec_M53B, axis=0)
calibrated_M53B = master_M53B / normalized_clean_non

plt.imshow(calibrated_M53B, cmap = 'gray', vmin=0, vmax=np.max(calibrated_M53B))
plt.savefig('calibrated_M53B.pdf')
plt.show()
#print(calibrated_M53B)
