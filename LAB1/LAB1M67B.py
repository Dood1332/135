from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

bias_data = []
#For loop creates list of bias.fits data
for i in range(100, 110):
    filename = f"bias/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist['PRIMARY'].data
        bias_data.append(data)

#Median of all bias data lists
master_bias = np.median(bias_data, axis=0)

#Graph of master bias
plt.hist2d(master_bias[:, 0], master_bias[:, 1], bins=25)
plt.savefig('master_bias.pdf')
plt.show()

flatB_data = []
for i in range(110, 115):
    filename = f"flat/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        flatB_data.append(data)

master_flatB = np.median(flatB_data, axis=0)

raw_M67B = []
for i in range(141, 146):
    filename = f"M67/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        raw_M67B.append(data)
master_M67B = np.median(raw_M67B, axis=0)

clean_flat = master_flatB - master_bias
clean_mean = clean_flat.mean()
normalized_clean = clean_flat / clean_mean
normalized_clean = np.where(normalized_clean == 0, 1e-10, normalized_clean)

#plt.hist2d(normalized_clean[:, 0], normalized_clean[:, 1], bins=30)
#plt.show()

exptimes = []
for i in range(141, 146):
    filename = f"M67/B/d{i}.fits"
    with fits.open(filename) as hdulist:
        exptime = hdulist[0].header['EXPTIME']
        exptimes.append(exptime)

clean_M67B = []
persec_M67B = []
for i in range(0, 5):
    cleandata = raw_M67B[i] - master_bias[i]
    clean_M67B.append(cleandata)
    persec = cleandata[i] / exptimes[i]
    persec_M67B.append(persec)

master_M67B = np.median(persec_M67B, axis=0)
calibrated_M67B = master_M67B / normalized_clean

plt.imshow(calibrated_M67B, cmap = 'gray', vmin=0, vmax=np.max(calibrated_M67B))
#plt.savefig('calibrated_M53B.pdf')
plt.show()
#print(calibrated_M53B)
