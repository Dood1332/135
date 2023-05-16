from astropy.io import fits as pyfits
from astropy import units as u

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from glob import glob

TestSpectrumLoResFiles = {"Sun" : "Low/Sun/X_560nm_Low_Sun.00000101.FIT",
                          "Hg"  : "Low/Hg/X_560nm_Low.00000100.FIT"}

TestSpectrumHiResFiles = {"Sun" : {4.0 : 'High/Sun/X_400nm_High_Sun.00000105.FIT',
                                   4.4 : 'High/Sun/X_440nm_High_Sun.00000107.FIT',
                                   5.2 : 'High/Sun/X_520nm_High_Sun.00000113.FIT',
                                   5.6 : 'High/Sun/X_560nm_High_Sun.00000115.FIT',
                                   6.0 : 'High/Sun/X_600nm_High_Sun.00000117.FIT',
                                   6.4 : 'High/Sun/X_640nm_High_Sun.00000119.FIT',
                                   6.8 : 'High/Sun/X_680nm_High_Sun.00000121.FIT'},
                          
                          "Hg"  : {4.0 : 'High/Hg/X_400nm_High_Hg.00000104.FIT',
                                   4.4 : 'High/Hg/X_440nm_High_Hg.00000106.FIT',
                                   5.2 : 'High/Hg/X_520nm_High_Hg.00000112.FIT',
                                   5.6 : 'High/Hg/X_560nm_High_Hg.00000114.FIT',
                                   6.0 : 'High/Hg/X_600nm_High_Hg.00000116.FIT',
                                   6.4 : 'High/Hg/X_640nm_High_Hg.00000118.FIT',
                                   6.8 : 'High/Hg/X_680nm_High_Hg.00000120.FIT'}}

Micrometer = [4.0, 4.4, 5.2, 5.6, 6.0, 6.4, 6.8]

hdu = pyfits.open(TestSpectrumLoResFiles['Hg'])

image_hdu = hdu[0]
spec_data = image_hdu.data
plt.imshow(spec_data)
plt.show()

#Get central 10 rows of the image
print(spec_data.shape)
img_cut = spec_data[250:260]
plt.imshow(img_cut)
plt.axis('off')
plt.show()

#Average rows to get accurate spectrum
spectrum = np.average(img_cut, axis=0)
plt.plot(spectrum)
plt.show()

def spectrum(img_data, row_min, row_max):
    return np.flip(np.average(img_data[row_min:row_max,:], axis=0))

#getting the mercury spectra
Hg_spectra = []
Sun_spectra = []
for microm in Micrometer:
    Hg_hdu = pyfits.open(TestSpectrumHiResFiles['Hg'][microm])
    Sun_hdu = pyfits.open(TestSpectrumHiResFiles['Sun'][microm])
    Hg_spec_data = Hg_hdu[0].data
    Sun_spec_data = Sun_hdu[0].data
    Hg_spectra.append(spectrum(Hg_spec_data, 250, 260))
    Sun_spectra.append(spectrum(Sun_spec_data, 250, 260))
    
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))
im1=ax1.plot(Hg_spectra[0])
im2=ax2.plot(Hg_spectra[1])
plt.suptitle("Mercury Spectra", fontsize=15)
plt.savefig("Mercury Spectra.png")
plt.show()
    
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))
im1=ax1.plot(Sun_spectra[0])
im2=ax2.plot(Sun_spectra[1])
plt.suptitle("Sun Spectra", fontsize=15)
plt.savefig("Sun Spectra.png")
plt.show()

lambda2 = 435.833
lambda1 = 404.656
p2 = 660
p1 = 368

a = (lambda2 - lambda1)/(p2 - p1)
print(a)