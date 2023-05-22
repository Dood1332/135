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
                                   4.8 : 'High/Sun/X_480nm_High_Sun.00000111.FIT',
                                   5.2 : 'High/Sun/X_520nm_High_Sun.00000113.FIT',
                                   5.6 : 'High/Sun/X_560nm_High_Sun.00000115.FIT',
                                   6.0 : 'High/Sun/X_600nm_High_Sun.00000117.FIT',
                                   6.4 : 'High/Sun/X_640nm_High_Sun.00000119.FIT',
                                   6.8 : 'High/Sun/X_680nm_High_Sun.00000121.FIT'},
                          
                          "Hg"  : {4.0 : 'High/Hg/X_400nm_High_Hg.00000104.FIT',
                                   4.4 : 'High/Hg/X_440nm_High_Hg.00000106.FIT',
                                   4.8 : 'High/Hg/X_480nm_High_Hg.00000110.FIT',
                                   5.2 : 'High/Hg/X_520nm_High_Hg.00000112.FIT',
                                   5.6 : 'High/Hg/X_560nm_High_Hg.00000114.FIT',
                                   6.0 : 'High/Hg/X_600nm_High_Hg.00000116.FIT',
                                   6.4 : 'High/Hg/X_640nm_High_Hg.00000118.FIT',
                                   6.8 : 'High/Hg/X_680nm_High_Hg.00000120.FIT'}}

Micrometer = [4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8]

hdu = pyfits.open(TestSpectrumLoResFiles['Hg'])

image_hdu = hdu[0]
spec_data = image_hdu.data
# plt.imshow(spec_data)
# plt.show()

#Get central 10 rows of the image
#print(spec_data.shape)
img_cut = spec_data[250:260]
# plt.imshow(img_cut)
# plt.axis('off')
# plt.show()

#Average rows to get accurate spectrum
spectrum = np.average(img_cut, axis=0)
# plt.plot(spectrum)
# plt.xlabel('Pixel')
# plt.ylabel('Intensity [a.u.]')
# plt.title('Low Resolution Mercury Spectrum')
# plt.savefig('Low Res Hg.png')
# plt.show()

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
    

# for i in range(0,7,2):
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))
#     plt.xlabel('Pixel')
#     plt.ylabel('Intensity [a.u.]')
#     im1=ax1.plot(Hg_spectra[i])
#     im2=ax2.plot(Hg_spectra[i+1])
#     plt.suptitle("High Resolution Mercury Spectra", fontsize=15)
#     plt.savefig(f"Mercury Spectra {i}, {i+1}.png")
#     plt.show()

# for i in range(0,7,2):
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))
#     plt.xlabel('Pixel')
#     plt.ylabel('Intensity [a.u.]')
#     im1=ax1.plot(Sun_spectra[i])
#     im2=ax2.plot(Sun_spectra[i+1])
#     plt.suptitle("Sun Spectra", fontsize=15)
#     plt.savefig(f"Sun Spectra{i}, {i+1}.png")
#     plt.show()

#Hg wavelengths and pixel #

lambdas = [404.656, 435.833, 546.074, 576.960, 579.066, 696.543, 706.722]
p       = [368, 660, 650, 191, 211, 583, 680]

hdulist = pyfits.open('High/Hg/X_600nm_High_Hg.00000116.FIT')
header = hdulist[0].data
print(len(header[0]))
#print(header[191])
#368, 660
a400 = (lambdas[0] - lambdas[1])/(p[0] - p[1])
b400 = lambdas[0] - (p[0]*a400)
print(a400)
print(b400)


xaxis = np.array(range(0,765))
xaxis = (a400 * xaxis) + b400
#print((xaxis))

table1 = [393.3682, 394.4016, 396.1535, 396.8492, 404.5825, 406.3605, 407.1749, 407.7724, 410.1748, 413.2067, 414.3878, 416.7277, 420.2040, 422.6740, 423.5949, 425.0130, 425.0797, 425.4346, 426.0486, 427.1774, 432.5775, 434.0475, 438.3557, 440.4761, 441.5135, 452.8627, 455.4036]

plt.figure(figsize=(10, 5))
plt.plot(xaxis,Sun_spectra[0])
plt.xlabel('Wavelength')
plt.ylabel('Intensity')
plt.show()

plt.scatter(range(0,len(table1)), table1)
plt.show()