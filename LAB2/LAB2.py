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

lambdas = [404.656, 435.833, 
           404.656, 435.833, 
           491.607, 
           579.066, 
           546.074, 579.066,
           576.960, 579.066,
           660.428, 667.728,
           696.543, 706.722]
p       = [368, 660, 
           73, 364, 
           508, 
           650, 
           281, 591,
           191, 211,
           682, 753,
           583, 680]

a = []
b = []

i = 0
# for i in range(0, len(lambdas), 2):
while i < len(lambdas):
    counter = 0
    if i == 4:
        a_value = a[-1]
        counter = 1
    elif i == 5:
        a_value = (lambdas[i+1] - lambdas[i+2]) / (p[i+1] - p[i+2])
        counter = 1
    else:
        a_value = (lambdas[i] - lambdas[i+1]) / (p[i] - p[i+1])
        counter = 2
    b_value = lambdas[i] - (p[i] * a_value)
    a.append(a_value)
    b.append(b_value)
    i += counter

print(len(a))
print(len(b))

axes=[]
xaxis = np.array(range(0,765))
for i in (range(len(a))):
    axis = (a[i] * xaxis) + b[i]
    axes.append(axis)

table1 = [393.3682, 394.4016, 396.1535, 396.8492, 404.5825, 406.3605, 407.1749, 407.7724, 410.1748, 
          413.2067, 414.3878, 416.7277, 420.2040, 422.6740, 423.5949, 425.0130, 425.0797, 425.4346, 
          426.0486, 427.1774, 432.5775, 434.0475, 438.3557, 440.4761, 441.5135, 452.8627, 455.4036]

for i in range(0,8):
    plt.figure(figsize=(10, 5))
    plt.plot(axes[i],Sun_spectra[i])
    plt.xlabel('Wavelength')
    plt.ylabel('Intensity')
    plt.savefig(f'Sun Wavelength{i}')
    plt.show()