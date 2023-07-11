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
plt.xlabel('Pixel')
plt.ylabel('Intensity [a.u.]')
plt.title('Low Resolution Mercury Spectrum')
plt.savefig('Low Res Hg.png')
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
axes=[]
xaxis = np.array(range(0,765))
for i in (range(len(a))):
    axis = (a[i] * xaxis) + b[i]
    axes.append(axis)

table1 = [393.3682, 394.4016, 396.1535, 396.8492, 404.5825, 406.3605, 407.1749, 407.7724, 410.1748, 
          413.2067, 414.3878, 416.7277, 420.2040, 422.6740, 423.5949, 425.0130, 425.0797, 425.4346, 
          426.0486, 427.1774, 432.5775, 434.0475, 438.3557, 440.4761, 441.5135, 452.8627, 455.4036,
          455.4036, 455.4036, 486.1342, 489.1502, 492.0514, 495.7613, 516.7327, 517.2698, 518.3619, 
          525.0216, 526.9550, 532.8051, 552.8418, 588.9973, 589.5940, 610.2727, 612.2226, 616.2180,
          630.2499, 656.2808]

# for i in range(0,8):
#     plt.figure(figsize=(10, 5))
#     plt.plot(axes[i],Sun_spectra[i])
#     plt.vlines(x=table1, ymin =np.min(Sun_spectra[i])-100, ymax=np.max(Sun_spectra[i])+100, color='red')
#     plt.xlim(((axes[i][0]), axes[i][-1]))
#     plt.xlabel('Wavelength')
#     plt.ylabel('Intensity')
#     plt.savefig(f'Sun Wavelength {i}')
#     plt.show()
print(393.3682-393.440)
# Sun400 = [393.440, Ca II,
#           396.861, Ca II,
#           404.549, Fe I,
#           406.364, Fe I,
#           410.200, H,
#           422.700, Ca I,
#           430.814, Fe I,
#           432.523, Fe I,
#           434.017, H,
#           438.395, Fe I,
#           440.424, Fe I]

# Sun480 = [438.252, Fe I,
#           452.716, Fe I,
#           466.643, Mg I,
#           486.142, H,
#           517.427, Mg I,
#           518.605, Mg I]

# Sun560 = [517.232, Mg I,
#           518.403, Mg I,
#           526.917, Fe I,
#           537.027, Fe I,
#           552.778, Mg I
#           589.070, Na I (D1),
#           ]
#wavelength, Intensity

# plt.plot(axes[i],Sun_spectra[i])
# plt.xlabel('Wavelength')
# plt.ylabel('Intensity')
# plt.show()

# A2_hdu = pyfits.open('A2.fits')
# G2_hdu = pyfits.open('G2_+0.0_Dwarf.fits')
# K2_hdu = pyfits.open('K2_+0.0_Dwarf.fits')
# F2_hdu = pyfits.open('F2_+0.0_Dwarf.fits')
# M2_hdu = pyfits.open('M2_+0.0_Dwarf.fits')
# A2_spec = A2_hdu[1].data
# G2_spec = G2_hdu[1].data
# K2_spec = K2_hdu[1].data
# F2_spec = F2_hdu[1].data
# M2_spec = M2_hdu[1].data
# print(len(A2_spec))
# A2x = []
# A2y = []
# G2x = []
# G2y = []
# K2x = []
# K2y = []
# F2x = []
# F2y = []
# M2x = []
# M2y = []
# for i in range(0,len(A2_spec)):
#     A2x.append((A2_spec[i][0]*100)+40)
#     A2y.append(A2_spec[i][1])
#     G2x.append((G2_spec[i][0]*100)+40)
#     G2y.append(G2_spec[i][1]+10)
#     K2x.append((K2_spec[i][0]*100)+40)
#     K2y.append(K2_spec[i][1]+15)
#     F2x.append((F2_spec[i][0]*100)+40)
#     F2y.append(F2_spec[i][1]+20)    
#     M2x.append((M2_spec[i][0]*100)+40)
#     M2y.append(M2_spec[i][1]+25)       
# plt.figure(figsize=(10, 7))
# plt.plot(A2x,A2y)
# plt.plot(G2x,G2y)
# plt.plot(K2x,K2y)
# plt.plot(F2x,F2y)
# plt.plot(M2x,M2y)
# plt.plot(axes[1],Sun_spectra[0])
# plt.plot(axes[1],Sun_spectra[5])
# plt.legend(['A2', 'G2', 'K2', 'F2', 'M2', 'The Sun 400nm', 'The Sun 600nm'])
# plt.show()
