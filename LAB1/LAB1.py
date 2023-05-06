from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sep
from matplotlib import rcParams
from matplotlib.patches import Ellipse


flat_filters = [(110, 115, "B"), (115, 120, "V"), (120, 125, "R")]
calibration_filters = [(126, 131, "cali/B"), (131, 136, "cali/V"), (136, 141, "cali/R")]
science_filters = [(171, 176, "M53/B"), (176, 181, "M53/V"), (181, 186, "M53/R"),
                   (141, 146, "M67/B"), (146, 151, "M67/V"), (151, 156, "M67/R")]


def calibrate_science_data(flat_filters, science_filters):
    bias_data = []
    for i in range(100, 110):
        filename = f"bias/d{i}.fits"
        with fits.open(filename) as hdulist:
            data = hdulist['PRIMARY'].data
            bias_data.append(data)

    #Median of all bias data lists
    master_bias = np.median(bias_data, axis=0)

    flat_data = []
    for flat_filter in flat_filters:
        for i in range(flat_filter[0], flat_filter[1]):
            filename = f"flat/{flat_filter[2]}/d{i}.fits"
            with fits.open(filename) as hdulist:
                data = hdulist['PRIMARY'].data
                flat_data.append(data)

    master_flat = np.median(flat_data, axis=0)

    clean_flat = master_bias - master_flat
    clean_mean = np.mean(clean_flat)
    normalized_clean = clean_flat / clean_mean

    science_data = []
    for science_filter in science_filters:
        exptimes = []
        science_filter_data = []
        for i in range(science_filter[0], science_filter[1]):
            filename = f"{science_filter[2]}/d{i}.fits"
            with fits.open(filename) as hdulist:
                exptime = hdulist[0].header['EXPTIME']
                exptimes.append(exptime)
                data = hdulist['PRIMARY'].data
                science_filter_data.append(data)

        clean_science_data = []
        persec_science_data = []
        for i in range(len(science_filter_data)):
            cleandata = science_filter_data[i] - master_bias
            clean_science_data.append(cleandata)
            persec = cleandata / (exptimes[i])
            persec_science_data.append(persec)

        master_science_data = np.median(persec_science_data, axis=0)

        calibrated_science_data = master_science_data / normalized_clean
        
        science_data.append(calibrated_science_data)

    return science_data


calibrated_science_data = calibrate_science_data(flat_filters, science_filters)
calibrated_calibration_data = calibrate_science_data(flat_filters, calibration_filters)

for i in range(0, 2):
        calibrated_science_data[i][calibrated_science_data[i] < 0] = 0
        calibrated_calibration_data[i][calibrated_calibration_data[i] < 0] = 0

plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(calibrated_science_data[3], interpolation='nearest', vmin=0, vmax=50)
plt.colorbar()
plt.savefig('calibrated_calibratedM53B.png')
plt.show()

#Making the HR Diagram

def bg_subtraction(data):
    bkg = sep.Background(data)
    
    bkg_image = bkg.back()
    bkg_rms = bkg.rms()

    data_sub = data - bkg

    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)

    
    print('Objects detected:',len(objects))

    fig, ax = plt.subplots()
    m, s = np.mean(data_sub), np.std(data_sub)
    im = ax.imshow(data_sub, interpolation='nearest', cmap='gray', vmin=0, vmax=50)
    
    #limit ellipse size
    min_size = 1
    max_size = 25

    for i in range(len(objects)):
        a = objects['a'][i]
        b = objects['b'][i]
        if a > min_size and a < max_size and b > min_size and b < max_size:
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
    plt.show()

bg_subtraction(calibrated_science_data[3])


def mag_atm_ext(filters, m_atm_ext):
    calibration_airmass = []
    for i in range(filters[0], filters[1]):
        filename = f"{filters[2]}/d{i}.fits"
        with fits.open(filename) as hdulist:
            data = hdulist[0].header['AIRMASS']
            calibration_airmass.append(data)
    mag_airmass = np.median(calibration_airmass) * m_atm_ext
    return mag_airmass

def flux(data):
    zp = 25.0
    gain = 1.0
    bkg = sep.Background(data)
    data_sub = data - bkg
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
    flux, flux_err, flag = sep.sum_circle(data_sub, objects['x'], objects['y'], 50, err=bkg.globalrms, gain=gain)
    return flux

F_sciB = flux(calibrated_science_data[0])
F_sciV = flux(calibrated_science_data[1])
F_sciR = flux(calibrated_science_data[2])

F_caliB = flux(calibrated_calibration_data[0])
F_caliV = flux(calibrated_calibration_data[1])
F_caliR = flux(calibrated_calibration_data[2])

print(len(F_sciB),len(F_sciV),len(F_sciR))


# Atm sci
m_atm_extB = .234 #B 478.5
m_atm_extV = .460 #V 386.2
m_atm_extR = .094 #R 710.0

m_atm_ext_sciB = mag_atm_ext(science_filters[0], m_atm_extB)
m_atm_ext_sciV = mag_atm_ext(science_filters[1], m_atm_extV)
m_atm_ext_sciR = mag_atm_ext(science_filters[2], m_atm_extR)

#print(m_atm_ext_sciB, m_atm_ext_sciV, m_atm_ext_sciR)

#Dust sci
m_dust_ext_sciB = 0.092
m_dust_ext_sciV = 0.083
m_dust_ext_sciR = 0.060

#Atm cali
m_atm_ext_caliB = mag_atm_ext(calibration_filters[0], m_atm_extB)
m_atm_ext_caliV = mag_atm_ext(calibration_filters[1], m_atm_extV)
m_atm_ext_caliR = mag_atm_ext(calibration_filters[2], m_atm_extR)
#print(m_atm_ext_caliB, m_atm_ext_caliV, m_atm_ext_caliR)

#Dust cali
m_dust_ext_caliB = 0.173
m_dust_ext_caliV = 0.148
m_dust_ext_caliR = 0.107

#Instrument Magnitudes
m_I_B = -2.5 * np.log(F_caliB) - m_atm_ext_caliB - m_dust_ext_caliB
m_I_V = -2.5 * np.log(F_caliV) - m_atm_ext_caliV - m_dust_ext_caliV
m_I_R = -2.5 * np.log(F_caliR) - m_atm_ext_caliR - m_dust_ext_caliR

m_I_sciB = -2.5 * np.log(F_sciB) - m_atm_ext_sciB - m_dust_ext_sciB
m_I_sciV = -2.5 * np.log(F_sciV) - m_atm_ext_sciV - m_dust_ext_sciV
m_I_sciR = -2.5 * np.log(F_sciR) - m_atm_ext_sciR - m_dust_ext_sciR

#Flux Magnitude 
m_B = -2.5 * np.log(F_caliB)
m_V = -2.5 * np.log(F_caliV)
m_R = -2.5 * np.log(F_caliR)

#Zero point
m_zpB = 13.056 - m_I_B
m_zpV = 13.327 - m_I_V
m_zpR = 13.456 - m_I_R

print(len(m_zpB), len(m_zpV), len(m_zpR))

# # plt.scatter()
# # plt.show()