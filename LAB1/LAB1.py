from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sep
from matplotlib import rcParams
from matplotlib.patches import Ellipse

m_dust_ext = 0.083

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


flat_filters = [(110, 115, "B"), (115, 120, "V"), (120, 125, "R")]
calibration_filters = [(126, 131, "cali/B"), (131, 136, "cali/V"), (136, 141, "cali/R")]
science_filters = [(171, 176, "M53/B"), (176, 181, "M53/V"), (181, 186, "M53/R"), (141, 146, "M67/B"), (146, 151, "M67/V"), (151, 156, "M67/R")]
calibrated_science_data = calibrate_science_data(flat_filters, science_filters)
calibrated_calibration_data = calibrate_science_data(flat_filters, calibration_filters)

# plt.figure("calibrated science")
# ax = plt.axes()
# ax.set_facecolor("red")
# plt.imshow(calibrated_science_data[0], vmin=0, vmax=50)
# plt.colorbar()
# plt.savefig('calibrated_current_index.pdf')
# plt.show()

plt.figure("calibrated science")
ax = plt.axes()
ax.set_facecolor("red")
plt.imshow(calibrated_calibration_data[0], vmin=0, vmax=50)
plt.colorbar()
plt.savefig('calibrated_calibratedB.pdf')
plt.show()

def bg_subtraction(data):
    bkg = sep.Background(data)
    
    bkg_image = bkg.back()
    bkg_rms = bkg.rms()
    
    data_sub = data - bkg

    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
    

    fig, ax = plt.subplots()
    m, s = np.mean(data_sub), np.std(data_sub)
    im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                    vmin=m-s, vmax=m+s, origin='lower')
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.show()

#bg_subtraction(calibrated_science_data[0])

def HR(data):
    zp = 25.0
    gain = 1.0
    bkg = sep.Background(data)
    data_sub = data - bkg
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
    flux, flux_err, flag = sep.sum_circle(data_sub, objects['x'], objects['y'], 3, err=bkg.globalrms, gain=gain)
    return flux
print(HR(calibrated_calibration_data[0]))

# cali_dataB = []
# for i in range(126, 131):
#     filename = f"cali/B/d{i}.fits"
#     with fits.open(filename) as hdulist:
#         data = hdulist['PRIMARY'].data
#         cali_dataB.append(data)
# master_caliB = np.median(cali_dataB)

# print(master_caliB)