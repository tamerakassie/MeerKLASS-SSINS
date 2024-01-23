   #Code for creating the masked array and perform the sky-subtraction.

"""

Connect to katcal_py3 (public) to run script on Jupyter-notebook/Slurm.

"""

#Imports
import numpy as np
import matplotlib.pylab as plt
import warnings
import time
import pickle
import sys
import katcali
import katcali.io as kio
import katcali.label_dump as kl
import katcali.diode as kd

def good_ant(fname):

    """
    Fuction to retrieve list of good antennas.
    """

    data=kio.load_data(fname)
    bad_ants=kio.check_ants(fname)

    ants_good = []
    for i in np.array(kio.ant_list(data)):
        if i not in bad_ants:
            ants_good.append(i)
    else:
        print (str(i) + ' is bad')

    return ants_good
"""
Define a function that will access the data block and return visibility data, flags, noise diodes vector. 
"""
def visData(fname, ant, pol='h'):
    """
    Parameters:
    ----------
    Fname: Path to observation block.

    Returns:
    -------
    vis, flags, noise_diodes vector (nd)

    Note : Current function looks at one pol and one dish.
    """

    data = kio.load_data(fname)
    target, c0, band_ants, flux_model = kio.check_ants(fname)
    ants_good = good_ant(fname)
    data.select(ants=ant, pol=pol)
    recv = ant + pol
    corr_id = kio.cal_corr_id(data, recv)

    assert(recv == data.corr_products[corr_id][0])
    assert(recv == data.corr_products[corr_id][1])

    print("Correlation ID:", corr_id, "Receiver:", recv)

    # Load visibilities and flags
    vis, flags = kio.call_vis(fname, recv)
    print("Shape of vis:", vis.copy().shape)
    vis_backup = vis.copy()

    ra, dec, az, el = kio.load_coordinates(data)
    ang_deg = kio.load_ang_deg(ra, dec, c0)
    ch_ref = 800
    timestamps, freqs = kio.load_tf(data)
    dp_tt, dp_ss, dp_f, dp_w, dp_t, dp_s, dp_slew, dp_stop = kl.cal_dp_label(data, flags, ant, pol, ch_ref, ang_deg)


    nd_on_time, nd_cycle, nd_set = kd.cal_nd_basic_para(fname)
    nd_on_edge, nd_off_edge = kd.cal_nd_edges(timestamps, nd_set, nd_cycle, nd_on_time)
    nd_ratio, nd_0, nd_1x = kd.cal_nd_ratio(timestamps, nd_on_time, nd_on_edge, data.dump_period)


    nd_t0, nd_t1x, nd_s0, nd_s1x, nd_t0_ca, nd_t0_cb, nd_t1x_ca, nd_t1x_cb = kl.cal_label_intersec(dp_tt, dp_ss, nd_0, nd_1x)

    return vis, flags,nd_s0


def MaskedArrayVisibilityFlags(vis, flags, nd_s0):
    """
    Parameters:
    ----------
    visibility, flags and nd_s0 from the visData Fuction.

    Returns:
    --------
    Visibility-Flags Masked Array
    """
    #vis, flags, nd_s0 = visData(fname)
    data0 = vis.copy


    nd_flags = np.ones_like(vis, dtype=bool)          # Empty mask, where all values are set to true. True is flagged data
    nd_flags[nd_s0, :] = False                        # Set the data with noise diodes removed to False so that this data is not flagged as bad data.  This is the scan only data.
    other_flags = flags                               # All other flags from the visData function. Boolean value is True.
    allflags =  np.logical_or(nd_flags, other_flags)  # Apply logical operator or to combine the noise diode flags, and visilibity data at a specific stage flags. 

    data_masked = np.ma.masked_array(vis, mask=allflags, fill_value = np.nan)

    return data_masked

def SkySubtraction(data_masked):

    """
    Function Returns differencing of the  visibility masked array.
    """

    vis_ss =data_masked[1:,:] - data_masked[0:-1,:]
    visSS=vis_ss.filled()

    return visSS




# This script will take an observation block and return the Sky-Subtractions, combine the Sky-subtractions as a stacked array, by looping over the antennas in the array. And averages over the antennas to form the Incoherent Noise Spectrum.
fname='1630519596'

ants_good = good_ant(fname)
vis_products = {}
for i in ants_good[0:len(ants_good)]:
    vis, flags, nd_s0 = visData(fname, i)
    vis_products[i] = (vis, flags, nd_s0)

masked_array = {}
for dish, (vis, flags,nd_s0) in vis_products.items():
    masked_data = MaskedArrayVisibilityFlags(vis, flags, nd_s0)
    masked_array[dish] = masked_data
visSS_results = {}
for ant_value, masked_data in masked_array.items():
    # Applying the function to create SkySubtraction
    visSS_data = SkySubtraction(masked_data)
    visSS_results[ant_value] = visSS_data

stacked_visSS = np.stack(list(visSS_results.values()), axis=0)
file_path = 'stacked_visSS60Dishes1.npy'

# Save the NumPy array as a pickle file
np.save(file_path, stacked_visSS)


print(f"Sky Subtraction of Dishes saved as pickle file: {pickle_file_path}")

                                                                                                               154,0-1       Bot


