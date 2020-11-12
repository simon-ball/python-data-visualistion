import json
import numpy as np
import scipy.optimize



d = r"C:\Users\simoba\Documents\_personal\Enduring Documents\Physics\2013 - PhD, Durham\Thesis Backup\Thesis_backup_2018-05-06_printed\Root\thesis\Figures\Single_channel_2\Rabi_Osc_SP\Figure Data\\"


"""
Data to store:
    main:
        time
        S recovery
        p recovery
        model time
        s model
        p model
    inset:
        recovery_time
        microwave
        blue
        recovery_b
        recovery_c
        recovery_c
"""

data = {"main": {},
        "inset": {}}

background = 0.000193/6.
time_adjust = 408.0/6.


data_file_name_S = d+"17_2016-07-12_2600-3000.txt"
data_file_name_P = d+"17_2016-07-12_3100-3300.txt"
data_file_name_all=d+"17_2016-07-12_2600-3300.txt"

S = np.loadtxt(data_file_name_S) - (background*4.)
P = np.loadtxt(data_file_name_P) - (background*4.)
all_ = np.loadtxt(data_file_name_all) - (background*7.)

norm = max([max(S), max(P)])

S = S/norm
P = P/norm
all_ = all_/norm
axis = np.linspace(450,0,46) / time_adjust

data["main"]["time"] = list(axis)
data["main"]["s"] = list(S)
data["main"]["p"] = list(P)
data["main"]["all"] = list(all_)



fit_axis = np.linspace(7, 0, 250)

def fit_cos_n(x, floor, amplitude, rate, phase, tau, n):
    # Fit the function cos^2n in the range
    #floor + (amplitude * (np.cos((rate * x) + phase)**(2*n))  )
    angle = np.cos( (rate*x) + phase)
    raised = abs(angle)**(2*n)          # Necessary to cope with scaling negative irregular numbers
    scaled = floor + (amplitude * np.exp(-1*x/tau) * raised)    
    return scaled

# S
guess_S = (0, 1, 1.57, 0, 100, 1)
bounds_S = [ (0, 0.8, 1.57, -0.1, 10, 0.8), (0.03, 1.2, 1.571, 0.1, 200, 1.5) ]
popt, pcov = scipy.optimize.curve_fit(fit_cos_n, axis, S, guess_S, bounds=bounds_S)
fit_S = fit_cos_n(fit_axis, *popt)

# P
guess_P = (0, 0.8, 1.57, 1.57, 100, 1)
bounds_P = [ (0, 0.4, 1.57, 1.57, 10, 1), (0.1, 1.2, 1.571, 1.571, 200, 1.5) ]
popt, pcov = scipy.optimize.curve_fit(fit_cos_n, axis, P, guess_P, bounds=bounds_P)
fit_P = fit_cos_n(fit_axis, *popt)



data["main"]["model_t"] = list(fit_axis)
data["main"]["model_s"] = list(fit_S)
data["main"]["model_p"] = list(fit_P)


### Create profiles for photon arrival times
# Blue
def sigmoid(x,floor=5, amp=750, x0=2750, rate=42.7):
    # 10-90 rise time is 4.35x the rate, not entirely sure why. 
    # Observed 10-90 rise time is 186ns -> rate=42.7
    return floor + ( amp / (1 + np.exp(-1 * (x-x0) / rate) ) )

blue_x = np.linspace(2500,3550, 1500)
blue_y = sigmoid(blue_x)

# MW
def square(x, floor=5, amp = 750, x0 = 3140, x1 = 3390):
    return (sigmoid(x,0,1,x0,rate=2) + sigmoid(x, 0, 1, x1, rate=-2) - 1) * amp
mw_y = square(blue_x)



data["inset"]["blue_x"] = list(blue_x)
data["inset"]["blue_y"] = list(blue_y)
data["inset"]["mw_y"] = list(mw_y)

data_file_SP_arr = d+"56_2016-07-12_arrivals.txt"
data_file_SP_bins = d+"56_2016-07-12_bins.txt"
data_file_P_arr = d+"55_2016-07-12_arrivals.txt"
data_file_P_bins = d+"55_2016-07-12_bins.txt"
data_file_S_arr = d+"54_2016-07-12_arrivals.txt"
data_file_S_bins = d+"54_2016-07-12_bins.txt"

data["inset"]["s"] = list(np.loadtxt(data_file_S_arr))
data["inset"]["s_t"] = list(np.loadtxt(data_file_S_bins)[:-1])
data["inset"]["p"] = list(np.loadtxt(data_file_P_arr))
data["inset"]["p_t"] = list(np.loadtxt(data_file_P_bins)[:-1])
data["inset"]["sp"] = list(np.loadtxt(data_file_SP_arr))
data["inset"]["sp_t"] = list(np.loadtxt(data_file_SP_bins)[:-1])





file = r"../thesis_fig_7-3.json"

with open(file, "w+") as f:
    json.dump(data, f)