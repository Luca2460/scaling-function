""" The script to process the dataset generated by the gen.py script
"""

from types import MethodType
import numpy as np
from scipy import optimize
import h5py
import matplotlib.pyplot as plt
from math import sqrt
import json

# dataset = "dataset.hdf5"
# dataset0 = "dataset0.hdf5" 

# #f = h5py.File(dataset, "r")
# # f0 = h5py.File(dataset0, "r") # To plot at H=0 and be used with retrieveData2() below

def MandSigmaFromSampleDict(d, err=False):
    Ms = []
    sigmas = []
    
    for sample in d.values():
        # AVERAGES OVER SAMPLES
        Ms.append(np.mean(sample[:, 8]))
        sigmas.append(np.std(sample[:, 8]))

    errs = None
    if err:
        errs = np.std(Ms) / sqrt(len(Ms))
        return np.mean(Ms), np.mean(sigmas), errs
    return np.mean(Ms), np.mean(sigmas)

def retrieveData():
    """Returns H, T (vectors) and M, sigma (rank 2 arrays)"""
    Hs = [f[Hname].attrs["H"] for Hname in f]
    H0group = f[list(f.keys())[0]]
    Ts = [H0group[Tname].attrs["T"] for Tname in H0group]

    MsAndSigmas = [[list(MandSigmaFromSampleDict(Tgroup)) for Tgroup in Hgroup.values()]
                   for Hgroup in f.values()]

    MsAndSigmas = np.array(MsAndSigmas)

    Ms = MsAndSigmas[:, :, 0]
    sigmas = MsAndSigmas[:, :, 1]
    
    return np.array(Hs), np.array(Ts), Ms, sigmas

def retrieveFromSampleDict(d):
    Ms = []
    Es = []
    sigmasE = []
    
    for sample in d.values():
        Ms.append(np.mean(sample[:, 8]))
        Es.append(np.mean(sample[:, 5]))
        sigmasE.append(np.std(sample[:, 5]))

    return np.mean(Ms), np.mean(Es), np.mean(sigmasE)
    

# def retrieveData2(f=f0):
#     """Returns H, T (vectors) and M, sigma (rank 2 arrays), and the error on M"""
#     H0group = f[list(f.keys())[0]]
#     Ts = [H0group[Tname].attrs["T"] for Tname in H0group]

#     MsAndSigmas = [[list(retrieveFromSampleDict(Tgroup)) for Tgroup in Hgroup.values()]
#                    for Hgroup in f.values()]

#     MsAndSigmas = np.array(MsAndSigmas)

#     Ms = MsAndSigmas[:, :, 0]
#     Es = MsAndSigmas[:, :, 1]
#     sigmasE = MsAndSigmas[:, :, 2]
    
#     return np.array(Ts), Ms, Es, sigmasE


def rescale(T, ki, H):
    Tc = 0.67

    delta = 4.9
    gamma = 1.39
    beta = 0.365

    Treduce = T / Tc - 1
    
    return Treduce ** (gamma + beta) / H, ki / H ** (1 / delta - 1)

def scaling():
    ki_fluct = 15 ** 3 * sigmas * sigmas / Ts

    ki_diff = (Ms[1::2] - Ms[::2]) / 0.1
    H_mean = (Hs[::2] + Hs[1::2]) / 2

    plt.figure()
    for H, ki in zip(Hs, ki_fluct):
        # plt.scatter(Ts, ki, label="H={}".format(H))
        x, y = rescale(Ts, ki, H)
        plt.scatter(x, y, label="H={}".format(H))

    plt.xscale("log")
    plt.xlabel("εˠ⁺ᵝ/H")
    plt.ylabel("χ/H^(1/̣δ-1)")
#    plt.xlim(1e-3, 1e1)
#    plt.ylim(0, 0.25)
    plt.legend()

    plt.figure()
    for H, ki in zip(H_mean, ki_diff):
        # plt.scatter(Ts, ki, label="H={}".format(H))
        x, y = rescale(Ts, ki, H)
        plt.scatter(x, y, label="H={}".format(H))

    plt.xscale("log")
    plt.xlabel("εˠ⁺ᵝ/H")
    plt.ylabel("χ/H^(1/̣δ-1)")
#    plt.xlim(1e-3, 1e1)
#    plt.ylim(0, 0.25)
    plt.legend()

    plt.show()

# Plot magnetisations vs Ts for various fields H
def MsvsTs():
    plt.figure()
    for i in range(len(Ms)): # len(Ms) = num of different H fields used
        plt.scatter(Ts, Ms[i], label="H={}".format(Hs[i]))

    plt.xlabel("T")
    plt.ylabel("M")
    plt.legend()

# To find Tc with only one value of H, H=0.
def MsvsTs0():
    plt.figure()
    for i in range(1): 
        plt.scatter(Ts, Ms[i], label="H={}".format(Hs[i]))

    plt.xlabel("T")
    plt.ylabel("M")
    plt.legend()



### MERGE DATASETS ###
# Mtot = [[], [], [], [], [], [], [], []] # 8 values of H
# sigmastot = [[], [], [], [], [], [], [], []]
Mtot = [[], [], [], [], [], []] # 6 values of H
sigmastot = [[], [], [], [], [], []]
Ttot = []

################

################
for i in range(2): ### INSERT NUMBER OF DATASETS HERE ###
    dataset = "dataset" + str(i+3) + ".hdf5" 
    f = h5py.File(dataset, "r")
    Hs, Ts, Ms, sigmas = retrieveData()
    Ttot.extend(Ts)
    for i in range(len(Ms)):
        Mtot[i].extend(Ms[i])
        sigmastot[i].extend(sigmas[i])
Ms = np.array(Mtot)
Ts = np.array(Ttot)
sigmas = np.array(sigmastot)

################

################

# data = np.array([[Hs],[Ms],[Ts],[sigmas]])
# np.savetxt("HsMsTsSigmas.txt", data, fmt="%s")

Hs = Hs.tolist()
Ts = Ts.tolist()
Ms = Ms.tolist()
sigmas = sigmas.tolist()
data = [Hs, Ts, Ms, sigmas]

with open('HsTsMsSigmasHighFields.txt', 'w') as f: # HsTsMsSigmasHighFields
    json.dump(data, f)

# # FOR HIGH FIELD DATASETS:
# Change the name of the txt file here
# Change the range in the dataset script and in the following line:    dataset = "dataset" + str(i+1) + ".hdf5" (Change (i+1))
# Change the numeber in lines:
#         # AVERAGES OVER SAMPLES
#         Ms.append(np.mean(sample[:, 8])) # 6 might need to be updated to the number of fields (?) # 6 USED TO BE 8
#         sigmas.append(np.std(sample[:, 8]))

#MsvsTs0()
#MsvsTs()
#scaling()
#plt.show()


# # PLOTING SCALING FUNCTION (requires plt.show(), actually at the bottom)
# Hs, Ts, Ms, sigmas = retrieveData()
# scaling()

#TEST#
# ki_fluct = 15 ** 3 * sigmas * sigmas / Ts
# ki_diff = (Ms[1::2] - Ms[::2]) / 0.1
# H_mean = (Hs[::2] + Hs[1::2]) / 2


# print(np.ndim(ki_diff))

# print(Ts)
# # plt.scatter(Ts, ki_diff)
# # plt.scatter(Ts, ki_fluct)
# # plt.scatter(Ts, Ms)

# plt.xlabel("T")
# plt.ylabel("Mz per site")
# plt.xlim(0.1, 1.2)
#END TEST#


# Used to be commented, from here...
# Ts, Ms, Es, sigmasE = retrieveData2()

# Ms = Ms.reshape((-1,))
# Es = Es.reshape((-1,))

# plt.figure()
# plt.scatter(Ts, Ms)
# plt.xlabel("T")
# plt.ylabel("Mz per site")
# plt.xlim(0.1, 1.2)
# # ...To here

def fit_M(Tmin = 0.5, Tmax = None):
    Tc = 0.67
    if not Tmax: Tmax = Tc
    mask = (Ts > Tmin) * (Ts < Tmax)
    
    epsilons = (Ts[mask] - Tc) / Tc
    Ms_mask = Ms[mask]

    x = np.log(np.abs(epsilons))
    y = np.log(Ms_mask)
    a, b = np.polyfit(x, y, 1)

    plt.figure()

    plt.scatter(Ts, Ms)
    plt.scatter(Ts[mask], np.exp(b) * (-epsilons) ** a, label="$\\beta$ = {}".format(a))
    
    plt.legend()
    plt.show()

def fit_M_nonlin(Tmin = 0.5, Tmax = 0.6):
    Tc = 0.67
    if not Tmax: Tmax = Tc
    mask = (Ts > Tmin) * (Ts < Tmax)
    
    epsilons = (Ts[mask] - Tc) / Tc
    Ms_mask = Ms[mask]

    fitfunc = lambda p, t: p[0] * (1 - t / p[1]) ** p[2]
    errfunc = lambda p, t, y: fitfunc(p, t) - y
    p0 = np.array([1, Tc, 0.37])
    p, out = optimize.leastsq(errfunc, p0[:], args=(Ts[mask], Ms_mask))
    print(p, out)
    
    plt.figure()

    plt.scatter(Ts, Ms, label="Data")
    plt.plot(Ts[mask], fitfunc(p, Ts[mask]), c="orange", label="Fit : $T_c = {:.2}, \\beta= {:.3}$".format(p[1], p[2]))

    plt.xlim(0.1, 1.2)
    plt.xlabel("T")
    plt.ylabel("Magnetization per site")
    plt.legend()
    plt.show()
    
def colormap():
    Tc = 0.67
    Tmaxs = np.linspace(0.1, Tc)
    deltas = np.linspace(0.05, 0.5)

    betas = np.zeros((len(Tmaxs), len(deltas)))
    for i, Tmax in enumerate(Tmaxs):
        for j, delta in enumerate(deltas):
            Tmin = Tmax - delta
            mask = (Ts > Tmin) * (Ts < Tmax)
            if mask.sum() >= 10:
                epsilons = (Ts[mask] - Tc) / Tc
                Ms_mask = Ms[mask]
                x = np.log(np.abs(epsilons))
                y = np.log(Ms_mask)
                a, b = np.polyfit(x, y, 1)

                fitfunc = lambda p, t: p[0] * (1 - t / p[1]) ** p[2]
                errfunc = lambda p, t, y: fitfunc(p, t) - y
                p0 = np.array([1, Tc, 0.37])
                p, out = optimize.leastsq(errfunc, p0[:], args=(Ts[mask], Ms_mask))
                
                betas[i, j] = p[2]

    plt.figure()

    im = plt.gca().imshow(betas)
    plt.xticks(5 * np.arange(len(Tmaxs[::5])), Tmaxs[::5])
    plt.yticks(5 * np.arange(len(deltas[::5])), (Tmaxs - deltas)[::5])

    plt.xlabel("Tmax")
    plt.ylabel("Tmin")
    
    plt.legend()
    plt.show()

    beta_expected = 0.37
    betas = np.abs(betas - beta_expected) / beta_expected

    plt.figure()

    im = plt.gca().imshow(betas, cmap="hot")
    plt.xticks(np.arange(len(Tmaxs)), Tmaxs)
    plt.yticks(np.arange(len(deltas)), deltas)

    plt.xlabel("Tmax")
    plt.ylabel("Tmin")
    
    plt.legend()
    plt.show()

# below here everything was commented excepted for the "plt.show()  #this was uncommented" line
# Tc = 0.67
# epsilons = (Ts[Ts < Tc] - Tc) / Tc
# plt.figure()
# x = np.log(np.abs(epsilons))
# y = np.log(Ms[Ts < Tc])
# a, b = np.polyfit(x[x > -1], y[x > -1], 1)
# print(a, b)
# plt.scatter(x, y)
# plt.plot(x, a * x + b, 'r')
# # plt.scatter(np.abs(epsilons), Ms[Ts < Tc])
# plt.xlabel("log(t)")
# plt.ylabel("log(M)")
# # plt.xscale("log")
# # plt.yscale("log")
# plt.grid(True)

# plt.figure()
# plt.scatter(Ts, Es)
# plt.xlabel("T")
# plt.ylabel("E per site")
# plt.xlim(0.1, 1.2)

# plt.figure()
# plt.scatter(Ts, 15 ** 3 * sigmasE * sigmasE / Ts / Ts, label="Fluctuation")
# plt.scatter(Ts, np.gradient(Es, Ts), label="Differentiation")
# plt.xlabel("T")
# plt.ylabel("dE/dT")
# plt.xlim(0.1, 1.2)

# Tmax = Ts[np.argmax(np.gradient(Es, Ts))]
# plt.axvline(x=Tmax, label="T = {:.2}".format(Tmax))

# plt.legend()

#plt.show()  #this was uncommented


# i=6
# ki_diff = (Ms[i+1] - Ms[i]) / 0.1
# ki_fluct = sigmas[i] * sigmas[i] / Ts

# plt.figure()

# plt.scatter(Ts, ki_diff, label="Diff")
# plt.scatter(Ts, 15 ** 3 * ki_fluct, label="Fluct")

# plt.xlabel("T")
# plt.ylabel("magnetic susceptibility")
# plt.legend()
# plt.title("With H = {}".format(Hs[i]))

# plt.show()

