"""A file to run simulations and generate the hdf5 dataset"""

import numpy as np
from os import system
import h5py

# used to generate dataset (?) (you need to manually change the name of the dataset below)
Hs = np.array([0.1, 0.2,
              1, 1.1,
              2, 2.1,
              3, 3.1])

# used to generate dataset0 (?). Zero field only needed for phase transition and not for scaling function
#Hs = np.array([0])

Tc = 0.49
Ts = np.linspace(0, 5, 40) # 20 # Ts = np.linspace(0.1, 2, 200) 
Ts = Ts*np.log(10)/1.33 # 1.33 should be gamma + beta
Ts = Tc + Tc*np.exp(Ts)

config = "big.in"
outdir = "" #"data/big/"
Nanneal = 5000 # to be adjusted with the total number of thermalization steps

outfile = outdir + "out"
statefile = outdir + "state"

Nsample = 13
samples = ["sample_{}".format(i+1) for i in range(Nsample)]

# CHANGE dataset.hdf5 WITH dataset0.hdf5 WHEN CHANGING THE VALUES FOR H
with h5py.File("dataset3.hdf5", "w") as f:
    f.attrs["config"] = config
    f.attrs["Nanneal"] = Nanneal

    # I would need to add a for here looping the rest of the code
    # to loop in the various ds values. I might want to launch sim
    # for differents values separately though.
    
    # grp.attrs["ds"] = ds would need to have ds in here like Hs

    for sample in samples:
        print("\n*** NEW SAMPLE # {} ***".format(sample))
        for H in Hs:
            print("*** # {} ***".format(sample)) # ADDED
            print("*** NEW H : Starting H={} ...".format(H)) # ADDED
            # print("Starting H={} ...".format(H)) # Removed
            grp = f.require_group("H={}".format(H))
            grp.attrs["H"] = H
            
            # First step
            grp2 = grp.require_group("T={}".format(Ts[0]))
            grp2.attrs["T"] = Ts[0]
            # Run the simulation
            system("./sim {} \"filename={}\" \"temperature={}\" \"H=(0 0 {})\" \"outstate={}\""
                   .format(config, outfile, Ts[0], H, statefile)) # could add \"ds={}\" and ds in format() to impose ds from gen like Ts
            # Store data
            data = np.loadtxt(outfile)
            dset = grp2.create_dataset(sample, data=data)
            dset.attrs["E"] = np.mean(data[:, 5])
            dset.attrs["varE"] = np.var(data[:, 5])
            dset.attrs["M"] = np.mean(data[:, 8])
            dset.attrs["varM"] = np.var(data[:, 8])
            # Clear temporary files
            system("rm {}".format(outfile))
            
            # Then annealing
            for T in Ts[1:]:
                print("*** # {} ***".format(sample)) # ADDED
                print("Starting H={} ...".format(H)) # ADDED
                print("Starting T={} ...".format(T))
                grp2 = grp.require_group("T={}".format(T))
                grp2.attrs["T"] = T
                
                # Run the simulation
                system("./sim {} \"filename={}\" \"temperature={}\" \"H=(0 0 {})\" \"outstate={}\" \"instate={}\" \"Nthermal={}\" "
                       .format(config, outfile, T, H, statefile, statefile, Nanneal))
                # Store data
                data = np.loadtxt(outfile)
                dset = grp2.create_dataset(sample, data=data)
                dset.attrs["E"] = np.mean(data[:, 5])
                dset.attrs["varE"] = np.var(data[:, 5])
                dset.attrs["M"] = np.mean(data[:, 8])
                dset.attrs["varM"] = np.var(data[:, 8])
                # Clear temporary files
                system("rm {}".format(outfile))

