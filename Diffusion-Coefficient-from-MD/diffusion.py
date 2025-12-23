import matplotlib.pyplot as plt 
import os
import numpy as np
import scipy as sp


# Specify number of frames and number of replicate simulations
num_frames = 1000
reps = 10

# get MSD across all 20 sims into a numpy array
files = os.listdir()
agr_files = []

for file in files:
    if '.agr' in file:
        agr_files.append(file)
    
# Initialize an empty numpy array with the proper dimensions
MSD_array = np.empty(shape = (num_frames, reps))


for file, column in zip(agr_files, range(MSD_array.shape[1])):

    time = []
    MSD = []

    # Extract data from each file and append to an empty list
    with open(file) as f:
        data = f.readlines()

    line_ind = 0
    while "a-D-Glc[X]" not in data[line_ind]:
        
        if "@" in data[line_ind]:
            line_ind += 1
        else:
            # convert from angstroms^2 
            MSD.append(float(data[line_ind].split()[1]))
            time.append(float(data[line_ind].split()[0]))
            line_ind += 1

    # Iterate through the current MSD list and append each value to the corresponding row and index in
    # the numpy array
    for value, row in zip(MSD,range(MSD_array.shape[0])):
        MSD_array[row,column] = value


# Calculate average MSD across all 20 sims and perform linear regression.
avg_MSD = np.mean(MSD_array, axis = 1)

# finally, convert time list to a numpy array for later
time = np.asarray(time)

# function to calculate diffusion coefficient
def DC(time, avg_MSD):
    """
    Takes in two numpy arrays for the time and avg MSD across replicates and computes the diffusion coefficient from them using linear regression.
    """
    reg = sp.stats.linregress(time, avg_MSD)
    D = ((reg.slope*(10**12))/(10**20))/6
    return D


# Now make the D v time plot
fig, ax = plt.subplots(1)

DCs = []

DC_times = range(50,10001,50)
for i in range(5,1001,5):
    DCs.append(DC(time[:i], avg_MSD[:i]))
  
# plot DC v time
ax.plot(DC_times, DCs)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('D (m$^{2}$/s)')
ax.set_ylim(0,15*10**-10)
ax.set_title("a-D-Glc, OPC, GLYCAM06 Charges")
fig.savefig("D_v_time.png")

print(f"\nFinal D: {DCs[-1]}\n")



