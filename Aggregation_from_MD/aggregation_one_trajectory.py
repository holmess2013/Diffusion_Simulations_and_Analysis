import aggregation
import statistics 
from collections import Counter
import numpy as np

# The input for this program should be an autoimaged trajectory of multiple sugars with water stripped, in PDB format.

sys_name = "B-D-Glc_OPC"

# Read in trajectory
with open(f"{sys_name}_wrapped.pdb") as f:
    traj = f.read().splitlines()

# Make a list of all frames, where each element is one frame in pdb format (also a list).
frames = []
current_frame = []

# This will separate each frame using ENDMDL as a delimiter.
for line in traj:
    if "ENDMDL" in line:
        frames.append(current_frame)
        current_frame=[]
    else:
        current_frame.append(line)


output_file = open(f"avg_agg_size_{sys_name}.txt","w")

avg_agg_size_per_frame = []
max_size_per_frame = []

for frame, i in zip(frames, range(len(frames))):
    print(f"\nFrame {i}")
    avg_agg_size, max_size = aggregation.calculate_aggregates(frame)
    avg_agg_size_per_frame.append(avg_agg_size)
    max_size_per_frame.append(max_size)
    output_file.write(f"{avg_agg_size}\n")

output_file.close()

print('\nAverage aggregate size over full simulation: ', statistics.mean(avg_agg_size_per_frame))
print('Max size over full simulation: ', max(max_size_per_frame))

print(Counter(max_size_per_frame))

