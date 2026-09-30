import json
import numpy as np
from core import FLOAT, NUM_WORKERS, R, G0, A_DISPOSAL, DRV_DMASS, MASS_DICT, DRV_THRUST, Object, load_objects



# Load data
targets = load_objects()
with open("cluster-solution.json", encoding="utf-8") as f:
    solution = json.load(f)

# Load cluster and insertion orbit
selected_objects = solution["cluster"]
cluster = [targets[i] for i in selected_objects]
orbit = solution["insertion_orbit"]
a0 = orbit["a"]         # (km)
e0 = orbit["e"]
i0 = orbit["i"]         # (rad)
raan0 = orbit["raan"]   # (rad)
#t_0 = # Time of insertion(at periapsis) (timestamp?)


# Calculate burn plan (with times, dv, (and fuel requirements)) for each DRV
# When need debris mass, just check with MASS_DICT[name]

