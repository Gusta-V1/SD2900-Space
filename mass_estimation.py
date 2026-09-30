import json
import numpy as np
from core import FLOAT, NUM_WORKERS, R, G0, DRV_ISP, DRV_DMASS, MASS_DICT, Object, load_objects, dv_a


# Assumed

# Deorbit
a_deorbit = R + 500 # ish (km) (check in matlab after mass estimation)


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

# Get max required deorbit Delta v
dv_deorbit = max([dv_a(obj.a, a_deorbit) for obj in cluster])
#GET MAX MASS

dv_rendevous = 0.7 # Add margin? (km/s)

# Calculate propellant masses
mp2 = (DRV_DMASS + 1435)*(np.exp(dv_deorbit/(DRV_ISP*G0)) - 1) #TODO, CALCULATE max per probe to support diff target vars

mp1 = (DRV_DMASS + mp2)*(np.exp(dv_deorbit/(DRV_ISP*G0)) - 1)

print(mp2)
print(mp1)
print(mp1+mp2+DRV_DMASS)

print(dv_deorbit)