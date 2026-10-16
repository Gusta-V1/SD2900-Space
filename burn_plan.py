import json
import numpy as np
from core import FLOAT, NUM_WORKERS, R, G0, A_DISPOSAL, DRV_DMASS, MASS_DICT, DRV_THRUST, Object, load_objects, burn_propellant, dv_a, dvs_rendezvous



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
aper0 = orbit["aper"]   # (rad)
#t_0 = # Time of insertion(at periapsis) (timestamp?)


# Calculate burn plan (with times, dv, (and fuel requirements)) for each DRV
# When need debris mass, just check with MASS_DICT[name]

# Calculate burn times and mission duration for the longest mission
times={}
for obj in cluster:
    dv_rendezvous = sum(dvs_rendezvous(a0, e0, i0, raan0, aper0, obj))
    dv_disposal = dv_a(obj.a, A_DISPOSAL)

    mp2 = burn_propellant(DRV_DMASS + MASS_DICT[obj.name], dv_disposal)
    mp1 = burn_propellant(DRV_DMASS + mp2, dv_rendezvous)

    # Approximate masses
    m1 = DRV_DMASS + mp2 + mp1/2
    m2 = DRV_DMASS + mp2/2 + MASS_DICT[obj.name]

    # Times
    t1 = m1 * dv_rendezvous / DRV_THRUST
    t2 = m2 * dv_disposal / DRV_THRUST
    t_tot = t1 + t2

    times[obj.norad_id] = (t_tot, t1, t2)

# Get longest mission and print info
norad_id, (t_tot, t1, t2) = max(times.items(), key=lambda item: item[1][0])

print(f"Longest mission: NORAD {norad_id}")
print(f"Total duration: {t_tot / 86400} days")
print(f"Rendezvous burn: {t1 / 86400} days")
print(f"Disposal burn: {t2 / 86400} days")
