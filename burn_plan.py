import json
import numpy as np
import target_selection as tg

with open("cluster-solution.json", encoding="utf-8") as f:
    solution = json.load(f)

# Load cluster and insertion orbit
selected_objects = solution["cluster"]
cluster = [tg.targets[i] for i in selected_objects]
orbit = solution["insertion_orbit"]
a0 = orbit["a"]         # (km)
e0 = orbit["e"]
i0 = orbit["i"]         # (rad)
raan0 = orbit["raan"]   # (rad)
#t_0 = # Time of insertion(at periapsis) (timestamp?)

