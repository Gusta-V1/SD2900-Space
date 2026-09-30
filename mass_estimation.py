import json
import numpy as np
from core import G0, A_DISPOSAL, DRV_ISP, DRV_DMASS, MASS_DICT, load_objects, dv_a, dvs_rendezvous



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
#t_0 = # Time of insertion(at periapsis) (timestamp(s???)?)


# Calculate dv_rendevouz for each DRV (maybe restructure everything and move a lot to burn_plan)
# This is currently just for data entry, but later good for final burn plan
dv_tot_max=0
mp_tot_max=0
dv_rendezvous_max=0
mp1_max=0
dv_disposal_max=0
mp2_max=0
for obj in cluster:
    a = obj.a
    e = obj.e
    i = obj.i
    raan = obj.raan
    
    dv_rendezvous = sum(dvs_rendezvous(a0, a, e0, e, i0, i, raan0, raan))
    dv_disposal = dv_a(obj.a, A_DISPOSAL)

    # Calculate required propellant masses
    # mp1: Rendevous
    # mp2: Disposal
    mp2 = (DRV_DMASS + MASS_DICT[obj.name])*(np.exp(dv_disposal/(DRV_ISP*G0)) - 1)
    mp1 = (DRV_DMASS + mp2)*(np.exp(dv_rendezvous/(DRV_ISP*G0)) - 1)

    print(f"{obj.norad_id}: Rendezvous dv_1={dv_rendezvous*1000} m/s")
    #print(f"{obj.norad_id}: Rendezvous m_p1={mp1} kg")

    print(f"{obj.norad_id}: Disposal dv_2={dv_disposal*1000} m/s")
    #print(f"{obj.norad_id}: Disposal m_p2={mp2} kg")


    # Clunky way of getting maxes - TODO: IF THIS IS EVER USED, FIX THIS SHIT
    dv_tot = dv_rendezvous + dv_disposal
    mp_tot = mp1 + mp2
    print(f"{obj.norad_id}: Total dv_tot={dv_tot*1000} m/s")
    print(f"{obj.norad_id}: Total m_p={mp_tot} kg")
    print(80*"-")
    if dv_tot > dv_tot_max:
        dv_tot_max = dv_tot
    if mp_tot > mp_tot_max:
        mp_tot_max = mp_tot

    if dv_rendezvous > dv_rendezvous_max:
        dv_rendezvous_max = dv_rendezvous
    if mp1 > mp1_max:
        mp1_max = mp1
    if dv_disposal > dv_disposal_max:
        dv_disposal_max = dv_disposal
    if mp2 > mp2_max:
        mp2_max = mp2

# Print requirements
print("Design values for propellant and vehicle mass")
print(f"Total propellant: {mp_tot_max} kg")
print(f"Total DRV mass at launch: {mp_tot_max+DRV_DMASS} kg")