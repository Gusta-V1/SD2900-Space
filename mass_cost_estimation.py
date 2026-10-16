import json
import numpy as np
from core import FLOAT, MU, R, R_GAS, A_DISPOSAL, DRV_T, MASS_DICT, load_objects, burn_propellant, dv_a, dvs_rendezvous


# Main options
estimation_iterations = 100

# Properties
q_sp = 720e3        # Specific charge (Ws/kg)
p_sp = 84.5         # Specific power (W/kg)
#rq_sp = # Radiator specific radiation (W/kg) TODO: ONLY DO THERMALS IF TIME
tm_sp = 195         # Propellant tank specific volume-mass (MT Aerospace S-XTA-60 Xenon Tank) (kg/m^3)
p_Tank = 180*1e5    # Propellant tank pressure (Pa)
mp_mol = 0.131293   # Propellant (Xenon) molar mass (kg/mol)
rho_DRV = 79        # Vehicle density (kg/m^3)

# Initial estimations
P_propulsion = 5000 # Propulsion power demand (W)
P_misc = 60         # Misc power demand (W)
p_margin = 0.1      # Percent power budget margin
m_Thruster = 0 # Thruster mass (kg)
eps_structure = 0.15# Structural factor
m_margin = 0.1      # Mass margin factor
eta_Battery = 0.8   # Battery depth of discharge


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
#t_0 = # Time of insertion(at periapsis) (timestamp(s???)?)


# Estimation functions
def propellant_volume(mp: FLOAT): # NOTE: Not very accurate at high pressures.
    """Calculate propellant volume (m^3) for a fuel mass (kg)"""
    return mp * R_GAS*DRV_T / (p_Tank * mp_mol)

def tank_mass(mp): # NOTE: A bit too simplistic.
    """Estimate fuel tank mass (kg) from fuel mass (kg)"""
    return tm_sp * propellant_volume(mp)

def req_battery_capacity(P: FLOAT, a: FLOAT):
    """
    Calculate required battery capacity (Ws), assumes quasi-circular orbit.
    
    Params:
        P (FLOAT): - Power consumption
        a (FLOAT): - Semi-major axis
    """
    te = 2*np.sqrt(a**3 / MU) * np.asin(R/a)
    q = P * te / eta_Battery

    return q

def req_power_generation(P: FLOAT, a: FLOAT):
    """Calculate required power generation (W), assumes quasi-circular orbit.
    
    Params:
        P (FLOAT): - Power consumption
        a (FLOAT): - Semi-major axis
    """
    te = 2*np.sqrt(a**3 / MU) * np.asin(R/a)
    tp = 2*np.pi*np.sqrt(a**3 / MU)
    p_gen = P * (1 + te / (tp - te))

    return p_gen

def radiator_mass(P: FLOAT, T: FLOAT): #TODO: IF TIME, RADIATOR MASS
    pass


# Power requirements
P_req = (P_propulsion + P_misc) * (1 + p_margin)


# Estimate DRV mass and cost
dv_tot_dict={}
mp_tot_dict={}
ms_dict={}
m_tot_dict={}
for obj in cluster:
    # Power requirements
    P_Solar = req_power_generation(P_req, obj.a)
    Q_Battery = req_battery_capacity(P_req, obj.a)

    # Power and thermal systems mass
    #m_Radiator = radiator_mass(P_req) #TODO: CALCULATE THERMALS IF TIME
    m_Solar = P_Solar / p_sp
    m_Battery = Q_Battery / q_sp

    print(f"Power generation: {P_Solar} W")
    print(f"Battery capacity: {Q_Battery/3600} Wh")
    print(f"Solar panel mass: {m_Solar} kg")
    print(f"Battery mass: {m_Battery} kg")

    # Initial structural mass estimation
    m_Tank=0 # Initial guess
    ms = m_Tank + m_Thruster + m_Solar + m_Battery
    ms += ms*eps_structure    # Structural mass
    m_Margin = ms * m_margin  # Margin mass

    # Iteratively update fuel requirement estimate
    for j in range(estimation_iterations):
        dv_rendezvous = sum(dvs_rendezvous(a0, e0, i0, raan0, aper0, obj))
        dv_disposal = dv_a(obj.a, A_DISPOSAL)


        # Calculate required propellant masses
        # mp1: Rendevous
        # mp2: Disposal
        mp2 = burn_propellant(ms + m_Margin + MASS_DICT[obj.name], dv_disposal)
        mp1 = burn_propellant(ms + m_Margin + mp2, dv_rendezvous)

        mp_tot = mp1 + mp2

        # Mass estimation
        m_Tank = tank_mass(mp_tot)
        ms = m_Tank + m_Thruster + m_Solar + m_Battery
        ms += ms*eps_structure                # Structural mass
        m_Margin = (ms + mp_tot) * m_margin   # Margin mass
    ms += m_Margin


    dv_tot = dv_rendezvous + dv_disposal
    print(f"{obj.norad_id}: Rendezvous dv_1={dv_rendezvous*1000} m/s")
    print(f"{obj.norad_id}: Rendezvous m_p1={mp1} kg")
    print(f"{obj.norad_id}: Disposal dv_2={dv_disposal*1000} m/s")
    print(f"{obj.norad_id}: Disposal m_p2={mp2} kg")
    print(f"{obj.norad_id}: Total dv_tot={dv_tot*1000} m/s")
    print(f"{obj.norad_id}: Total m_p={mp_tot} kg")
    print(80*"-")

    dv_tot_dict[obj.norad_id] = dv_tot
    mp_tot_dict[obj.norad_id] = mp_tot
    ms_dict[obj.norad_id] = ms
    m_tot_dict[obj.norad_id] = ms + mp_tot

# Print requirements
print("Design values for propellant and vehicle mass")
print(f"Total propellant: {max(mp_tot_dict.values())} kg")
print(f"Total DRV dry mass (with margin): {max(ms_dict.values())} kg")
print(f"Total DRV mass at launch (with margin): {max(m_tot_dict.values())} kg")

print("Design values for DRV power")
print(f"Power requirement (with margin): {P_req} W")

print("Design values for DRV volume")
print(f"Estimated volume: {max(m_tot_dict.values())/rho_DRV} m^3")