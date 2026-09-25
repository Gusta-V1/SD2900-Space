import json
import random
import numpy as np
from dataclasses import dataclass
from scipy.optimize import minimize
from scipy.special import comb


# Precision
FLOAT = np.float64

# Earth Constants: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MU = 0.39860*1e6  # km^3/s^2
R = 6371.000      # km
J2 = 1082.63*1e-6

# Vehicle properties
# mass, thrust...

# Utilities
# {"OBJECT_NAME":"COSMOS 2251",
# "OBJECT_ID":"1993-036A",
# "EPOCH":"2026-09-10T06:23:36.172896", # Timestamp of data
# "MEAN_MOTION":14.33268558,            # Mean motion (rot/day)
# "ECCENTRICITY":0.0025219,             # Eccentricity
# "INCLINATION":74.0403,                # Inclination (deg)
# "RA_OF_ASC_NODE":179.7474,            # RAAN (deg)
# "ARG_OF_PERICENTER":200.8184,         # Argument of periapsis (deg)
# "MEAN_ANOMALY":159.1942,              # Mean anomaly (deg)
# "EPHEMERIS_TYPE":0,                   # Ephemeris type (0 = SGP4/SDP4)
# "CLASSIFICATION_TYPE":"U",            # Classification (unclassified, classified, secret)
# "NORAD_CAT_ID":22675,                 # NORAD Catalog Number
# "ELEMENT_SET_NO":999,                 # Element set number
# "REV_AT_EPOCH":73663,                 # Number of revolutions completed at epoch
# "BSTAR":2.2572e-5,                    # Inverse ballistic coefficient (1/R)
# "MEAN_MOTION_DOT":3.6e-7,             # Mean motion first time derivative divided by two (rot/day²)
# "MEAN_MOTION_DDOT":0}                 # Mean motion second time derivative divided by six (rot/day³)
@dataclass
class Object:
    """
    Satellite Object class

    Params:
        name (str): - Plaintext name
        object_id (str): - International identifier
        norad_id (str): - NORAD Catalog Number
        epoch (str): - Timestamp of observation
        rev_at_epoch (int) - Number of revolutions completed at epoch
        a (FLOAT): - Semi-major axis (km)
        n (FLOAT): - Mean motion (rad/s)
        nd (FLOAT): - Mean motion first time derivative (rad/s²)
        ndd (FLOAT): - Mean motion second time derivative (rad/s³)
        e (FLOAT): - Eccentricity
        i (FLOAT): - Inclination (rad)
        raan (FLOAT): - RAAN (rad)
        aper (FLOAT): - Argument of periapsis (rad)
        M (FLOAT): - Mean anomaly (rad)
        Bs (FLOAT): - Inverse ballistic coeffient BSTAR (1/R), assumes rho=0.15696615 kg/(m^2*R) (https://celestrak.org/NORAD/documentation/spacetrk.pdf)
    """
    name: str
    object_id: str
    norad_id: str
    epoch: str
    rev_at_epoch: int
    a: FLOAT
    n: FLOAT
    nd: FLOAT
    ndd: FLOAT
    e: FLOAT
    i: FLOAT
    raan: FLOAT
    aper: FLOAT
    M: FLOAT
    Bs: FLOAT

def mean_motion_to_a(n: FLOAT) -> tuple[FLOAT, FLOAT]:
    """
    Calculate semi-major axis (km) from mean motion (rad/s).
    """
    return np.cbrt(MU/n**2)

def TLE_mean_motion_to_si(n_ra: FLOAT, nd_ra: FLOAT, ndd_ra: FLOAT) -> tuple[FLOAT, FLOAT, FLOAT]:
    """
    Converts TLE mean motion and its time derivatives from rotations per day radians per second.
    """
    n = n_ra * 2*np.pi / 86400
    nd = 2* nd_ra * 2*np.pi / 86400**2 # Multiplied by 2 to counteract TLE scaling
    ndd = 6 * ndd_ra * 2*np.pi / 86400**3 # Multiplied by 6 to counteract TLE scaling
    return (n, nd, ndd)

def read_tle(entry):
    """
    Read TLE data for a JSON entry, and return a satellite Object object.
    """
    n, nd, ndd = TLE_mean_motion_to_si(FLOAT(entry["MEAN_MOTION"]), FLOAT(entry["MEAN_MOTION_DOT"]), FLOAT(entry["MEAN_MOTION_DDOT"]))
    a = mean_motion_to_a(n)

    target = Object(
        str(entry["OBJECT_NAME"]),
        str(entry["OBJECT_ID"]),
        str(entry["NORAD_CAT_ID"]),
        str(entry["EPOCH"]),
        int(entry["REV_AT_EPOCH"]),
        a,
        n,
        nd,
        ndd,
        FLOAT(entry["ECCENTRICITY"]),
        np.deg2rad(entry["INCLINATION"], dtype=FLOAT),
        np.deg2rad(entry["RA_OF_ASC_NODE"], dtype=FLOAT),
        np.deg2rad(entry["ARG_OF_PERICENTER"], dtype=FLOAT),
        np.deg2rad(entry["MEAN_ANOMALY"], dtype=FLOAT),
        FLOAT(entry["BSTAR"])
    )
    return target

def filter_candidate(target: Object, base: Object, a_lim: float, i_lim: float, raan_lim: float):
    Da = np.abs(target.a - base.a)
    Di = np.abs(target.i - base.i)
    Draan = np.abs(target.raan - base.raan)
    return Da < a_lim and Di < i_lim and Draan < raan_lim


# Delta v (not including drifts at the moment, quasi-circular)
def dv_a(a0: FLOAT, at: FLOAT) -> FLOAT:
    return np.abs(np.sqrt(MU/a0) - np.sqrt(MU/at))

def dv_e(e0: FLOAT, et: FLOAT, a: FLOAT) -> FLOAT:
    return np.sqrt(MU/a) * 2/3 * np.abs(np.asin(e0) - np.asin(et))

def dv_i(i0: FLOAT, it: FLOAT, a: FLOAT) -> FLOAT:
    v = np.sqrt(MU/a) # Assume quasi-circular
    return v * np.sqrt(2 - 2*np.cos((np.pi*(it - i0))/2))

def dv_raan(raan0: FLOAT, raant: FLOAT, a: FLOAT, i: FLOAT) -> FLOAT:
    return np.pi/2 * np.sqrt(MU/a) * np.abs(raant - raan0) * np.sin(i)


# Oblateness effects
def d_raan(a: FLOAT, e: FLOAT, i: FLOAT) -> FLOAT: # RAAN
    n = np.sqrt(MU/a**3)
    return -(3*n*J2*R**2) / (2*a**2*(1 - e**2))**2 * np.cos(i)

def d_aper(draan: FLOAT, i: FLOAT) -> FLOAT: # Argument of the periapsis
    return draan / np.cos(i) * (5/2 * np.sin(i)**2 - 2)

# Calculate maximum delta-v required for cluster
def cluster_dv_max(x, cluster: tuple[Object]):
    """
    Function to calculate maximum delta-v required for a cluster from given initial parameters.

    Params:
        x (np.ndarray([a, e, i, raan])): - Initial orbit parameters
        cluster (tuple[Object]): - List of targets
    """
    a0, e0, i0, raan0 = x # NOTE: These get modified with some maneuvers

    dvtot=np.empty(len(cluster))
    for j in range(len(cluster)):
        obj = cluster[j]
        a = obj.a
        e = obj.e
        i = obj.i
        raan = obj.raan

        # Burn plan logic
        # If final orbit apoapsis lower than current orbit: 
        # If final orbit inclination lower than current orbit:
        #   1. Change inclination to match target
        #   2. Change RAAN to match target (add lead, drift differs until eccentricity and semi major is identical)
        # If final orbit inclination higher than current orbit:
        #   1. Change RAAN to match target (add lead, drift differs until inclination, eccentricity and semi major is identical))
        #   2. Change inclination to match target
        # 3. Change eccentricity to match target (set argument of periapsis, include lead, drift differs until semi-major is identical)
        # 4. Change altitude to match semi-major axis and phase (position in orbit relative target).
        # 
        # If final orbit apoapsis higher than current orbit:
        # 1. Change altitude to match semi-major axis and phase (position in orbit relative target).
        # 2. Change eccentricity to match target (set argument of periapsis, include lead, drift differs until inclination and raan is identical)
        # If final orbit inclination lower than current orbit:
        #   3. Change inclination to match target
        #   4. Change RAAN to match target
        # If final orbit inclination higher than current orbit:
        #   3. Change RAAN to match target (add lead, drift differs until inclination is identical)
        #   4. Change inclination to match target

        if a<a0:
            if i<i0:
                dvi = dv_i(i0, i, a0)
                dvraan = dv_raan(raan0, raan, a0, i)
            else:
                dvraan = dv_raan(raan0, raan, a0, i0)
                dvi = dv_i(i0, i, a0)
            dve = dv_e(e0, e, a0)
            dva = dv_a(a0, a)
        else:
            dva = dv_a(a0, a)
            dve = dv_e(e0, e, a)
            if i<i0:
                dvi = dv_i(i0, i, a)
                dvraan = dv_raan(raan0, raan, a, i)
            else:
                dvraan = dv_raan(raan0, raan, a, i0)
                dvi = dv_i(i0, i, a)

        dvtot[j] = dva + dve + dvi + dvraan

    return np.max(dvtot)


# Process data
targets=dict()
with open('cosmos-2251-debris.json') as f:
    data = json.load(f)

    for entry in data:
        target = read_tle(entry)
        targets[target.norad_id] = target

num_objects = len(targets)
print(f"Total objects: {num_objects}")


# Main program
#   Note: That leads require thrust and weight data
num_targets = 8 # Number of targets
# Do random search for combinations of n sats with i iterations (cap to nCr(num_objects, n))
search_iters = 10000
search_iters = min(search_iters, int(comb(num_objects, num_targets)))

clusters=[]
for i in range(search_iters):
    base_key = random.choice(list(targets.keys()))
    base = targets[base_key]
    candidates_keys = list(targets.keys() - set(base_key))

    # Filter remaining targets based on inclination and raan (estimated dv costs at a=600 km i=74 deg)
    inclination_limit = np.deg2rad(5.0) # (deg) (about 3.4 km/s)
    raan_limit = np.deg2rad(5.0)        # (deg) (about 3.53 km/s, orbital precession only 0.45 deg per month)
    a_limit = 200                       # (km) (about 3.45 km/s when raising)

    candidates = [targets[key] for key in candidates_keys if filter_candidate(targets[key], base, a_limit, inclination_limit, raan_limit)]
    if len(candidates) < num_targets:
        continue

    # Sample candidates and minimize
    cluster = list([base] + random.choices(candidates, k=num_targets-1))

    #Initial params
    a0 = np.mean([t.a for t in candidates])
    e0 = np.mean([t.e for t in candidates])
    i0 = np.mean([t.i for t in candidates])
    raan0 = np.mean([t.raan for t in candidates])
    x0 = np.array([a0, e0, i0, raan0], dtype=FLOAT)

    bounds=[
        (R, None),      # Semi-major axis
        (0, 1),         # Eccentricity
        (0, np.pi/2),   # Inclination
        (0, 2*np.pi)    # RAAN
    ]
    result = minimize(cluster_dv_max, x0, args=(candidates,), method='SLSQP', bounds=bounds)

    clusters.append((cluster, result))

# Sort results
# (cluster, OptimizeResult) https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
clusters = sorted(clusters, key=lambda c: c[1].fun)
    
print(f"Optimal cluster: {clusters[0][0]}")
print(f"Optimal Parameters: a0 = {clusters[0][1].x[0]} km, e0 = {clusters[0][1].x[1]}, i0 = {clusters[0][1].x[2]} deg, RAAN0 = {clusters[0][1].x[3]}")
print(f"Minimum Total Delta V (km/s): {clusters[0][1].fun} km/s")

print(f"50:th Optimal cluster: {clusters[50][0]}")
print(f"50:th Optimal Parameters: a0 = {clusters[50][1].x[0]} km, e0 = {clusters[50][1].x[1]}, i0 = {clusters[50][1].x[2]} deg, RAAN0 = {clusters[50][1].x[3]}")
print(f"50:th Minimum Total Delta V (km/s): {clusters[50][1].fun} km/s")