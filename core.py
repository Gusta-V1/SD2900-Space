import json
import numpy as np
from dataclasses import dataclass


# Technical
FLOAT = np.float64
NUM_WORKERS = None

# Main options
OBJECT_DATASET = r"rb-750-1000-debris.json"

# Vehicle properties
DRV_DMASS = 200 # Dry mass (kg)
#DRV_GMASS = # Gross mass (kg)
#DRV_THRUST = 250*1e-3 # (N)
DRV_ISP = 10000 # (s)

# Debris mass dictionary (kg)
MASS_DICT = {
    "SL-8 R/B": 1435
}

# Earth Constants: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MU = 0.39860*1e6    # (km^3/s^2)
R = 6371.000        # (km)
J2 = 1082.63*1e-6
G0 = 9.80665*1e-3   # (km/s^2)


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
@dataclass(frozen=True)
class Object: # TODO: ADD WAY TO GET LOCATION/PARAMETERS AT TIMESTAMP
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

    def __repr__(self):
        return (
            "Object(\n"
            f"    name={self.name!r},\n"
            f"    object_id={self.object_id!r},\n"
            f"    norad_id={self.norad_id!r},\n"
            f"    epoch={self.epoch!r},\n"
            f"    rev_at_epoch={self.rev_at_epoch},\n"
            f"    a={self.a:.3f} km,\n"
            f"    n={self.n:.6e} rad/s,\n"
            f"    nd={self.nd:.6e} rad/s²,\n"
            f"    ndd={self.ndd:.6e} rad/s³,\n"
            f"    e={self.e:.8f},\n"
            f"    i={self.i:.6f} rad,\n"
            f"    raan={self.raan:.6f} rad,\n"
            f"    aper={self.aper:.6f} rad,\n"
            f"    M={self.M:.6f} rad,\n"
            f"    Bs={self.Bs:.6e} 1/R\n"
            ")"
        )

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
        np.deg2rad(FLOAT(entry["INCLINATION"])),
        np.deg2rad(FLOAT(entry["RA_OF_ASC_NODE"])),
        np.deg2rad(FLOAT(entry["ARG_OF_PERICENTER"])),
        np.deg2rad(FLOAT(entry["MEAN_ANOMALY"])),
        FLOAT(entry["BSTAR"])
    )
    return target

def load_objects(a_min: FLOAT=0) -> dict[str, Object]:
    objects=dict()
    with open(OBJECT_DATASET) as f:
        data = json.load(f)
        for entry in data:
            obj = read_tle(entry)

            # Apply target constraints
            if obj.a >= a_min:
                objects[obj.norad_id] = obj

    return objects


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