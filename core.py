import json
import numpy as np
from dataclasses import dataclass
from datetime import datetime


# Technical
FLOAT = np.float64
NUM_WORKERS = None

# Earth Constants: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MU = 0.39860*1e6    # (km^3/s^2)
R = 6371.000        # (km)
J2 = 1082.63*1e-6
G0 = 9.80665*1e-3   # (km/s^2)

# Universal Constants
T_CMBR = 2.72548    # CMBR Temperature (K): https://lambda.gsfc.nasa.gov/resources/lambda_graphics/cmb_monopole.html
R_GAS = 8.314462618 # Universal gas constant (J/(K*mol)): https://physics.nist.gov/cgi-bin/cuu/Value?r


# Main options
OBJECT_DATASET = r"rb-750-1000-debris.json"

# Vehicle properties
DRV_DMASS = 150             # Dry mass (kg)
DRV_THRUST = 250*1e-3*1e-3  # Thrust (kN)
DRV_ISP = 10000             # Isp (s)
DRV_T = 273.15 + 20         # Vehicle temperature (K)

# Pre-calculated values
A_DISPOSAL = R + 500        # Disposal orbit semi-major axis (km): calculated in MATLAB

# Debris mass dictionary (kg)
MASS_DICT = {
    "SL-8 R/B": 1435
}



# Object class
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

    def at_epoch(self, epoch: datetime) -> tuple[FLOAT, FLOAT, FLOAT]:
        """
        Calculates orbital paramers (raan, aper, M) at a given epoch, based on the last observation.
        """
        observation_epoch = datetime.fromisoformat(self.epoch)
        elapsed = (epoch - observation_epoch).total_seconds()

        TAU = 2*np.pi

        # Calculate parameters
        #   NOTE: n, nd are not updated as they are tuned to fit SGP4
        draan = d_raan(self.a, self.e, self.i)
        daper = d_aper(draan, self.i)

        raan = (self.raan + draan*elapsed) % TAU
        aper = (self.aper + daper*elapsed) % TAU
        m = (self.M + self.n*elapsed) % TAU

        return (raan, aper, m)


# Utilities
def mean_motion_to_a(n: FLOAT) -> tuple[FLOAT, FLOAT]:
    """
    Calculate semi-major axis (km) from mean motion (rad/s).
    """
    return np.cbrt(MU/n**2)

def burn_time(dv: FLOAT, m: FLOAT):
    """
    Calculates the burn time (s) for a given change in velocity (km/s) for a DRV.

    Params:
        dv (FLOAT): - Desired change in velocity (km/s)
        m (FLOAT): - Mean propelled mass (kg)
    """
    return dv * m / DRV_THRUST

def TLE_mean_motion_to_si(n_ra: FLOAT, nd_ra: FLOAT, ndd_ra: FLOAT) -> tuple[FLOAT, FLOAT, FLOAT]:
    """
    Converts TLE mean motion and its time derivatives from rotations per day radians per second.
    """
    n = n_ra * 2*np.pi / 86400
    nd = 2* nd_ra * 2*np.pi / 86400**2      # Multiplied by 2 to counteract TLE scaling
    ndd = 6 * ndd_ra * 2*np.pi / 86400**3   # Multiplied by 6 to counteract TLE scaling
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

def load_objects(a_min: FLOAT=R) -> dict[str, Object]:
    objects=dict()
    with open(OBJECT_DATASET) as f:
        data = json.load(f)
        for entry in data:
            obj = read_tle(entry)

            # Apply target constraints
            if obj.a >= a_min:
                objects[obj.norad_id] = obj

    return objects

def wrap_angle_delta(delta: FLOAT):
    """Calculate the wrapped angle difference (-pi, pi) for a given delta."""
    return np.arctan2(np.sin(delta), np.cos(delta))


# Delta v for parameter changes (quasi-circular)
def dv_a(a0: FLOAT, at: FLOAT) -> FLOAT:
    return np.abs(np.sqrt(MU/a0) - np.sqrt(MU/at))

def dv_e(e0: FLOAT, et: FLOAT, a: FLOAT) -> FLOAT:
    return np.sqrt(MU/a) * 2/3 * np.abs(np.asin(e0) - np.asin(et))

def dv_i(i0: FLOAT, it: FLOAT, a: FLOAT) -> FLOAT:
    v = np.sqrt(MU/a) # Assume quasi-circular
    return v * np.sqrt(2 - 2*np.cos((np.pi*(it - i0))/2))

def dv_raan(raan0: FLOAT, raant: FLOAT, a: FLOAT, i: FLOAT) -> FLOAT:
    return np.pi/2 * np.sqrt(MU/a) * np.abs(wrap_angle_delta(raant - raan0)) * np.sin(i)

def dv_aper(aper0: FLOAT, apert: FLOAT, a: FLOAT, e: FLOAT) -> FLOAT:
    return np.sqrt(MU/a) * 2/3 * e/np.sqrt(1 - e) * np.abs(wrap_angle_delta(apert - aper0))


# Oblateness effects
def d_raan(a: FLOAT, e: FLOAT, i: FLOAT) -> FLOAT: # RAAN
    n = np.sqrt(MU/a**3)
    return -(3*n*J2*R**2) / (2*a**2*(1 - e**2))**2 * np.cos(i)

def d_aper(draan: FLOAT, i: FLOAT) -> FLOAT: # Argument of the periapsis
    return draan / np.cos(i) * (5/2 * np.sin(i)**2 - 2)


# Delta v requirements of maneuvers
def dvs_rendezvous(a0: FLOAT, at: FLOAT, e0: FLOAT, et: FLOAT, i0: FLOAT, it: FLOAT, raan0: FLOAT, raant: FLOAT, aper0: FLOAT, apert: FLOAT) -> tuple[FLOAT, FLOAT, FLOAT, FLOAT, FLOAT]:
    """
    Returns required maneuver Delta v:s for a rendezvous from an orbit with inital parameters (a0, e0, i0, raan0, aper0) to a target orbit (at, et, it, raant, apert).

    Params:
        a0 (FLOAT): - Initial semi-major axis
        at (FLOAT): - Target semi-major axis
        e0 (FLOAT): - Initial eccentricity
        et (FLOAT): - Target eccentricity
        i0 (FLOAT): - Initial inclination
        it (FLOAT): - Target inclination
        raan0 (FLOAT): - Initial RAAN
        raant (FLOAT): - Target RAAN
        aper0 (FLOAT): - Initial argument of periapsis
        apert (FLOAT): - Target argument of periapsis
    Returns:
        tuple (FLOAT): Tuple of dvs for the maneuver (dva, dve, dvi, dvraan, dvaper)
    """

    # Burn plan logic
    # If final orbit semi-major axis lower than current orbit: 
    #   If final orbit inclination lower than current orbit:
    #     1. Change inclination to match target
    #     2. Change RAAN to match target (add lead, drift differs until eccentricity and semi major is identical)
    #   If final orbit inclination higher than current orbit:
    #     1. Change RAAN to match target (add lead, drift differs until inclination, eccentricity and semi major is identical))
    #     2. Change inclination to match target
    #   If final orbit eccentricity lower than final orbit:
    #     3. Change eccentricity to match target
    #     4. Change argument of periapsis to match target (include lead, drift differs until semi-major is identical)
    #   If final orbit eccentricity higher than final orbit:
    #     3. Change argument of periapsis to match target (include lead, drift differs until eccentricity and semi-major is identical)
    #     4. Change eccentricity to match target
    #   5. Change altitude to match semi-major axis and phase (position in orbit relative target).
    # 
    # If final orbit semi-major axis higher than current orbit:
    #   1. Change altitude to match semi-major axis and phase (position in orbit relative target).
    #   # If final orbit eccentricity lower than final orbit:
    #     2. Change eccentricity to match target
    #     3. Change argument of periapsis to match target (include lead, drift differs until inclination and raan is identical)
    #   If final orbit eccentricity higher than final orbit:
    #     2. Change argument of periapsis to match target (include lead, drift differs until eccentricity, inclination and raan is identical)
    #     3. Change eccentricity to match target
    #   If final orbit inclination lower than current orbit:
    #     4. Change inclination to match target
    #     5. Change RAAN to match target
    #   If final orbit inclination higher than current orbit:
    #     4. Change RAAN to match target (add lead, drift differs until inclination is identical)
    #     5. Change inclination to match target

    if at<a0:
        if it<i0:
            dvi = dv_i(i0, it, a0)
            dvraan = dv_raan(raan0, raant, a0, it)
        else:
            dvraan = dv_raan(raan0, raant, a0, i0)
            dvi = dv_i(i0, it, a0)
        if et<e0:
            dve = dv_e(e0, et, a0)
            dvaper = dv_aper(aper0, apert, a0, et)
        else:
            dvaper = dv_aper(aper0, apert, a0, e0)
            dve = dv_e(e0, et, a0)
        dva = dv_a(a0, at)
    else:
        dva = dv_a(a0, at)
        if et<e0:
            dve = dv_e(e0, et, at)
            dvaper = dv_aper(aper0, apert, at, et)
        else:
            dvaper = dv_aper(aper0, apert, at, e0)
            dve = dv_e(e0, et, at)
        if it<i0:
            dvi = dv_i(i0, it, at)
            dvraan = dv_raan(raan0, raant, at, it)
        else:
            dvraan = dv_raan(raan0, raant, at, i0)
            dvi = dv_i(i0, it, at)

    return (dva, dve, dvi, dvraan, dvaper)

# Delta v requirements of maneuvers
def dvs_rendezvous_WIP(a0: FLOAT, e0: FLOAT, i0: FLOAT, raan0: FLOAT, aper0: FLOAT, epoch0: datetime, target: Object) -> tuple[FLOAT, FLOAT, FLOAT, FLOAT]:
    """
    Returns required maneuver Delta v:s, for a rendezvous from an orbit with inital parameters (a0, e0, i0, raan0, aper0) and insertion time epoch0 to rendezvous with a target.

    Params:
        a0 (FLOAT): - Initial semi-major axis
        e0 (FLOAT): - Initial eccentricity
        i0 (FLOAT): - Initial inclination
        raan0 (FLOAT): - Initial RAAN
        aper0 (FLOAT): - Initial Argument of periapsis
        target (Object): - Target to rendezvous with
    Returns:
        tuple (FLOAT): Tuple of dvs for the maneuver (dva, dve, dvi, dvraan)
    """

    M0 = 0 # Initial mean anomaly TODO: WHAT, HOW TO DO WITH MEAN ANOMALY? - assumes inserted into periapsis, which makes sense
        # GOAL IS TO HAVE SAME MEAN ANOMALY WHEN MATCHED ORBITS.. calculate n (FLOAT): - Mean motion (rad/s), as it changes with alt change?
        # THINK THIS OUT IN MAFYNETTI, IS ALSO EASIER TO WRITE IN REPORT THEN
        # TODO ADD TUPLE (or individually) BURN DELAYS, LOOK INTO WHERE DELAYS ARE REASONABLE (RAAN aper M drift)
        # During burns, integrate from one to the other (for linear function just halfway)

    # ALSO SET ARG OF PERIAPSIS EVENTUALLY WHERE EFFICIENT (see mafynetti)

    at = target.a
    et = target.e
    it = target.i
    raant, apert, M = target.at_epoch(epoch0)

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

    if at<a0:
        if it<i0:
            dvi = dv_i(i0, it, a0)
            dvraan = dv_raan(raan0, raant, a0, it)
        else:
            dvraan = dv_raan(raan0, raant, a0, i0)
            dvi = dv_i(i0, it, a0)
        dve = dv_e(e0, et, a0)
        dva = dv_a(a0, at)
    else:
        dva = dv_a(a0, at)
        dve = dv_e(e0, et, at)
        if it<i0:
            dvi = dv_i(i0, it, at)
            dvraan = dv_raan(raan0, raant, at, it)
        else:
            dvraan = dv_raan(raan0, raant, at, i0)
            dvi = dv_i(i0, it, at)

    return (dva, dve, dvi, dvraan)