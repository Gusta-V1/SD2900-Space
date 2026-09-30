import json
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import minimize
from core import FLOAT, NUM_WORKERS, R, Object, load_objects, dv_a, dv_e, dv_i, dv_raan



# Main options
num_targets = 8         # Number of targets
suboptimal_targets = 3  # N:th best cluster

# Insertion and target orbit parameters
at_min = R + 800             # Minimum target altitude
i0_min = np.deg2rad(28.5)    # Minimum insertion inclination

# Cluster selection limits
#   Filter cluster candidates on semi-major axis, inclination and raan (estimated dv costs at a=600 km i=74 deg)
inclination_limit = np.deg2rad(5.0) # (deg) (about 3.4 km/s)
raan_limit = np.deg2rad(5.0)        # (deg) (about 3.53 km/s, orbital precession only 0.45 deg per month)
a_limit = 200                       # (km) (about 3.45 km/s when raising)

# Utilities
def filter_candidate(target: Object, base: Object, a_lim: float, i_lim: float, raan_lim: float):
    Da = np.abs(target.a - base.a)
    Di = np.abs(target.i - base.i)
    Draan = np.abs(target.raan - base.raan)
    return Da < a_lim and Di < i_lim and Draan < raan_lim


# Calculate maximum delta-v required for cluster
def cluster_dv_max(x, cluster: tuple[Object]):
    """
    Function to calculate maximum delta-v required for a cluster from given initial parameters.

    Params:
        x (np.ndarray([a, e, i, raan])): - Initial orbit parameters
        cluster (tuple[Object]): - List of targets
    """
    a0, e0, i0, raan0 = x

    dvtot=np.empty(len(cluster))
    for j in range(len(cluster)):
        obj = cluster[j]
        # TODO: Calculate for a given point in time? (add time parameter/constraint?), time in seconds since reference time(stamp) t_ref
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

# Process cluster
def _process_clusters(base_key):
    local_clusters=[]
    base = targets[base_key]
    candidates_keys = list(targets.keys() - {base_key})

    candidates = [targets[key] for key in candidates_keys if filter_candidate(targets[key], base, a_limit, inclination_limit, raan_limit)]
    if len(candidates) < num_targets:
        return local_clusters

    # Sample combinations of candidates and minimize
    cluster_list = [[base] + list(s) for s in itertools.combinations(candidates, num_targets-1)]
    for cluster in cluster_list:
        #Initial params
        a0 = np.mean([t.a for t in candidates])
        e0 = np.mean([t.e for t in candidates])
        i0 = np.mean([t.i for t in candidates])
        raan0 = np.mean([t.raan for t in candidates])
        x0 = np.array([a0, e0, i0, raan0], dtype=FLOAT)

        bounds=[
            (R, None),          # Semi-major axis
            (0, 1),             # Eccentricity
            (i0_min, np.pi/2),  # Inclination
            (0, 2*np.pi)        # RAAN
        ]
        result = minimize(cluster_dv_max, x0, args=(candidates,), method='SLSQP', bounds=bounds)

        local_clusters.append((cluster, result))

    return local_clusters


# Load data
targets = load_objects(at_min)
num_objects = len(targets)


# Main program
#   Note: That leads require thrust and weight data
#   TODO: Take orbital precession into consideration, 
#     (also create other function to build burn plan)
#     Add time dependency?
if __name__ == '__main__':
    print(f"Total objects: {num_objects}")

    clusters=[]
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        for local_clusters in executor.map(_process_clusters, targets.keys()):
            clusters.extend(local_clusters)

    # Sort results
    clusters = sorted(clusters, key=lambda c: c[1].fun)

    # Remove objects from n best cluster iterations
    suboptimal_clusters = list(clusters)
    for i in range(suboptimal_targets):
        removed = {o for o in suboptimal_clusters[0][0]}
        suboptimal_clusters = [c for c in suboptimal_clusters if removed.isdisjoint(set(c[0]))]
        suboptimal_clusters = sorted(suboptimal_clusters, key=lambda c: c[1].fun)

    print(f"Computations: {len(clusters)}")   
    print(f"Optimal cluster: {clusters[0][0]}")
    print(f"Optimal Parameters: \n\ta0 = {clusters[0][1].x[0]:.3f} km, \n\te0 = {clusters[0][1].x[1]:.8f}, \n\ti0 = {np.rad2deg(clusters[0][1].x[2]):.6f} deg, \n\tRAAN0 = {np.rad2deg(clusters[0][1].x[3]):.6f} deg")
    print(f"Minimum Total Delta V (km/s): {clusters[0][1].fun} km/s")


    # Save optimal cluster and insertion orbit
    data = {
        "cluster": [str(target.norad_id) for target in clusters[0][0]],
        "insertion_orbit": {
            "a": clusters[0][1].x[0],
            "e": clusters[0][1].x[1],
            "i": clusters[0][1].x[2],
            "raan": clusters[0][1].x[3],
        },
    }
    with open("cluster-solution.json", "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    
    if len(suboptimal_clusters) == 0:
        print(f"ERROR: Could not get {suboptimal_targets}:nd optimal cluster, try lowering 'suboptimal_targets'")
    else:
        print(80*"-")
        print(f"{suboptimal_targets}:th Optimal cluster: {suboptimal_clusters[0][0]}")
        print(f"{suboptimal_targets}:th Optimal Parameters: \n\ta0 = {suboptimal_clusters[0][1].x[0]:.3f} km, \n\te0 = {suboptimal_clusters[0][1].x[1]:.8f}, \n\ti0 = {np.rad2deg(suboptimal_clusters[0][1].x[2]):.6f} deg, \n\tRAAN0 = {np.rad2deg(suboptimal_clusters[0][1].x[3]):.6f} deg")
        print(f"{suboptimal_targets}:th Minimum Total Delta V (km/s): {suboptimal_clusters[0][1].fun} km/s")