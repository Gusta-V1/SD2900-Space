import json, math

with open('SD2600\cosmos-2251-debris.json') as f:
    data = json.load(f)

print(f"Total objects: {len(data)}")

mu = 398600.4418  # km^3/s^2
Re = 6378.137

def sma_from_mm(mm_rev_per_day):
    n = mm_rev_per_day * 2*math.pi / 86400.0  # rad/s
    a = (mu / n**2) ** (1/3)
    return a, n

for d in data:
    a, n = sma_from_mm(d['MEAN_MOTION'])
    d['sma'] = a
    d['n'] = n
    d['alt'] = a - Re

incs = [d['INCLINATION'] for d in data]
print("inclination range:", min(incs), max(incs))

# Tight cluster: near-identical inclination AND near-circular (good docking candidates)
cand = [d for d in data if 74.03 <= d['INCLINATION'] <= 74.06 and d['ECCENTRICITY'] < 0.003]
cand.sort(key=lambda d: d['RA_OF_ASC_NODE'])
print(f"\nCandidates in tight band: {len(cand)}")

# pick 8 spread across RAAN for a representative "cluster"
import numpy as np
if len(cand) >= 8:
    idxs = np.linspace(0, len(cand)-1, 8).astype(int)
    picks = [cand[i] for i in idxs]
else:
    picks = cand[:8]

print(f"\n{'NORAD':>6} {'OBJECT_ID':14} {'i (deg)':>8} {'RAAN (deg)':>10} {'e':>9} {'alt (km)':>9} {'n (rev/day)':>12}")
for d in picks:
    print(f"{d['NORAD_CAT_ID']:>6} {d['OBJECT_ID']:14} {d['INCLINATION']:8.4f} {d['RA_OF_ASC_NODE']:10.4f} {d['ECCENTRICITY']:9.5f} {d['alt']:9.2f} {d['MEAN_MOTION']:12.6f}")

with open('cluster8.json','w') as f:
    json.dump(picks, f, indent=2)