# -*- coding: utf-8 -*-
"""卤键 PerFrame vs TwoPass 性能对比。"""
import time
import tracemalloc
import numpy as np
from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HalogenBondDetectorPerFrame, HalogenBondDetectorTwoPass)
import MDAnalysis as mda

TPR = "Tests/test_MD_case/md.tpr"
XTC = "Tests/test_MD_case/md1ns.xtc"

reader = GmxTprReader()
sd = reader.read(TPR)
groups = AmberFFGroupIdentifier().identify(sd)
relevant = [g for g in groups
            if g.group_type in ("halogen_donor", "halogen_acceptor")]

donors = [g for g in relevant if g.group_type == "halogen_donor"]
acceptors = [g for g in relevant if g.group_type == "halogen_acceptor"]
print(f"Donors: {len(donors)}, Acceptors: {len(acceptors)}, "
      f"Candidates: {len(donors) * len(acceptors)}")
print()


def bench(label, detector_cls, groups, trajectory):
    tracemalloc.start()
    t0 = time.time()
    results = detector_cls().detect(groups, trajectory=trajectory)
    t1 = time.time()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = t1 - t0
    peak_mb = peak / 1024 / 1024
    n_pairs = results[0].n_pairs if results else 0
    print(f"  {label:<12} {elapsed:>8.2f}s   {peak_mb:>8.1f} MB   pairs={n_pairs}")
    return results


# PerFrame
u1 = mda.Universe(TPR, XTC)
print("[PerFrame]")
res_per = bench("PerFrame", HalogenBondDetectorPerFrame, relevant, u1.trajectory)

# TwoPass 分阶段
u2 = mda.Universe(TPR, XTC)
det_two = HalogenBondDetectorTwoPass()

print("\n[TwoPass 分阶段]")
tracemalloc.start()
t0 = time.time()
sparse = det_two.run_pass1(relevant, u2.trajectory)
t1 = time.time()
_, peak1 = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"  Pass1        {t1-t0:>8.2f}s   {peak1/1024/1024:>8.1f} MB   sparse_pairs={sparse.n_pairs}")

u3 = mda.Universe(TPR, XTC)
tracemalloc.start()
t0 = time.time()
res_two = det_two.run_pass2(sparse, u3.trajectory)
t1 = time.time()
_, peak2 = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"  Pass2        {t1-t0:>8.2f}s   {peak2/1024/1024:>8.1f} MB   pairs={res_two[0].n_pairs if res_two else 0}")

# TwoPass full
u4 = mda.Universe(TPR, XTC)
print(f"\n[TwoPass full]")
res_two_full = bench("TwoPass", HalogenBondDetectorTwoPass, relevant, u4.trajectory)

# 结果对比
it_per = res_per[0] if res_per else None
it_two = res_two_full[0] if res_two_full else None

print(f"\n{'='*55}")
print("[结果对比]")
if it_per and it_two:
    per_pairs = {(g1.group_id, g2.group_id) for g1, g2 in it_per.groups}
    two_pairs = {(g1.group_id, g2.group_id) for g1, g2 in it_two.groups}
    common = per_pairs & two_pairs

    print(f"  PerFrame: {len(per_pairs)} pairs")
    print(f"  TwoPass:  {len(two_pairs)} pairs")
    print(f"  Common:   {len(common)}")
    print(f"  PerFrame only: {len(per_pairs - two_pairs)}")
    print(f"  TwoPass only:  {len(two_pairs - per_pairs)}")

    per_occ = {}
    for i, (g1, g2) in enumerate(it_per.groups):
        per_occ[(g1.group_id, g2.group_id)] = float(it_per.occupancy()[i])
    two_occ = {}
    for i, (g1, g2) in enumerate(it_two.groups):
        two_occ[(g1.group_id, g2.group_id)] = float(it_two.occupancy()[i])

    occ_diffs = [abs(per_occ[k] - two_occ[k]) for k in common]
    print(f"  Occupancy max diff (common): {max(occ_diffs) if occ_diffs else 'N/A'}")
