# -*- coding: utf-8 -*-
"""氢键全量基准测试：含水-水氢键。

测试策略二（PerFrame）和策略三（TwoPass）在全体系下的时间和内存消耗。
"""

import gc
import time
import tracemalloc
from pathlib import Path

import numpy as np
import MDAnalysis as mda

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorPerFrame,
    HydrogenBondDetectorTwoPass,
)


TPR_FILE = Path(__file__).parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent / "test_MD_case" / "md1ns.xtc"


def load_all_groups():
    """加载全部基团（含水-水）。"""
    reader = GmxTprReader()
    system_data = reader.read(str(TPR_FILE))
    return AmberFFGroupIdentifier().identify(system_data)


def measure(func, *args, **kwargs):
    """测量执行时间和内存峰值。"""
    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    t1 = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, t1 - t0, peak / 1024 / 1024


def run_case(label, groups, tuple_filter):
    """运行一组测试。"""
    print(f"\n{'='*60}")
    print(f"场景: {label}")
    print(f"  基团数: {len(groups):,}")
    donors = sum(1 for g in groups if g.group_type == 'H_donor')
    acceptors = sum(1 for g in groups if g.group_type == 'H_acceptor')
    print(f"  H_donor: {donors:,}  H_acceptor: {acceptors:,}")
    print(f"  候选对: {donors * acceptors:,}")
    print(f"{'='*60}")

    # PerFrame
    print(f"\n  [PerFrame] 运行中...")
    u1 = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    pf_detector = HydrogenBondDetectorPerFrame()
    pf_result, pf_time, pf_mem = measure(
        pf_detector.detect, groups,
        trajectory=u1.trajectory, tuple_filter=tuple_filter,
    )
    n_pf = pf_result[0].n_pairs if pf_result else 0
    print(f"    耗时: {pf_time:.2f}s | 内存峰值: {pf_mem:.1f}MB | 对数: {n_pf:,}")

    # TwoPass
    print(f"\n  [TwoPass] 运行中...")
    u2 = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    tp_detector = HydrogenBondDetectorTwoPass()
    tp_result, tp_time, tp_mem = measure(
        tp_detector.detect, groups,
        trajectory=u2.trajectory, tuple_filter=tuple_filter,
    )
    n_tp = tp_result[0].n_pairs if tp_result else 0
    print(f"    耗时: {tp_time:.2f}s | 内存峰值: {tp_mem:.1f}MB | 对数: {n_tp:,}")

    # 比值
    speed = tp_time / pf_time if pf_time > 0 else float('inf')
    mem = pf_mem / tp_mem if tp_mem > 0 else float('inf')
    print(f"\n  速度比 (TwoPass/PerFrame): {speed:.1f}x")
    print(f"  内存比 (PerFrame/TwoPass): {mem:.1f}x")

    return {
        "label": label,
        "n_groups": len(groups),
        "n_donors": donors,
        "n_acceptors": acceptors,
        "n_pairs_pf": n_pf, "n_pairs_tp": n_tp,
        "time_pf": pf_time, "time_tp": tp_time,
        "mem_pf": pf_mem, "mem_tp": tp_mem,
    }


def main():
    print("=" * 60)
    print("氢键全量基准测试：PerFrame vs TwoPass")
    print("101 帧 | RBD + KRAS + D927 + GNP + 水 + 离子")
    print("=" * 60)

    all_groups = load_all_groups()

    # 场景 1：全体系（含水-水），过滤同分子内
    tuple_filter = lambda gt: gt[0].molecule != gt[1].molecule
    r1 = run_case("全体系（蛋白-蛋白-水-水）", all_groups, tuple_filter)

    # 汇总
    print(f"\n\n{'='*60}")
    print("汇总")
    print(f"{'='*60}")
    print(f"  基团数: {r1['n_groups']:,}")
    print(f"  H_donor: {r1['n_donors']:,}  H_acceptor: {r1['n_acceptors']:,}")
    print(f"  PerFrame 对数: {r1['n_pairs_pf']:,}")
    print(f"  TwoPass  对数: {r1['n_pairs_tp']:,}")
    print(f"  PerFrame 耗时: {r1['time_pf']:.2f}s  内存: {r1['mem_pf']:.1f}MB")
    print(f"  TwoPass  耗时: {r1['time_tp']:.2f}s  内存: {r1['mem_tp']:.1f}MB")
    if r1['time_pf'] > 0:
        print(f"  速度比: {r1['time_tp']/r1['time_pf']:.1f}x")
    if r1['mem_tp'] > 0:
        print(f"  内存比: {r1['mem_pf']/r1['mem_tp']:.1f}x")


if __name__ == "__main__":
    main()
