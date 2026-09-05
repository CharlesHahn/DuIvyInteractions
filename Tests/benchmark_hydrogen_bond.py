# -*- coding: utf-8 -*-
"""氢键检测器性能对比：PerFrame vs TwoPass。

真实数据测试（不含水），比较结果一致性、内存、时间。
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


# ============================================================
# 配置
# ============================================================

TPR_FILE = Path(__file__).parent / "test_MD_case" / "md.tpr"
XTC_FILE = Path(__file__).parent / "test_MD_case" / "md1ns.xtc"

MOL_RBD = "seg_0_RBD_pro"
MOL_KRAS = "seg_2_KRAS_pro"
WATER_RESIDUES = {"SOL", "HOH", "WAT"}
PROTEIN_MOLECULES = {MOL_RBD, MOL_KRAS}


# ============================================================
# 数据准备
# ============================================================

def load_groups(with_water=False, water_fraction=1.0):
    """加载基团。with_water=True 时包含水。water_fraction 控制水分子采样比例。"""
    reader = GmxTprReader()
    system_data = reader.read(str(TPR_FILE))
    groups = AmberFFGroupIdentifier().identify(system_data)
    if with_water:
        # 蛋白 + 配体
        protein_groups = [g for g in groups
                          if g.residue_name not in WATER_RESIDUES
                          and g.molecule in PROTEIN_MOLECULES]
        # 水基团（按残基采样）
        water_groups = [g for g in groups if g.residue_name in WATER_RESIDUES]
        if water_fraction < 1.0:
            # 按残基ID采样
            water_resids = sorted(set(g.residue_id for g in water_groups))
            n_sample = max(1, int(len(water_resids) * water_fraction))
            sampled_resids = set(water_resids[:n_sample])
            water_groups = [g for g in water_groups if g.residue_id in sampled_resids]
            print(f"  水分子采样: {n_sample}/{len(water_resids)} 残基")
        return protein_groups + water_groups
    return [g for g in groups
            if g.residue_name not in WATER_RESIDUES
            and g.molecule in PROTEIN_MOLECULES]


def create_trajectory():
    """创建 MDAnalysis Universe。"""
    return mda.Universe(str(TPR_FILE), str(XTC_FILE))


# ============================================================
# 性能测量
# ============================================================

def measure(func, *args, **kwargs):
    """测量函数执行时间和内存。"""
    gc.collect()
    tracemalloc.start()

    start_time = time.perf_counter()
    result = func(*args, **kwargs)
    end_time = time.perf_counter()

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return result, end_time - start_time, peak / 1024 / 1024  # MB


# ============================================================
# 主测试
# ============================================================

def run_benchmark(groups, label, tuple_filter):
    """运行一次完整基准测试。"""
    print(f"\n{'='*60}")
    print(f"测试场景: {label}")
    print(f"{'='*60}")
    print(f"  基团数量: {len(groups)}")

    # PerFrame
    print(f"\n  [PerFrame] 运行中...")
    detector_pf = HydrogenBondDetectorPerFrame()
    u1 = create_trajectory()
    results_pf, time_pf, mem_pf = measure(
        detector_pf.detect, groups,
        trajectory=u1.trajectory, tuple_filter=tuple_filter,
    )
    n_pairs_pf = results_pf[0].n_pairs if results_pf else 0
    print(f"    耗时: {time_pf:.3f} s | 内存: {mem_pf:.2f} MB | 对数: {n_pairs_pf}")

    # TwoPass
    print(f"\n  [TwoPass] 运行中...")
    detector_tp = HydrogenBondDetectorTwoPass()
    u2 = create_trajectory()
    results_tp, time_tp, mem_tp = measure(
        detector_tp.detect, groups,
        trajectory=u2.trajectory, tuple_filter=tuple_filter,
    )
    n_pairs_tp = results_tp[0].n_pairs if results_tp else 0
    print(f"    耗时: {time_tp:.3f} s | 内存: {mem_tp:.2f} MB | 对数: {n_pairs_tp}")

    # 一致性检查
    if results_pf and results_tp:
        pairs_pf = {(g1.group_id, g2.group_id) for g1, g2 in results_pf[0].groups}
        pairs_tp = {(g1.group_id, g2.group_id) for g1, g2 in results_tp[0].groups}
        match = pairs_pf == pairs_tp
        print(f"\n  结果一致性: {'✅ 完全一致' if match else '❌ 不一致'}")

    # 比值
    speed_ratio = time_tp / time_pf if time_pf > 0 else float('inf')
    mem_ratio = mem_pf / mem_tp if mem_tp > 0 else float('inf')

    return {
        "label": label,
        "n_groups": len(groups),
        "n_pairs": n_pairs_pf,
        "time_pf": time_pf, "time_tp": time_tp,
        "mem_pf": mem_pf, "mem_tp": mem_tp,
        "speed_ratio": speed_ratio, "mem_ratio": mem_ratio,
    }


def main():
    print("=" * 60)
    print("氢键检测器性能对比：PerFrame vs TwoPass")
    print("101 帧 | RBD + KRAS 蛋白体系")
    print("=" * 60)

    tuple_filter = lambda gt: gt[0].molecule != gt[1].molecule

    # 场景 1：不含水
    groups_no_water = load_groups(with_water=False)
    r1 = run_benchmark(groups_no_water, "不含水（蛋白-蛋白）", tuple_filter)

    # 场景 2：含水（只取部分水，10%采样）
    groups_with_water = load_groups(with_water=True, water_fraction=0.1)
    r2 = run_benchmark(groups_with_water, "含水 10%采样（蛋白-蛋白-水）", tuple_filter)

    # 汇总
    print(f"\n\n{'='*60}")
    print("汇总对比")
    print(f"{'='*60}")
    print(f"\n{'指标':<18} {'不含水':>14} {'含水10%':>14}")
    print("-" * 48)
    print(f"{'基团数':<18} {r1['n_groups']:>14,} {r2['n_groups']:>14,}")
    print(f"{'氢键对数':<18} {r1['n_pairs']:>14,} {r2['n_pairs']:>14,}")
    print(f"{'PerFrame 耗时':<18} {r1['time_pf']:>13.3f}s {r2['time_pf']:>13.3f}s")
    print(f"{'TwoPass 耗时':<18} {r1['time_tp']:>13.3f}s {r2['time_tp']:>13.3f}s")
    print(f"{'PerFrame 内存':<18} {r1['mem_pf']:>13.2f}MB {r2['mem_pf']:>13.2f}MB")
    print(f"{'TwoPass 内存':<18} {r1['mem_tp']:>13.2f}MB {r2['mem_tp']:>13.2f}MB")
    print(f"{'速度比 (TP/PF)':<18} {r1['speed_ratio']:>13.1f}x {r2['speed_ratio']:>13.1f}x")
    print(f"{'内存比 (PF/TP)':<18} {r1['mem_ratio']:>13.1f}x {r2['mem_ratio']:>13.1f}x")


if __name__ == "__main__":
    main()
