# -*- coding: utf-8 -*-
"""生成所有相互作用检测器的 h5 结果文件，并验证往返一致性。

基于真实数据（test_MD_case），运行所有 TwoPass 策略的检测器，
将结果保存到 interaction_h5data 目录，并验证 h5 存储的无损性。
"""

import sys
import time
import numpy as np
from pathlib import Path
from typing import List

# 添加项目根目录到 Python 路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.core.datas import Group, Interaction
from DuIvyInteractions.interaction_detectors import (
    HydrogenBondDetectorTwoPass,
    PiStackingDetectorTwoPass,
    SaltBridgeDetectorTwoPass,
    HydrophobicDetectorTwoPass,
    HalogenBondDetectorTwoPass,
    MetalCoordinationDetectorTwoPass,
    WaterBridgeDetectorTwoPass,
    PiCationDetectorTwoPass,
)
from DuIvyInteractions.io import save_interactions, load_interactions
import MDAnalysis as mda


# 路径配置
TPR_FILE = project_root / "Tests" / "test_MD_case" / "md.tpr"
XTC_FILE = project_root / "Tests" / "test_MD_case" / "md1ns.xtc"
OUTPUT_DIR = project_root / "Tests" / "interaction_h5data"

# 水分子残基名
WATER_RESIDUES = {"SOL", "HOH", "WAT"}


def filter_groups_by_type(groups: List[Group], required_types: List[str],
                          exclude_water: bool = True) -> List[Group]:
    """按类型过滤基团，可选择排除水分子。"""
    filtered = []
    for g in groups:
        # 检查类型是否匹配
        if g.group_type not in required_types:
            continue
        # 排除水分子
        if exclude_water and g.residue_name in WATER_RESIDUES:
            continue
        filtered.append(g)
    return filtered


def compare_interactions(original: List[Interaction], loaded: List[Interaction],
                         name: str) -> bool:
    """比较原始数据和加载的数据是否一致。"""
    print(f"\n  验证 {name}:")
    
    # 检查数量
    if len(original) != len(loaded):
        print(f"    ❌ 数量不一致: {len(original)} vs {len(loaded)}")
        return False
    
    if len(original) == 0:
        print(f"    ✅ 无相互作用，跳过验证")
        return True
    
    all_ok = True
    for i, (orig, load) in enumerate(zip(original, loaded)):
        print(f"    Interaction {i}:")
        
        # 比较 interaction_type
        if orig.interaction_type != load.interaction_type:
            print(f"      ❌ interaction_type 不一致: '{orig.interaction_type}' vs '{load.interaction_type}'")
            all_ok = False
        else:
            print(f"      ✅ interaction_type 一致")
        
        # 比较 existence
        if not np.array_equal(orig.existence, load.existence):
            print(f"      ❌ existence 不一致")
            all_ok = False
        else:
            print(f"      ✅ existence 一致 (shape={orig.existence.shape})")
        
        # 比较 times
        if not np.array_equal(orig.times, load.times):
            print(f"      ❌ times 不一致")
            all_ok = False
        else:
            print(f"      ✅ times 一致 (len={len(orig.times)})")
        
        # 比较 metrics
        if set(orig.metrics.keys()) != set(load.metrics.keys()):
            print(f"      ❌ metrics keys 不一致")
            all_ok = False
        else:
            metrics_ok = True
            for key in orig.metrics:
                orig_val = orig.metrics[key]
                load_val = load.metrics[key]
                
                # 检查类型
                if hasattr(orig_val, 'dtype') and orig_val.dtype.kind in ('U', 'S', 'O'):
                    # 字符串类型
                    orig_list = [str(s) for s in orig_val.flatten()]
                    load_list = [str(s) for s in np.array(load_val).flatten()]
                    if orig_list != load_list:
                        print(f"      ❌ metrics['{key}'] 不一致")
                        metrics_ok = False
                else:
                    # 数值类型
                    if not np.array_equal(orig_val, load_val):
                        print(f"      ❌ metrics['{key}'] 不一致")
                        metrics_ok = False
            
            if metrics_ok:
                print(f"      ✅ metrics 一致 (keys={list(orig.metrics.keys())})")
            else:
                all_ok = False
        
        # 比较 groups
        if len(orig.groups) != len(load.groups):
            print(f"      ❌ groups 数量不一致: {len(orig.groups)} vs {len(load.groups)}")
            all_ok = False
        else:
            groups_ok = True
            for j, (orig_tuple, load_tuple) in enumerate(zip(orig.groups, load.groups)):
                if len(orig_tuple) != len(load_tuple):
                    print(f"      ❌ groups[{j}] tuple 长度不一致")
                    groups_ok = False
                    break
                
                for k, (orig_group, load_group) in enumerate(zip(orig_tuple, load_tuple)):
                    # 比较 group 基本信息
                    if (orig_group.group_id != load_group.group_id or
                        orig_group.group_type != load_group.group_type or
                        orig_group.molecule != load_group.molecule or
                        orig_group.residue_name != load_group.residue_name or
                        orig_group.residue_id != load_group.residue_id):
                        print(f"      ❌ groups[{j}][{k}] 基本信息不一致")
                        groups_ok = False
                        break
                    
                    # 比较 atoms
                    if len(orig_group.atoms) != len(load_group.atoms):
                        print(f"      ❌ groups[{j}][{k}] atoms 数量不一致")
                        groups_ok = False
                        break
                    
                    for l, (orig_atom, load_atom) in enumerate(zip(orig_group.atoms, load_group.atoms)):
                        if (orig_atom.atom_global_idx != load_atom.atom_global_idx or
                            orig_atom.atom_idx_in_residue != load_atom.atom_idx_in_residue or
                            orig_atom.atom_name != load_atom.atom_name or
                            orig_atom.atom_type != load_atom.atom_type or
                            orig_atom.atom_element != load_atom.atom_element or
                            abs(orig_atom.atom_charge - load_atom.atom_charge) > 1e-10 or
                            abs(orig_atom.atom_mass - load_atom.atom_mass) > 1e-10):
                            print(f"      ❌ groups[{j}][{k}].atoms[{l}] 不一致")
                            groups_ok = False
                            break
                    
                    if not groups_ok:
                        break
                
                if not groups_ok:
                    break
            
            if groups_ok:
                print(f"      ✅ groups 一致 (n_pairs={len(orig.groups)})")
            else:
                all_ok = False
        
        # 比较 metadata
        for j, (orig_tuple, load_tuple) in enumerate(zip(orig.groups, load.groups)):
            for k, (orig_group, load_group) in enumerate(zip(orig_tuple, load_tuple)):
                if orig_group.metadata != load_group.metadata:
                    print(f"      ❌ groups[{j}][{k}].metadata 不一致")
                    all_ok = False
                    break
            if not all_ok:
                break
        
        if all_ok:
            print(f"      ✅ metadata 一致")
    
    return all_ok


def main():
    """主函数：运行所有检测器并保存结果，验证往返一致性。"""
    print("=" * 60)
    print("相互作用检测器 h5 结果生成与验证")
    print("=" * 60)
    
    # 1. 读取 SystemData
    print("\n[1/4] 读取 SystemData...")
    start = time.time()
    reader = GmxTprReader()
    system_data = reader.read(str(TPR_FILE))
    print(f"  完成，耗时 {time.time() - start:.2f}s")
    print(f"  残基数量: {system_data.n_residues}")
    
    # 2. 识别基团
    print("\n[2/4] 识别基团...")
    start = time.time()
    identifier = AmberFFGroupIdentifier()
    all_groups = identifier.identify(system_data)
    print(f"  完成，耗时 {time.time() - start:.2f}s")
    print(f"  基团数量: {len(all_groups)}")
    
    # 3. 加载轨迹
    print("\n[3/4] 加载轨迹...")
    start = time.time()
    u = mda.Universe(str(TPR_FILE), str(XTC_FILE))
    print(f"  完成，耗时 {time.time() - start:.2f}s")
    print(f"  帧数: {u.trajectory.n_frames}")
    
    # 4. 定义检测器和对应的基团类型
    detector_configs = [
        {
            "detector": HydrogenBondDetectorTwoPass(),
            "required_types": ["H_donor", "H_acceptor"],
            "exclude_water": True,  # 排除水分子，避免组合爆炸
        },
        {
            "detector": PiStackingDetectorTwoPass(),
            "required_types": ["aromatic_ring"],
            "exclude_water": False,  # 无水分子
        },
        {
            "detector": SaltBridgeDetectorTwoPass(),
            "required_types": ["charged_positive", "charged_negative"],
            "exclude_water": False,  # 无水分子
        },
        {
            "detector": HydrophobicDetectorTwoPass(),
            "required_types": ["hydrophobic"],
            "exclude_water": False,  # 无水分子
        },
        {
            "detector": HalogenBondDetectorTwoPass(),
            "required_types": ["halogen_donor", "halogen_acceptor"],
            "exclude_water": False,  # 无水分子
        },
        {
            "detector": MetalCoordinationDetectorTwoPass(),
            "required_types": ["metal", "metal_binding"],
            "exclude_water": False,  # 无水分子
        },
        {
            "detector": WaterBridgeDetectorTwoPass(),
            "required_types": ["H_donor", "H_acceptor", "water"],
            "exclude_water": False,  # 需要水分子
        },
        {
            "detector": PiCationDetectorTwoPass(),
            "required_types": ["aromatic_ring", "charged_positive"],
            "exclude_water": False,  # 无水分子
        },
    ]
    
    # 5. 运行检测器、保存结果并验证
    print("\n[4/4] 运行检测器、保存结果并验证...")
    results = []
    for config in detector_configs:
        detector = config["detector"]
        required_types = config["required_types"]
        exclude_water = config["exclude_water"]
        
        print(f"\n  运行 {detector.name}...")
        print(f"    需要基团类型: {required_types}")
        print(f"    排除水分子: {exclude_water}")
        
        # 过滤基团
        filtered_groups = filter_groups_by_type(all_groups, required_types, exclude_water)
        print(f"    过滤后基团数量: {len(filtered_groups)}")
        
        start = time.time()
        
        # 检测
        interactions_original = detector.detect(filtered_groups, u.trajectory)
        
        # 保存
        output_path = OUTPUT_DIR / f"{detector.name}.h5"
        save_interactions(interactions_original, str(output_path), compress=True)
        
        # 加载
        interactions_loaded = load_interactions(str(output_path))
        
        elapsed = time.time() - start
        
        # 验证
        ok = compare_interactions(interactions_original, interactions_loaded, detector.name)
        results.append((detector.name, ok))
        
        # 输出统计
        if interactions_original:
            it = interactions_original[0]
            print(f"    基团对数: {it.n_pairs}")
            print(f"    帧数: {it.n_frames}")
            print(f"    文件大小: {output_path.stat().st_size / 1024:.2f} KB")
        else:
            print(f"    无相互作用")
        print(f"    耗时: {elapsed:.2f}s")
    
    # 6. 总结
    print("\n" + "=" * 60)
    print("验证结果总结")
    print("=" * 60)
    all_ok = True
    for name, ok in results:
        status = "✅ 通过" if ok else "❌ 失败"
        print(f"  {name}: {status}")
        if not ok:
            all_ok = False
    
    if all_ok:
        print("\n✅ 所有检测器的 h5 存储往返验证通过！")
    else:
        print("\n❌ 部分检测器验证失败！")
    
    print("=" * 60)


if __name__ == "__main__":
    main()
