# -*- coding: utf-8 -*-
"""生成所有相互作用检测器的 h5 结果文件。

基于真实数据（test_MD_case），运行所有 TwoPass 策略的检测器，
将结果保存到 interaction_h5data 目录。
"""

import sys
import time
from pathlib import Path
from typing import List

# 添加项目根目录到 Python 路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from DuIvyInteractions.system_readers import GmxTprReader
from DuIvyInteractions.group_identifiers import AmberFFGroupIdentifier
from DuIvyInteractions.core.datas import Group
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
from DuIvyInteractions.io import save_interactions
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


def main():
    """主函数：运行所有检测器并保存结果。"""
    print("=" * 60)
    print("相互作用检测器 h5 结果生成")
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
    
    # 5. 运行检测器并保存结果
    print("\n[4/4] 运行检测器并保存结果...")
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
        interactions = detector.detect(filtered_groups, u.trajectory)
        
        # 保存
        output_path = OUTPUT_DIR / f"{detector.name}.h5"
        save_interactions(interactions, str(output_path), compress=True)
        
        elapsed = time.time() - start
        
        # 输出统计
        if interactions:
            it = interactions[0]
            print(f"    基团对数: {it.n_pairs}")
            print(f"    帧数: {it.n_frames}")
            print(f"    文件大小: {output_path.stat().st_size / 1024:.2f} KB")
        else:
            print(f"    无相互作用")
        print(f"    耗时: {elapsed:.2f}s")
    
    print("\n" + "=" * 60)
    print("完成！")
    print("=" * 60)


if __name__ == "__main__":
    main()
