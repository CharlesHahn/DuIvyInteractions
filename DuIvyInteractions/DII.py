# -*- coding: utf-8 -*-
"""DuIvyInteraction 命令行入口。"""

import argparse
import os

import numpy as np

from .group_identifiers import IDENTIFIER_CLASSES

# 概览显示的 Top N 个占位率
OVERVIEW_TOP_N = 5


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="dii", description="DuIvyInteraction - 基于 MD 拓扑力场参数的相互作用判定工具")
    sub = parser.add_subparsers(dest="command", required=True)

    # run 子命令
    run_p = sub.add_parser("run", help="运行相互作用检测并保存 h5")
    run_p.add_argument("-t", "--tpr", required=True, help="GROMACS 拓扑文件")
    run_p.add_argument("-f", "--xtc", required=True, help="轨迹文件")
    run_p.add_argument("-o", "--output", required=True, help="输出目录")
    run_p.add_argument("--ff", required=True,
                       choices=list(IDENTIFIER_CLASSES),
                       help=f"力场（当前支持: {', '.join(IDENTIFIER_CLASSES)}）")
    run_p.add_argument("--interactions", default="all",
                       help="相互作用类型（逗号分隔），默认 all")
    run_p.add_argument("--strategy", default="two_pass",
                       choices=["two_pass", "per_frame", "per_tuple"],
                       help="检测策略，默认 two_pass")

    # export 子命令
    export_p = sub.add_parser("export", help="导出 h5 结果为 xvg/xpm/csv 并打印概览")
    export_p.add_argument("-i", "--input", required=True, help="h5 文件路径")
    export_p.add_argument("-o", "--output", required=True, help="输出目录（自动创建）")

    return parser


def _run_export(args) -> None:
    """导出 h5 结果为 xvg/xpm/csv，并打印概览。"""
    from .io.h5 import load_interactions
    from .io import (
        HydrogenBondExporter, SaltBridgeExporter, PiStackingExporter, PiCationExporter,
        HalogenBondExporter, HydrophobicExporter, MetalCoordinationExporter, WaterBridgeExporter,
    )

    # h5 里的 interaction_type（检测器名）→ Exporter 子类
    exporter_classes = {
        "hydrogen_bond": HydrogenBondExporter,
        "pi_stacking": PiStackingExporter,
        "salt_bridge": SaltBridgeExporter,
        "hydrophobic": HydrophobicExporter,
        "halogen_bond": HalogenBondExporter,
        "metal_coordination": MetalCoordinationExporter,
        "water_bridge": WaterBridgeExporter,
        "pi_cation": PiCationExporter,
    }

    interactions = load_interactions(args.input)
    if not interactions:
        raise SystemExit(f"h5 中无相互作用结果: {args.input}")
    it = interactions[0]

    if it.interaction_type not in exporter_classes:
        raise SystemExit(
            f"未知相互作用类型: '{it.interaction_type}'。可用: {', '.join(exporter_classes)}")
    exporter = exporter_classes[it.interaction_type]()

    _print_overview(it, exporter)

    os.makedirs(args.output, exist_ok=True)
    prefix = f"{args.output}/{it.interaction_type}"

    # xvg: 每帧活跃数 + 数值 metric 时间序列
    exporter.save_xvg_count(it, f"{prefix}_count.xvg")
    for metric, arr in it.metrics.items():
        if arr.dtype.kind in ('U', 'S', 'O'):
            continue
        exporter.save_xvg(it, metric, f"{prefix}_{metric}.xvg")

    # xpm: existence 热力图；π-stacking 额外输出类型图
    exporter.save_xpm(it, f"{prefix}_existence.xpm")
    if it.interaction_type == "pi_stacking":
        exporter.save_xpm_stacking_type(it, f"{prefix}_type.xpm")

    # csv: 每对汇总
    exporter.to_csv_summary(it, f"{prefix}_summary.csv")


def _print_overview(it, exporter) -> None:
    """打印概览信息。"""
    occ = it.occupancy()
    order = np.argsort(occ)[::-1][:OVERVIEW_TOP_N]
    print(f"\n===== {exporter.name} 概览 =====")
    print(f"类型:     {it.interaction_type}")
    print(f"基团对数: {it.n_pairs}")
    print(f"帧数:     {it.n_frames}")
    print(f"时间范围: {it.times[0]:.1f} ~ {it.times[-1]:.1f} ps")
    print(f"\nTop {OVERVIEW_TOP_N} 占位率:")
    for i, idx in enumerate(order, 1):
        print(f"  {i}. {exporter.get_pair_label(it, int(idx))}  {occ[idx]:.1%}")


def main():
    """命令行入口。"""
    args = build_parser().parse_args()

    if args.command == "run":
        from .pipeline import Pipeline, ALL_INTERACTIONS
        if args.interactions.lower() == "all":
            interactions = None
        else:
            interactions = [s.strip().lower() for s in args.interactions.split(",")]
            if "all" in interactions:
                raise SystemExit("'all' 不能与其他类型混用")
            for name in interactions:
                if name not in ALL_INTERACTIONS:
                    raise SystemExit(
                        f"未知相互作用类型: '{name}'。可用: {', '.join(ALL_INTERACTIONS)}")
        Pipeline(args.ff, args.strategy).run(
            args.tpr, args.xtc, args.output, interactions)
    elif args.command == "export":
        _run_export(args)


if __name__ == "__main__":
    main()
