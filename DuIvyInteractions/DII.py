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
        prog="dii", description="DuIvyInteraction - molecular interaction detection based on MD topology force-field parameters")
    sub = parser.add_subparsers(dest="command", required=True)

    # run 子命令
    run_p = sub.add_parser("run", help="run interaction detection and save h5")
    run_p.add_argument("-t", "--tpr", required=True, help="GROMACS topology file")
    run_p.add_argument("-f", "--xtc", required=True, help="trajectory file")
    run_p.add_argument("-o", "--output", required=True, help="output directory")
    run_p.add_argument("--ff", required=True,
                       choices=list(IDENTIFIER_CLASSES),
                       help=f"force field (supported: {', '.join(IDENTIFIER_CLASSES)})")
    run_p.add_argument("--interactions", default="all",
                       help="interaction types (comma-separated), default all")
    run_p.add_argument("--strategy", default="two_pass",
                       choices=["two_pass", "per_frame", "per_tuple"],
                       help="detection strategy, default two_pass")

    # export 子命令
    export_p = sub.add_parser("export", help="export h5 results to xvg/xpm/csv and print overview")
    export_p.add_argument("-i", "--input", required=True, help="h5 file path")
    export_p.add_argument("-o", "--output", required=True, help="output directory (auto-created)")

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

    try:
        interactions = load_interactions(args.input)
    except OSError:
        raise SystemExit(
            f"cannot read h5 file: {args.input} (file may be corrupted or missing)")
    if not interactions:
        raise SystemExit(f"no interaction results in h5: {args.input}")

    os.makedirs(args.output, exist_ok=True)

    # 同类型重复时追加序号（如 salt_bridge_2_*）
    seen_types = {}
    for it in interactions:
        if it.interaction_type not in exporter_classes:
            raise SystemExit(
                f"unknown interaction type: '{it.interaction_type}'. Available: {', '.join(exporter_classes)}")
        exporter = exporter_classes[it.interaction_type]()

        # 空数据防护：0 对或 0 帧 → 只提示，跳过导出
        if it.n_pairs == 0 or it.n_frames == 0:
            print(f"[skip] {it.interaction_type}: nothing to export"
                  f" (pairs={it.n_pairs}, frames={it.n_frames})")
            continue

        _print_overview(it, exporter)

        seen_types[it.interaction_type] = seen_types.get(it.interaction_type, 0) + 1
        n = seen_types[it.interaction_type]
        prefix = f"{args.output}/{it.interaction_type}"
        if n > 1:
            prefix += f"_{n}"

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
    print(f"\n===== {exporter.name} overview =====")
    print(f"Type:     {it.interaction_type}")
    print(f"Pairs:    {it.n_pairs}")
    print(f"Frames:   {it.n_frames}")
    if it.n_frames > 0:
        print(f"Time range: {it.times[0]:.1f} ~ {it.times[-1]:.1f} ps")
    print(f"\nTop {OVERVIEW_TOP_N} occupancies:")
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
                raise SystemExit("'all' cannot be mixed with other types")
            for name in interactions:
                if name not in ALL_INTERACTIONS:
                    raise SystemExit(
                        f"unknown interaction type: '{name}'. Available: {', '.join(ALL_INTERACTIONS)}")
        Pipeline(args.ff, args.strategy).run(
            args.tpr, args.xtc, args.output, interactions)
    elif args.command == "export":
        _run_export(args)


if __name__ == "__main__":
    main()
