# -*- coding: utf-8 -*-
"""DuIvyInteraction 命令行入口。"""

import argparse

from .group_identifiers import IDENTIFIER_CLASSES


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

    return parser


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


if __name__ == "__main__":
    main()
