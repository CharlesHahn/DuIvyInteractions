# -*- coding: utf-8 -*-
"""DuIvyInteraction 命令行入口。"""

import argparse
import ipaddress
import os
import secrets
import socket

import numpy as np

from .group_identifiers import IDENTIFIER_CLASSES

# 概览显示的 Top N 个占位率
OVERVIEW_TOP_N = 5
# web：token 字节数、通配绑定地址、回环主机名、UDP 探测目标
WEB_TOKEN_BYTES = 16
WILDCARD_HOSTS = ("0.0.0.0", "::")
LOOPBACK_HOSTS = ("localhost", "::1")
LOCAL_HOST = "127.0.0.1"
IPV4_PROBE_TARGET = ("8.8.8.8", 80)


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
    run_p.add_argument("--select", action="append", nargs="+",
                       metavar="SELEXP", default=None,
                       help="MDAnalysis selection string defining one participant set "
                            "(repeatable or space-separated multi-value; "
                            "each group must be fully covered by one select)")
    run_p.add_argument("--pairing", default="any", choices=["any", "between"],
                       help="pairing rule between selects: any (default) or between "
                            "(keep only cross-select pairs; requires >=2 --select)")
    run_p.add_argument("--mark-redundant", action="store_true",
                       help="mark hydrophobic pair-frames overlapping pi-stacking "
                            "(mark only, not remove; auto-runs pi_stacking if needed)")

    # web 子命令
    web_p = sub.add_parser("web", help="launch local web UI to explore h5 results")
    web_p.add_argument("--port", type=int, default=8888, help="listen port, default 8888 (auto-advance if taken)")
    web_p.add_argument("--host", default="127.0.0.1", help="bind host, default 127.0.0.1")
    web_p.add_argument("--token", default="",
                       help="access token; auto-generated when binding to a "
                            "non-loopback host")
    web_p.add_argument("--root", default=None,
                       help="restrict file browsing to DIR (default: current "
                            "directory for non-loopback bindings; unrestricted "
                            "on loopback)")

    # export 子命令
    export_p = sub.add_parser("export", help="export h5 results to xvg/xpm/csv and print overview")
    export_p.add_argument("-i", "--input", required=True, help="h5 file path")
    export_p.add_argument("-o", "--output", required=True, help="output directory (auto-created)")

    return parser


def _run_export(args) -> None:
    """导出 h5 结果为 xvg/xpm/csv，并打印概览。"""
    from .io.h5 import load_interactions
    from .io import EXPORTER_CLASSES

    # h5 里的 interaction_type（检测器名）→ Exporter 子类（注册表见 io/__init__.py）
    exporter_classes = EXPORTER_CLASSES

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
    for name, flags in it.flags.items():
        active_flags = flags & it.existence
        n_marked = int(active_flags.any(axis=1).sum())
        n_frames = int(active_flags.sum())
        print(f"Flags:    {name}: {n_marked}/{it.n_pairs} pairs "
              f"({n_frames} pair-frames)")
    print(f"\nTop {OVERVIEW_TOP_N} occupancies:")
    for i, idx in enumerate(order, 1):
        print(f"  {i}. {exporter.get_pair_label(it, int(idx))}  {occ[idx]:.1%}")


def _free_port(host: str, port: int, tries: int = 10) -> int:
    """端口占用时 +1 顺延（最多 tries 个），返回首个可用端口。"""
    import socket
    for candidate in range(port, port + tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, candidate))
                return candidate
            except OSError:
                continue
    raise SystemExit(f"error: no free port in {port}..{port + tries - 1}")


def _is_loopback_host(host: str) -> bool:
    """判断主机名是否为回环地址（localhost/::1/127.0.0.0/8）。"""
    name = (host or "").strip().lower()
    if name in LOOPBACK_HOSTS:
        return True
    try:
        return ipaddress.ip_address(name).is_loopback
    except ValueError:
        return False


def _browse_root(root_arg: str, is_local: bool) -> str:
    """浏览根目录：显式 --root 优先；本地默认不限制；远程默认当前目录。"""
    if root_arg:
        return os.path.abspath(root_arg)
    return "" if is_local else os.getcwd()


def _is_usable_ipv4(ip: str) -> bool:
    """判断字符串是否为可对外使用的 IPv4（排除回环与 0.0.0.0）。"""
    try:
        addr = ipaddress.IPv4Address(ip)
    except ipaddress.AddressValueError:
        return False
    return not (addr.is_loopback or addr.is_unspecified)


def _probe_ipv4() -> str:
    """通过 UDP 探测本机出网 IPv4（失败返回空串）。"""
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
            sock.connect(IPV4_PROBE_TARGET)
            return sock.getsockname()[0]
    except Exception:
        return ""


def _local_ipv4s() -> list[str]:
    """本机非回环 IPv4 地址列表（去重、稳定排序）。"""
    found: set[str] = set()
    try:
        found.update(socket.gethostbyname_ex(socket.gethostname())[2])
    except Exception:
        pass
    found.add(_probe_ipv4())
    return sorted(ip for ip in found if _is_usable_ipv4(ip))


def _url(host: str, port: int) -> str:
    """拼接 http URL。"""
    return "http://%s:%d" % (host, port)


def _lan_urls(host: str, port: int, ips: list[str]) -> list[str]:
    """按绑定主机生成可访问 URL 列表（回环在前，去重保序）。"""
    if host in WILDCARD_HOSTS:
        urls = [_url(LOCAL_HOST, port)] + [_url(ip, port) for ip in ips]
    elif _is_loopback_host(host):
        urls = [_url(LOCAL_HOST, port)]
    else:
        urls = [_url(host, port)]
    return list(dict.fromkeys(urls))


def _print_web_banner(urls: list[str], token: str, host: str,
                      is_local: bool, root: str) -> None:
    """打印访问地址、root 范围、token、cookie 提示与绑定警告。"""
    suffix = f"/?token={token}" if token and not is_local else ""
    for url in urls:
        print(f"[INFO] dii web serving at {url}{suffix}")
    if root:
        print(f"[INFO] file access restricted to: {root}")
    else:
        print("[INFO] file access: unrestricted (loopback only)")
    if token:
        print(f"[INFO] token: {token}")
        print("[INFO] open the URL above once; a cookie keeps you "
              "signed in afterwards")
    if not is_local:
        print(f"[WARN] binding to {host}: reachable from the network; "
              "token protects access")
    print("[INFO] press Ctrl+C to stop")


def _allowed_hosts(host: str) -> list[str]:
    """Web Host 头白名单：回环地址 + 本机 IPv4（+ 显式绑定主机）。"""
    hosts = [LOCAL_HOST, "localhost", "::1"] + _local_ipv4s()
    if host and host not in WILDCARD_HOSTS:
        hosts.append(host)
    return list(dict.fromkeys(hosts))


def _open_browser(url: str) -> None:
    """延时用默认浏览器打开访问地址。"""
    import threading
    import webbrowser
    threading.Timer(0.5, webbrowser.open, args=(url,)).start()


def _run_web(args: argparse.Namespace) -> None:
    """启动 Web 服务并打印访问地址。"""
    from .web.app import create_app
    from .web.picker import picker_available
    port = _free_port(args.host, args.port)
    is_local = _is_loopback_host(args.host)
    token = args.token or ("" if is_local else secrets.token_urlsafe(WEB_TOKEN_BYTES))
    browse_mode = "picker" if (is_local and picker_available()) else "server"
    root = _browse_root(args.root, is_local)
    allowed_hosts = _allowed_hosts(args.host)
    urls = _lan_urls(args.host, port, _local_ipv4s())
    _print_web_banner(urls, token, args.host, is_local, root)
    if is_local:
        _open_browser(urls[0])
    import uvicorn
    uvicorn.run(create_app(token=token, browse_mode=browse_mode,
                           root=root, allowed_hosts=allowed_hosts),
                host=args.host, port=port, log_level="warning")


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
        try:
            # append+nargs='+' → [[sel, ...], ...] 拍平为 [sel, ...]
            selects = (sum(args.select, []) if args.select else None)
            Pipeline(args.ff, args.strategy).run(
                args.tpr, args.xtc, args.output, interactions,
                select=selects, pairing=args.pairing,
                mark_redundant=args.mark_redundant)
        except ValueError as e:
            # 用户输入/选择歧义 → 友好报错（不暴露 traceback）
            raise SystemExit(f"error: {e}")
    elif args.command == "web":
        _run_web(args)
    elif args.command == "export":
        _run_export(args)


if __name__ == "__main__":
    main()
