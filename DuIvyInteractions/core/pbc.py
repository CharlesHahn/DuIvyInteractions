# -*- coding: utf-8 -*-
"""PBC（周期边界）处理工具：统一盒尺寸获取与防护。

设计文档：doc/pbc_design.md
"""

from typing import Optional

import numpy as np


def get_box(ts) -> Optional[np.ndarray]:
    """取帧盒尺寸；无盒/退化盒返回 None（退化为欧氏距离，PBC 不生效）。

    Args:
        ts: MDAnalysis 时间步（Timestep），需有 dimensions 属性

    Returns:
        box 数组（[lx, ly, lz, alpha, beta, gamma]，与 ts.dimensions 一致）
        或 None（无盒）。
    """
    dims = getattr(ts, "dimensions", None)
    if dims is None or len(dims) == 0 or dims[0] <= 0:
        return None
    return dims