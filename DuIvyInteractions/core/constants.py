# -*- coding: utf-8 -*-
"""核心常量定义。"""

# 基团类型
GROUP_TYPES = frozenset({
    "H_donor", "H_acceptor",           # H键
    "aromatic_ring",                    # π 相关
    "charged_positive", "charged_negative",  # 电荷
    "halogen", "metal", "water", "hydrophobic",  # 其他
})

# 键类型
BOND_TYPES = frozenset({
    "single", "double", "triple", "aromatic",  # 化学键
    "constrained", "virtual",                   # 特殊
})
