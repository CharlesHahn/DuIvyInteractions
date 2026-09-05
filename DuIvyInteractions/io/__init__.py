# -*- coding: utf-8 -*-
"""io - 结果文件读写（Interaction 数据的序列化/反序列化）。"""

from .h5 import save_interactions, load_interactions

__all__ = ["save_interactions", "load_interactions"]
