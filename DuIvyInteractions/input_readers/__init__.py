# -*- coding: utf-8 -*-
"""input_readers - 输入文件读取器（Reader 接口的实现）。"""

from .gmx_tpr_dump_reader import GmxTprDumpReader
from .gmx_tpr_reader import GmxTprReader

__all__ = ["GmxTprDumpReader", "GmxTprReader"]
