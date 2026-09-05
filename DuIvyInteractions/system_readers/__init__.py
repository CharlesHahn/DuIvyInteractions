# -*- coding: utf-8 -*-
"""system_readers - 系统数据读取器（Reader 接口的实现）。"""

from .gmx_tpr_dump_reader import GmxTprDumpReader
from .gmx_tpr_reader import GmxTprReader

__all__ = ["GmxTprDumpReader", "GmxTprReader"]
