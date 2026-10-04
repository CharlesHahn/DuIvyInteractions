# -*- coding: utf-8 -*-
"""group_identifiers - 基团识别器（GroupIdentifier 接口的实现）。"""

from .amber_ff_identifier import AmberFFGroupIdentifier
from .gromos_ff_identifier import GromosFFGroupIdentifier
from .charmm_ff_identifier import CharmmFFGroupIdentifier

# 力场 → 识别器类（单一维护点，供 pipeline 和 DII 使用）
IDENTIFIER_CLASSES = {
    "amber": AmberFFGroupIdentifier,
    "gromos": GromosFFGroupIdentifier,   # GROMOS 53A6 / 54A7 家族
    "charmm": CharmmFFGroupIdentifier,   # CHARMM36 / C36m 家族
}

__all__ = ["AmberFFGroupIdentifier", "GromosFFGroupIdentifier",
           "CharmmFFGroupIdentifier", "IDENTIFIER_CLASSES"]
