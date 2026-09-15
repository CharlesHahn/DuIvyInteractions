# -*- coding: utf-8 -*-
"""group_identifiers - 基团识别器（GroupIdentifier 接口的实现）。"""

from .amber_ff_identifier import AmberFFGroupIdentifier

# 力场 → 识别器类（单一维护点，供 pipeline 和 DII 使用）
IDENTIFIER_CLASSES = {
    "amber": AmberFFGroupIdentifier,
}

__all__ = ["AmberFFGroupIdentifier", "IDENTIFIER_CLASSES"]
