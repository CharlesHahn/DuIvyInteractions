# -*- coding: utf-8 -*-
"""web - 本地 Web 传输层（FastAPI 薄胶水；逻辑在 visualizers/）。"""

from .app import create_app

__all__ = ["create_app"]
