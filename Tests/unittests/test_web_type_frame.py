# -*- coding: utf-8 -*-
"""类型页帧选择器：/type/<name>?frame=N 后端钳制 + 模板 + matrix.js 静态断言。

方案 A：删除独立帧滑块，矩阵列头即帧选择器；帧号进入 URL 与矩阵配置。
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

from DuIvyInteractions.core.datas import AtomData, Group, Interaction
from DuIvyInteractions.io.h5 import save_interactions

FRAMES = 12
REPO_ROOT = Path(__file__).resolve().parents[2]
STATIC_DIR = REPO_ROOT / "DuIvyInteractions" / "web" / "static"
TEMPLATE_DIR = REPO_ROOT / "DuIvyInteractions" / "web" / "templates"
MATRIX_CONFIG_RE = re.compile(
    r'<script id="matrixConfig" type="application/json">(.*?)</script>', re.S)
FRAME_INPUT_RE = re.compile(
    r'<input[^>]*id="frameInput"[^>]*value="([^"]*)"')


@pytest.fixture(autouse=True)
def isolate_recent(monkeypatch, tmp_path):
    monkeypatch.setattr("DuIvyInteractions.visualizers.data.RECENT_PATH",
                        tmp_path / "recent.json")


def _atom(gidx: int, name: str, element: str) -> AtomData:
    return AtomData(atom_global_idx=gidx, atom_idx_in_residue=0, atom_name=name,
                    atom_type="", atom_element=element, atom_charge=0.0,
                    atom_mass=0.0)


def _group(gid: int, gtype: str, gidx: int, res: str, rid: int,
           with_h: bool = False) -> Group:
    atoms = [_atom(gidx, "C", "C")]
    if with_h:
        atoms.append(_atom(gidx + 1, "H", "H"))
    return Group(group_id=gid, group_type=gtype, molecule="test",
                 residue_name=res, residue_id=rid, atoms=atoms)


def _interaction(frames: int = FRAMES) -> Interaction:
    existence = np.zeros((2, frames), dtype=bool)
    if frames:
        existence[:, 0] = True
        existence[0, min(4, frames - 1)] = True
    metrics = {"distance": np.where(existence, 3.0, np.nan)}
    groups = [(_group(2 * i, "H_donor", 10 * i, "LYS", i, with_h=True),
               _group(2 * i + 1, "H_acceptor", 10 * i + 1, "ASP", 100 + i))
              for i in range(2)]
    return Interaction(interaction_type="hydrogen_bond", groups=groups,
                       existence=existence, metrics=metrics,
                       times=np.arange(frames, dtype=float) * 10.0)


def _client_for(directory: Path):
    from fastapi.testclient import TestClient
    from DuIvyInteractions.web.app import create_app
    client = TestClient(create_app())
    client.post("/select", data={"path": str(directory)})
    return client


@pytest.fixture
def h5_dir(tmp_path):
    directory = tmp_path / "out"
    directory.mkdir()
    save_interactions([_interaction()], str(directory / "hydrogen_bond.h5"))
    return directory


@pytest.fixture
def client(h5_dir):
    return _client_for(h5_dir)


def _type_page(client, query: str = ""):
    return client.get("/type/hydrogen_bond" + query)


def _html_frame_input(html: str) -> int:
    match = FRAME_INPUT_RE.search(html)
    assert match, "frameInput element missing from type page"
    return int(match.group(1))


def _html_matrix_cfg(html: str) -> dict:
    match = MATRIX_CONFIG_RE.search(html)
    assert match, "matrixConfig script missing from type page"
    return json.loads(match.group(1))


def _js_source() -> str:
    return (STATIC_DIR / "matrix.js").read_text(encoding="utf-8")


# ============================================================
# 后端：frame 查询参数与钳制
# ============================================================

class TestFrameQuery:

    def test_default_frame_is_zero(self, client):
        response = _type_page(client)
        assert response.status_code == 200
        assert _html_frame_input(response.text) == 0
        assert _html_matrix_cfg(response.text)["frame"] == 0

    def test_explicit_frame_is_rendered(self, client):
        response = _type_page(client, "?frame=5")
        assert response.status_code == 200
        assert _html_frame_input(response.text) == 5
        assert _html_matrix_cfg(response.text)["frame"] == 5

    def test_frame_above_range_clamped_to_last(self, client):
        response = _type_page(client, "?frame=99999")
        assert response.status_code == 200
        assert _html_frame_input(response.text) == FRAMES - 1
        assert _html_matrix_cfg(response.text)["frame"] == FRAMES - 1

    def test_negative_frame_clamped_to_zero(self, client):
        response = _type_page(client, "?frame=-3")
        assert response.status_code == 200
        assert _html_frame_input(response.text) == 0
        assert _html_matrix_cfg(response.text)["frame"] == 0

    def test_non_integer_frame_rejected_with_422(self, client):
        """非整数值走 FastAPI 参数校验 → 422。

        与既有 `page: int = 1` 保持同一约定，不改成 str 自行解析，
        避免同页参数出现两套校验方式。
        """
        assert _type_page(client, "?frame=abc").status_code == 422

    def test_frame_survives_other_query_params(self, client):
        response = _type_page(
            client, "?frame=7&search=LYS&min_occupancy=0.2"
                    "&sort_by=index&sort_dir=asc&page=1")
        assert response.status_code == 200
        cfg = _html_matrix_cfg(response.text)
        assert _html_frame_input(response.text) == 7
        assert cfg["frame"] == 7
        assert cfg["search"] == "LYS"
        assert cfg["min_occupancy"] == 0.2
        assert cfg["sort_by"] == "index" and cfg["sort_dir"] == "asc"

    def test_zero_frame_h5_clamps_to_zero(self, tmp_path):
        directory = tmp_path / "empty_frames"
        directory.mkdir()
        save_interactions([_interaction(frames=0)],
                          str(directory / "hydrogen_bond.h5"))
        response = _type_page(_client_for(directory), "?frame=99")
        assert response.status_code == 200
        assert _html_frame_input(response.text) == 0
        assert _html_matrix_cfg(response.text)["frame"] == 0


# ============================================================
# 模板：控件 4 → 1，帧检查器紧跟矩阵
# ============================================================

class TestTypeTemplate:

    def test_frame_slider_removed(self, client):
        html = _type_page(client).text
        assert "frameSlider" not in html
        assert "inspect frame" not in html

    def test_frame_inspector_elements_present(self, client):
        html = _type_page(client).text
        assert 'id="frameInput"' in html
        assert 'id="framePrev"' in html and 'id="frameNext"' in html
        assert 'id="frameInfo"' in html
        assert 'id="frameBody"' in html

    def test_frame_inspector_follows_matrix(self, client):
        html = _type_page(client).text
        assert html.index('id="matrixCanvas"') < html.index('id="frameBody"')

    def test_matrix_hint_text_present(self, client):
        html = _type_page(client).text
        assert ("click a column header to inspect that frame"
                in html)
        assert "click a cell to open the pair at that frame" in html

    def test_matrix_window_controls_preserved(self, client):
        html = _type_page(client).text
        for element_id in ("viewStart", "viewEnd", "applyWindow"):
            assert 'id="%s"' % element_id in html

    def test_matrix_config_carries_frame(self, client):
        cfg = _html_matrix_cfg(_type_page(client).text)
        assert cfg["frame"] == 0
        assert cfg["n_frames"] == FRAMES
        assert cfg["type"] == "hydrogen_bond"

    def test_template_declares_frame_input(self):
        template = (TEMPLATE_DIR / "type.html").read_text(encoding="utf-8")
        assert '<input id="frameInput"' in template
        assert '"frame": frame' in template


# ============================================================
# matrix.js：滑块可选、列头选帧、格点带帧跳转、URL 同步
# ============================================================

class TestMatrixJs:

    def test_slider_is_optional(self):
        source = _js_source()
        assert "if (slider) slider.value = frame;" in source
        assert "if (slider) {" in source

    def test_pair_link_carries_frame(self):
        source = _js_source()
        assert '"/pair/" + encodeURIComponent(cfg.type)' in source
        assert '"?frame=" + frame' in source

    def test_url_sync_uses_replace_state(self):
        source = _js_source()
        assert "history.replaceState" in source
        assert 'url.searchParams.set("frame", frame)' in source

    def test_initial_frame_comes_from_config(self):
        source = _js_source()
        assert "cfg.frame" in source
        assert "loadFrame(initialFrame())" in source

    def test_frame_input_is_synced(self):
        source = _js_source()
        assert "if (frameInput) frameInput.value = frame;" in source

    def test_axis_acts_as_column_header_hot_zone(self):
        source = _js_source()
        assert "axisCol" in source
        assert "var inAxis = rows.length > 0 && y >= plotHeight()" in source
        assert "AXIS_H" in source

    def _open_pair_body(self, source):
        match = re.search(r"  function openPair\(.*?\n  \}\n", source, re.S)
        assert match
        return match.group(0)

    def test_pair_navigation_uses_pjax_router(self):
        """格点/打开 pair 走 pjax 路由，不再整页跳转。"""
        source = _js_source()
        assert "function navigateTo(url)" in source
        assert "navigateTo(url);" in self._open_pair_body(source)
        assert "window.DiiRouter.navigate(url);" in source

    def test_open_pair_does_not_assign_location(self):
        body = self._open_pair_body(_js_source())
        assert "navigateTo(url);" in body
        assert "window.location.href" not in body

    def test_navigate_to_keeps_location_fallback(self):
        source = _js_source()
        match = re.search(r"  function navigateTo\(.*?\n  \}\n", source, re.S)
        assert match
        body = match.group(0)
        assert "window.DiiRouter && typeof window.DiiRouter.navigate" in body
        assert "window.location.href = url;" in body

    def test_no_direct_location_jump_to_pair(self):
        source = _js_source()
        assert source.count("window.location.href =") == 1
        assert not re.search(
            r'window\.location\.href\s*=\s*[^;]*"/pair/"', source)
        assert '"/pair/" + encodeURIComponent(cfg.type)' in source
        assert "openPair(pair.idx)" in source
        assert "openPair(rows[hit.row].idx, selectedFrame)" in source
