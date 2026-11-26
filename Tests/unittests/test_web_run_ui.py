# -*- coding: utf-8 -*-
"""W3 Web 发起计算 UI 测试：run/job 模板渲染（jinja2 直载，不依赖 app）+ 静态资源断言。"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WEB_DIR = REPO_ROOT / "DuIvyInteractions" / "web"
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"
RUN_JS = STATIC_DIR / "run.js"
FILEPICK_JS = STATIC_DIR / "filepick.js"
APP_CSS = STATIC_DIR / "app.css"
INDEX_HTML = TEMPLATE_DIR / "index.html"

INTERACTIONS = ["halogen_bond", "hydrogen_bond", "hydrophobic",
                "metal_coordination", "pi_cation", "pi_stacking",
                "salt_bridge", "water_bridge"]
FFS = ["amber", "gromos", "charmm", "opls"]
STRATEGIES = ["two_pass", "per_frame", "per_tuple"]
PAIRINGS = ["any", "between"]
CHECKED = ["hydrogen_bond", "pi_stacking"]


@pytest.fixture(scope="module")
def env():
    """独立 jinja2 环境：只加载模板目录，不依赖 FastAPI app。"""
    from jinja2 import Environment, FileSystemLoader, select_autoescape

    return Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)),
                       autoescape=select_autoescape(["html", "xml"]))


def _form(**overrides):
    form = {
        "topology": "",
        "trajectory": "",
        "output": "",
        "ff": "amber",
        "strategy": "two_pass",
        "interactions": list(CHECKED),
        "select": "",
        "pairing": "any",
        "mark_redundant": False,
    }
    form.update(overrides)
    return form


def _run_context(**overrides):
    context = {
        "title": "new calculation",
        "browse_mode": "picker",
        "root": "/tmp/root",
        "ffs": list(FFS),
        "strategies": list(STRATEGIES),
        "pairings": list(PAIRINGS),
        "interactions": list(INTERACTIONS),
        "form": _form(),
        "error": "",
        "recent_jobs": [],
    }
    context.update(overrides)
    return context


def _job(**overrides):
    job = {
        "id": "jobabc",
        "status": "running",
        "summary": "dii run -t /data/md.tpr",
        "params_preview": "dii run -t /data/md.tpr -f /data/md.xtc",
        "lines": ["line one", "line two"],
        "lines_total": 2,
        "returncode": None,
        "elapsed": 1.5,
        "output_dir": "/data/out",
        "h5_files": [],
        "error": None,
    }
    job.update(overrides)
    return job


def _job_context(job, **overrides):
    context = {
        "title": "job",
        "browse_mode": "picker",
        "root": "/tmp/root",
        "job": job,
        "recent_jobs": [],
    }
    context.update(overrides)
    return context


def _render(env, name, context):
    return env.get_template(name).render(**context)


def _element(html, tag, element_id):
    pattern = r"<%s[^>]*id=\"%s\"[^>]*>" % (tag, re.escape(element_id))
    match = re.search(pattern, html)
    assert match, f"missing <{tag} id=\"{element_id}\">"
    return match.group(0)


def _select_block(html, select_id):
    pattern = r"<select[^>]*id=\"%s\"[^>]*>(.*?)</select>" % re.escape(select_id)
    match = re.search(pattern, html, re.S)
    assert match, f"missing <select id=\"{select_id}\">"
    return match.group(1)


# ============================================================
# run.html
# ============================================================


class TestRunTemplate:
    def test_page_has_three_browse_buttons_with_mode(self, env):
        html = _render(env, "run.html", _run_context())
        topology = _element(html, "button", "browseTopology")
        trajectory = _element(html, "button", "browseTrajectory")
        output = _element(html, "button", "browseOutput")
        for button, target in ((topology, "topologyInput"),
                               (trajectory, "trajectoryInput"),
                               (output, "outputInput")):
            assert 'data-browse="picker"' in button, button
            assert 'data-target="%s"' % target in button, button
        assert 'data-kind="dir"' in output
        assert 'data-kind="dir"' not in topology
        assert html.count('data-browse="picker"') == 3

    def test_page_renders_server_browse_mode(self, env):
        html = _render(env, "run.html", _run_context(browse_mode="server"))
        assert html.count('data-browse="server"') == 3
        assert 'id="filepickModal"' in html
        assert "/static/filepick.js" in html
        assert 'id="filepickUseFolder"' in html

    def test_page_lists_force_fields(self, env):
        html = _render(env, "run.html", _run_context())
        options = re.findall(r'<option value="([^"]+)"', _select_block(html, "ffSelect"))
        assert options == FFS
        assert 'value="amber" selected' in html

    def test_page_lists_strategies_and_pairings(self, env):
        html = _render(env, "run.html", _run_context())
        strategies = re.findall(r'<option value="([^"]+)"',
                                _select_block(html, "strategySelect"))
        pairings = re.findall(r'<option value="([^"]+)"',
                              _select_block(html, "pairingSelect"))
        assert strategies == STRATEGIES
        assert pairings == PAIRINGS
        assert 'value="two_pass" selected' in html

    def test_page_lists_eight_interaction_checkboxes(self, env):
        html = _render(env, "run.html", _run_context())
        ids = re.findall(r'id="pt_([a-z_]+)"', html)
        assert ids == INTERACTIONS
        assert html.count('name="interactions"') == 8
        checked = re.findall(r'<input[^>]*value="([a-z_]+)"[^>]*checked[^>]*>', html)
        assert checked == CHECKED

    def test_page_has_required_field_ids_and_names(self, env):
        html = _render(env, "run.html", _run_context())
        fields = {
            "topologyInput": "topology",
            "trajectoryInput": "trajectory",
            "outputInput": "output",
            "ffSelect": "ff",
            "strategySelect": "strategy",
            "selectInput": "select",
            "pairingSelect": "pairing",
            "markRedundant": "mark_redundant",
        }
        for element_id, name in fields.items():
            tag = "select" if element_id.endswith("Select") else (
                "textarea" if element_id == "selectInput" else "input")
            element = _element(html, tag, element_id)
            assert 'name="%s"' % name in element, element
        for button_id in ("interactionsAll", "interactionsNone", "startRun"):
            _element(html, "button", button_id)
        assert 'id="commandPreview"' in html
        assert 'action="/run"' in html and 'method="post"' in html

    def test_page_error_renders_banner(self, env):
        html = _render(env, "run.html", _run_context(error="bad topology path"))
        banner = _element(html, "div", "runError")
        assert 'class="banner' in banner
        assert "bad topology path" in html
        empty = _render(env, "run.html", _run_context())
        assert 'id="runError"' not in empty

    def test_page_recent_jobs_rendered_only_when_present(self, env):
        recent = [{"id": "job-1", "status": "done", "summary": "run ok",
                   "elapsed": 12.3, "output": "/data/out"}]
        html = _render(env, "run.html", _run_context(recent_jobs=recent))
        assert 'href="/run/job-1"' in html
        assert "run ok" in html
        assert "12.3s" in html
        assert "/data/out" in html
        assert "job-state-done" in html
        empty = _render(env, "run.html", _run_context())
        assert "panel run-recent" not in empty

    def test_page_has_notes_and_default_hidden(self, env):
        html = _render(env, "run.html", _run_context())
        for note_id in ("perTupleNote", "markNote"):
            element = _element(html, "p", note_id)
            assert "hidden" in element
        assert "no longer maintained" in html
        assert "hydrophobic" in html


# ============================================================
# job.html
# ============================================================


class TestJobTemplate:
    def test_page_has_log_status_controls(self, env):
        html = _render(env, "job.html", _job_context(_job()))
        elements = [("div", "jobStatus"), ("span", "jobState"),
                    ("span", "jobElapsed"), ("button", "cancelRun"),
                    ("pre", "jobLog"), ("pre", "jobCommand"),
                    ("div", "jobResult"), ("p", "jobStatusError")]
        for tag, element_id in elements:
            _element(html, tag, element_id)
        assert 'data-job-id="jobabc"' in html
        assert "dii run -t /data/md.tpr" in html

    def test_page_renders_server_lines_with_total_cursor(self, env):
        job = _job(lines=["line one", "line two", "line three"],
                   lines_total=800)
        html = _render(env, "job.html", _job_context(job))
        log = _element(html, "pre", "jobLog")
        assert 'data-since="800"' in log
        assert "line one" in html and "line two" in html and "line three" in html

    def test_page_since_defaults_to_zero_without_total(self, env):
        job = _job()
        job.pop("lines_total")
        html = _render(env, "job.html", _job_context(job))
        assert 'data-since="0"' in _element(html, "pre", "jobLog")

    def test_job_template_uses_total_log_cursor(self):
        template = (TEMPLATE_DIR / "job.html").read_text(encoding="utf-8")
        assert 'data-since="{{ job.lines_total | default(0) }}"' in template
        assert "job.lines | length" not in template

    def test_page_done_branch_lists_h5_files(self, env):
        job = _job(status="done", returncode=0,
                   h5_files=["hydrogen_bond.h5", "pi_stacking.h5"])
        html = _render(env, "job.html", _job_context(job))
        assert "hydrogen_bond.h5" in html
        assert "pi_stacking.h5" in html
        _element(html, "button", "openResults")
        _element(html, "a", "downloadLog")
        assert "/download/run-log/jobabc" in html
        assert "disabled" in _element(html, "button", "cancelRun")

    def test_page_failed_branch_shows_error(self, env):
        job = _job(status="failed", returncode=1, error="Traceback: boom",
                   h5_files=[])
        html = _render(env, "job.html", _job_context(job))
        assert "Traceback: boom" in html
        assert "reproduced in a terminal" in html
        assert "disabled" in _element(html, "button", "cancelRun")

    def test_page_result_block_is_hidden_by_default(self, env):
        html = _render(env, "job.html", _job_context(_job()))
        result = _element(html, "div", "jobResult")
        assert "hidden" in result

    def test_page_recent_jobs_rendered_only_when_present(self, env):
        recent = [{"id": "job-2", "status": "failed", "summary": "bad ff",
                   "elapsed": 3.0, "output": "/data/out2"}]
        html = _render(env, "job.html", _job_context(_job(), recent_jobs=recent))
        assert 'href="/run/job-2"' in html
        assert "bad ff" in html
        empty = _render(env, "job.html", _job_context(_job()))
        assert "panel run-recent" not in empty


# ============================================================
# 静态资源
# ============================================================


class TestStaticFrontend:
    def test_run_js_contract(self):
        text = RUN_JS.read_text(encoding="utf-8")
        for token in ("commandPreview", "/api/run/", "since", "800",
                      "openResults", "/select", "cancelRun", "/api/pick",
                      "DiiFilePick", "POLL_INTERVAL_MS = 800"):
            assert token in text, token
        assert "innerHTML" not in text
        assert "[--select" not in text

    def test_run_js_builds_equivalent_command(self):
        text = RUN_JS.read_text(encoding="utf-8")
        for token in ("dii run", "--interactions", "--select", "--pairing",
                      "--mark-redundant", "<topology>", "<output>",
                      "per_tuple", "hydrophobic"):
            assert token in text, token
        assert 'replace(/"/g' in text, "select lines must escape double quotes"

    def test_filepick_js_supports_dir_mode_without_breaking_file_mode(self):
        text = FILEPICK_JS.read_text(encoding="utf-8")
        for token in ("KIND_DIR", "KIND_FILE", "filepickUseFolder",
                      "use this folder", "DATASET_ENDPOINT", "kind=dataset",
                      "onPick", 'ENDPOINT = "/api/list?kind=structure&path="',
                      "isStructureFile", "onPick = opts.onPick"):
            assert token in text, token
        assert "innerHTML" not in text
        assert "state.kind !== KIND_DIR" in text

    def test_css_has_run_and_job_rules(self):
        text = APP_CSS.read_text(encoding="utf-8")
        assert ".run-grid" in text
        assert ".run-input-row" in text
        assert ".run-command" in text
        assert ".job-status-row" in text
        assert ".job-log" in text
        assert ".job-state-failed" in text
        assert re.findall(r"^\.job-", text, re.M)
        assert re.findall(r"^\.run-", text, re.M)

    def test_css_braces_balanced(self):
        css = APP_CSS.read_text(encoding="utf-8")
        assert css.count("{") == css.count("}")

    def test_index_has_run_entry(self):
        html = INDEX_HTML.read_text(encoding="utf-8")
        assert 'href="/run"' in html
        assert "Run job" in html
        assert "New calculation" not in html
