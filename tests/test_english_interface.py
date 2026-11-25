"""English interface regressions; no GUI window, fitting, or sampling is started."""
import ast
import json
from pathlib import Path
import re
import sys
from threading import Event
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ce3d_xps_fitter_gui as gui


HAN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")


class Value:
    """A variable substitute that never initializes Tk."""
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value

    def set(self, value):
        self.value = value


@pytest.mark.parametrize("filename", ["ce3d_xps_fitter_gui.py", "ce3d_xps_fitter_core.py"])
def test_interface_source_string_literals_have_no_chinese(filename):
    tree = ast.parse((ROOT / filename).read_text(encoding="utf-8-sig"))
    violations = [(node.lineno, node.value) for node in ast.walk(tree)
                  if isinstance(node, ast.Constant) and isinstance(node.value, str)
                  and HAN.search(node.value)]
    assert not violations, f"Non-English interface strings in {filename}: {violations}"


def test_file_dialogs_and_vamas_prompt_use_english(monkeypatch):
    app = object.__new__(gui.Ce3dXpsFitterApp)
    app.root = object()
    app._add_files = Mock()
    app.output_dir = Value("")
    choose_files = Mock(return_value=())
    choose_folder = Mock(return_value="")
    show_info = Mock()
    monkeypatch.setattr(gui.filedialog, "askopenfilenames", choose_files)
    monkeypatch.setattr(gui.filedialog, "askdirectory", choose_folder)
    monkeypatch.setattr(gui.messagebox, "showinfo", show_info)
    app._browse_files()
    app._browse_output()
    assert choose_files.call_args.kwargs["title"] == "Select one or more Ce 3d data files"
    assert choose_folder.call_args.kwargs["title"] == "Select output folder"
    assert [label for label, _ in choose_files.call_args.kwargs["filetypes"]] == [
        "XPS data (automatic format detection)", "Excel workbooks", "VAMAS spectra", "Text data", "All files"]
    app.active = False
    app._inspecting_vamas = False
    app.input_files = []
    app.file_list = SimpleNamespace(curselection=lambda: ())
    app._choose_vamas()
    assert show_info.call_args.args == (
        "Select a VAMAS file", "Add a VAMAS file, then select the file to inspect in the input list.")


@pytest.mark.parametrize("method,args,message", [
    ("_integer", ("invalid", "Walkers", 22), "Walkers must be an integer."),
    ("_integer", ("21", "Walkers", 22), "Walkers must be at least 22."),
    ("_number", ("invalid", "KE offset"), "KE offset must be a number."),
    ("_number", ("nan", "KE offset"), "KE offset must be a finite number."),
    ("_number", ("0", "Satellite gamma", True), "Satellite gamma must be a finite positive number."),
])
def test_numeric_error_messages_are_complete_english_sentences(method, args, message):
    with pytest.raises(ValueError) as error:
        getattr(gui.Ce3dXpsFitterApp, method)(*args)
    assert str(error.value) == message


@pytest.fixture
def argument_app(tmp_path):
    app = object.__new__(gui.Ce3dXpsFitterApp)
    spectrum = tmp_path / "input.csv"
    spectrum.write_text("Energy,Intensity\n870,1\n925,1\n", encoding="utf-8")
    app.input_files = [spectrum]
    for name, value in dict(output_dir=str(tmp_path / "results"), model="LA", satellites=False,
                            mode="quick", walkers="80", quick_steps="8000", burn="4000", seed="42",
                            background_mode="fixed_shirley", satellite_gamma="0.8", vamas_block="",
                            vamas_ordinate="", kinetic_offset="").items():
        setattr(app, name, Value(value))
    return app


@pytest.mark.parametrize("field,value,message", [
    ("output_dir", "", "Select an output folder."),
    ("walkers", "23", "Walkers must be an even number."),
    ("quick_steps", "3000", "Quick-fit steps must be a multiple of 2000."),
    ("quick_steps", "4000", "Burn-in steps must be less than the total steps for a quick fit."),
    ("burn", "3000", "Burn-in steps must be a multiple of 2000."),
    ("seed", "4294967296", "Random seed must be between 0 and 4294967295."),
    ("vamas_block", "-1", "VAMAS block index must be at least 0."),
    ("vamas_ordinate", "invalid", "Signal channel index must be an integer."),
    ("kinetic_offset", "inf", "KE offset must be a finite number."),
])
def test_advanced_validation_messages_are_english(argument_app, field, value, message):
    getattr(argument_app, field).set(value)
    with pytest.raises(ValueError) as error:
        argument_app._snapshot_arguments()
    assert str(error.value) == message


def test_missing_input_validation_uses_english(argument_app):
    argument_app.input_files = []
    with pytest.raises(ValueError, match="^Select at least one input file[.]$"):
        argument_app._snapshot_arguments()
    argument_app.input_files = [Path("missing_spectrum.csv")]
    with pytest.raises(ValueError, match="^Some input files no longer exist; select them again[.]$"):
        argument_app._snapshot_arguments()


@pytest.mark.parametrize("mode,phrase", [
    ("quick", "Completion does not establish convergence."),
    ("convergence", "no preset step limit. Stop and save is available."),
])
def test_mode_status_is_english_and_preserves_convergence_meaning(mode, phrase):
    app = object.__new__(gui.Ce3dXpsFitterApp)
    app.mode = Value(mode)
    app.active = False
    app.quick_steps_entry = Mock()
    app.step_status = Value("")
    app._sync_mode()
    assert phrase in app.step_status.get()
    assert not HAN.search(app.step_status.get())


def test_progress_message_uses_english_without_sampling():
    app = object.__new__(gui.Ce3dXpsFitterApp)
    app.step_status = Value("")
    app.status = Value("")
    app.progress_bar = Mock()
    app.stop_event = Event()
    app._update_progress(dict(file_index=1, file_count=2, steps=2000, total_steps=8000,
                              elapsed_seconds=90, Ce3_quantiles_percent=[10, 11, 12]))
    assert app.step_status.get() == "File 1/2 · 2000 steps sampled · Elapsed: 1.5 min · Ce³⁺ 11.00% [10.00, 12.00]"
    assert app.status.get() == "Sampling..."
    app._update_progress(dict(event="complete"))
    assert app.status.get() == "Preparing results..."


@pytest.mark.parametrize("saved_status,display_status,phrase", [
    ("cancelled", "cancelled", "Convergence was not established."),
    ("failed", "failed", "Run finished with failed files"),
    ("budget_exhausted", "incomplete", "before the engineering convergence criteria were met"),
    ("converged", "converged", "Engineering convergence criteria met."),
    ("quick_complete", "quick_complete", "review diagnostics to assess reliability"),
    ("unknown", "unknown", "review the actual status"),
])
def test_saved_result_status_messages_use_english(tmp_path, saved_status, display_status, phrase):
    (tmp_path / "fit_result.json").write_text(json.dumps(dict(status=saved_status)), encoding="utf-8")
    status, message = gui.inspect_result(tmp_path, 1)
    assert status == display_status
    assert phrase in message
    assert not HAN.search(message)


@pytest.mark.parametrize("batch_status,display_status,phrase", [
    ("cancelled", "cancelled", "Convergence was not established."),
    ("completed_with_errors", "failed", "Batch finished with failed files"),
])
def test_batch_status_messages_use_english(tmp_path, batch_status, display_status, phrase):
    (tmp_path / "batch_summary.json").write_text(json.dumps(dict(status=batch_status)), encoding="utf-8")
    status, message = gui.inspect_result(tmp_path, 2)
    assert status == display_status
    assert phrase in message
    assert not HAN.search(message)
