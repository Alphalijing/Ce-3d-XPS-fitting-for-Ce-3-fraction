"""GUI result routing tests; no window, fitting or sampling is started."""
from pathlib import Path
import json
import queue
import sys
import tempfile
from threading import Event, Lock
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ce3d_xps_fitter_gui as gui


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


class ResultRoutingTests(unittest.TestCase):
    def test_organized_single_result_is_not_double_counted_with_legacy_copy(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            save(root / "fit_result.json", {"status": "failed"})
            save(root / "details" / "fit_result.json", {"status": "quick_complete"})
            self.assertEqual(len(gui._saved_result_records(root)), 1)
            self.assertEqual(gui.inspect_result(root, 1)[0], "quick_complete")

    def test_batch_details_summary_is_preferred_and_child_results_are_deduplicated(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            save(root / "batch_summary.json", {"status": "completed_with_errors"})
            save(root / "details" / "batch_summary.json", {
                "status": "batch_complete", "results": [{"status": "quick_complete"}, {"status": "converged"}]})
            self.assertEqual(gui.inspect_result(root, 2)[0], "quick_complete")
            (root / "batch_summary.json").unlink()
            (root / "details" / "batch_summary.json").unlink()
            for sample in ("one", "two"):
                save(root / sample / "fit_result.json", {"status": "failed"})
                save(root / sample / "details" / "fit_result.json", {"status": "converged"})
            self.assertEqual(len(gui._saved_result_records(root)), 2)
            self.assertEqual(gui.inspect_result(root, 2)[0], "converged")

    def test_open_results_chooses_report_then_legacy_folder(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            app = object.__new__(gui.Ce3dXpsFitterApp)
            app.last_output = root
            app.root = object()
            with patch.object(gui.os, "startfile") as opened:
                app._open_results()
                opened.assert_called_once_with(str(root))
            (root / "Results.html").write_text("<html></html>", encoding="utf-8")
            with patch.object(gui.os, "startfile") as opened:
                app._open_results()
                opened.assert_called_once_with(str(root / "Results.html"))

    def test_final_worker_queues_mean_sd_summary_before_report_path(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            run = root / "sample"
            save(run / "details" / "fit_result.json", {
                "status": "quick_complete", "posterior_available": True,
                "Ce3_fraction_mean": .123, "Ce3_fraction_std": .004,
                "Ce3_percent": {"q025": 11.5, "median": 12.29, "q975": 13.1}})
            app = object.__new__(gui.Ce3dXpsFitterApp)
            app.events = queue.Queue()
            app._log_lock = Lock()
            kwargs = {"model": "LA", "satellites": False, "mode": "quick", "outdir": root}
            with patch("ce3d_xps_fitter_core.run_fit", return_value=run) as fit:
                app._fit_worker([root / "input.csv"], kwargs, root / "gui.log", Event())
            fit.assert_called_once()
            events = list(app.events.queue)
            summaries = [text for kind, text in events if kind == "result_summary"]
            self.assertEqual(len(summaries), 1)
            self.assertIn("12.30", summaries[0])
            self.assertIn("0.40", summaries[0])
            self.assertIn("±", summaries[0])
            self.assertIn("12.29", summaries[0])
            lines = [str(text) for kind, text in events if kind == "log"]
            mean_line = next(i for i, line in enumerate(lines) if "±" in line and "12.30" in line)
            report_line = next(i for i, line in enumerate(lines) if "Results.html" in line)
            self.assertLess(mean_line, report_line)
            self.assertFalse(any(kind == "error" for kind, _ in events))

    def test_summary_queue_changes_only_existing_summary_label(self):
        app = object.__new__(gui.Ce3dXpsFitterApp)
        app.events = queue.Queue()
        app.events.put(("result_summary", "Ce3+ = 12.30 ± 0.40%"))
        app.step_status = Mock()
        app.step_status_label = Mock()
        app.close_when_finished = False
        app.root = Mock()
        app._drain_events()
        app.step_status.set.assert_called_once_with("Ce3+ = 12.30 ± 0.40%")
        font = app.step_status_label.configure.call_args.kwargs["font"]
        self.assertEqual(font[1:], (11, "bold"))
        app.root.after.assert_called_once()


if __name__ == "__main__":
    unittest.main()
