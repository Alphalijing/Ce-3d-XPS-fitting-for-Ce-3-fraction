"""Exercise real output writers without running an optimizer or MCMC.

All spectra and chain coordinates in this check are deterministic fixtures.
They are not measured data, posterior samples, or a convergence validation.
The same check runs in the source tree and in the frozen Windows executable.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import struct
from tempfile import TemporaryDirectory
import xml.etree.ElementTree as ET

import emcee
import numpy as np

from ce3d_io import SpectrumInput
from ce3d_xps_fitter_core import ChainStore, _export, _save_background, write_json
from ce3d_xps_model import Ce3dModel
from ce3d_results import organize_outputs


def _csv_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError(f"Empty exported CSV: {path.name}")
    return rows


def _png_dimensions(path: Path) -> list[int]:
    header = path.read_bytes()[:24]
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise AssertionError(f"Invalid exported PNG: {path.name}")
    dimensions = list(struct.unpack(">II", header[16:24]))
    if min(dimensions) <= 0:
        raise AssertionError(f"Invalid PNG dimensions: {path.name}")
    return dimensions


def _check_case(folder: Path, lineshape: str, satellites: bool) -> dict:
    folder.mkdir()
    energy = np.linspace(870.0, 925.0, 111)
    seed_signal = 100.0 + 300.0 * np.exp(-0.5 * ((energy - 900.0) / 5.0) ** 2)
    model = Ce3dModel(energy, seed_signal, lineshape=lineshape, satellites=satellites)
    theta = model.initial()
    theta[2:4] = [0.2, 10000.0]
    background = 20.0 + 0.1 * (energy - energy[0])
    raw = background + model.predict(theta)[0] + 0.01 * np.sin(energy)
    spectrum = SpectrumInput(
        energy=energy, raw=raw, background=background,
        metadata={"fixture": True, "measured_data": False, "description": "Deterministic output-writer check"},
        notes=["Synthetic fixture; no MCMC sampling"], source=folder / "fixture.csv",
    )
    config = {"background": "fixed_shirley", "model": lineshape, "satellites": satellites,
              "fixture": True, "sampling_performed": False}
    write_json(folder / "experiment_metadata.json", spectrum.metadata)
    write_json(folder / "run_configuration.json", config)
    _save_background(folder, spectrum, label="Deterministic background fixture")

    # Small, finite coordinates exercise the real on-disk chain/export contract.
    # They are generated directly, not obtained by sampling a posterior.
    rng = np.random.default_rng(147)
    chain = theta + rng.normal(size=(8, 4, model.ndim)) * model.scales * 0.0005
    log_prob = np.stack([model.log_prob(block) for block in chain])
    if not np.isfinite(log_prob).all():
        raise AssertionError("Output-check fixture must have finite log probabilities")
    store = ChainStore(folder)
    random_state = np.random.RandomState(147).get_state()
    for start, stop in ((0, 4), (4, 8)):
        state = emcee.State(chain[stop - 1], log_prob=log_prob[stop - 1], random_state=random_state)
        store.append(chain[start:stop], log_prob[start:stop], state)
    saved_chain, saved_log_prob = store.read()
    np.testing.assert_array_equal(saved_chain, chain)
    np.testing.assert_array_equal(saved_log_prob, log_prob)
    with np.load(folder / "chain" / "last_state.npz", allow_pickle=False) as saved_state:
        np.testing.assert_array_equal(saved_state["coords"], chain[-1])
        np.testing.assert_array_equal(saved_state["log_prob"], log_prob[-1])
    manifest = json.loads((folder / "chain" / "manifest.json").read_text(encoding="utf-8"))
    if manifest["steps"] != 8 or len(manifest["blocks"]) != 2:
        raise AssertionError("Export-check chain manifest does not describe both saved blocks")
    state_json = json.loads((folder / "chain" / "random_state.json").read_text(encoding="utf-8"))
    if state_json["algorithm"] != random_state[0] or state_json["keys"] != random_state[1].tolist():
        raise AssertionError("Saved random state did not round-trip")

    burn = 2
    result = _export(folder, spectrum, model, store, burn, "export_check_fixture", config,
                     {"evaluated": False, "passed": False, "fixture": True}, 0.0)
    saved_result = json.loads((folder / "fit_result.json").read_text(encoding="utf-8"))
    if saved_result != result or result["status"] != "export_check_fixture" or not result["posterior_available"]:
        raise AssertionError("Exported result JSON did not round-trip")
    flat = chain[burn:].reshape(-1, model.ndim)
    expected_quantiles = 100.0 * np.quantile(flat[:, 2], [0.025, 0.5, 0.975])
    np.testing.assert_allclose([result["Ce3_percent"][key] for key in ("q025", "median", "q975")],
                               expected_quantiles, rtol=1e-12, atol=1e-12)
    parameters = np.median(flat, axis=0)

    samples = _csv_rows(folder / "posterior_samples_ce3d.csv")
    parameter_rows = _csv_rows(folder / "posterior_parameters.csv")
    peak_rows = _csv_rows(folder / "peak_table.csv")
    curve_rows = _csv_rows(folder / "fit_curves.csv")
    prepared = _csv_rows(folder / "prepared_spectrum.csv")
    if len(samples) != len(flat) or len(parameter_rows) != model.ndim or len(peak_rows) != model.n_components:
        raise AssertionError("Exported table dimensions do not match the fixture")
    np.testing.assert_allclose([[float(row[name]) for name in model.param_names] for row in samples], flat,
                               rtol=1e-12, atol=1e-12)
    if len(curve_rows) != len(energy) or len(prepared) != len(energy):
        raise AssertionError("Exported spectrum length changed")
    component_curves = np.array([[float(row[name]) for name in model.component_names] for row in curve_rows])
    net_fit = np.array([float(row["net_fit"]) for row in curve_rows])
    total_fit = np.array([float(row["total_fit"]) for row in curve_rows])
    residual = np.array([float(row["residual"]) for row in curve_rows])
    np.testing.assert_allclose(component_curves.sum(axis=1), net_fit, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(total_fit, background + net_fit, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(residual, raw - total_fit, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose([float(row["background"]) for row in prepared], background, rtol=1e-12)
    areas = np.trapz(component_curves, energy, axis=0)
    np.testing.assert_allclose([float(row["area"]) for row in peak_rows], areas, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(areas.sum(), parameters[3], rtol=1e-12)
    np.testing.assert_allclose(areas[:4].sum() / areas.sum(), parameters[2], rtol=1e-12)

    image_sizes = {name: _png_dimensions(folder / name) for name in ("background.png", "fit_ce3d.png")}
    svg = ET.parse(folder / "fit_ce3d.svg").getroot()
    if svg.tag != "{http://www.w3.org/2000/svg}svg" or not svg.findall(".//{http://www.w3.org/2000/svg}path"):
        raise AssertionError("Exported SVG is not a nonempty SVG document")
    report = folder / "Results.html"
    if not report.is_file() or report.stat().st_size < 1000:
        raise AssertionError("Primary results report was not generated")
    organize_outputs(folder)
    if {p.name for p in folder.iterdir()} != {"Results.html", "figures", "tables", "details"}:
        raise AssertionError("Output directory is not organized around the primary report")
    files = sorted(path.relative_to(folder).as_posix() for path in folder.rglob("*") if path.is_file())
    return {"model": lineshape, "satellites": satellites, "parameter_count": model.ndim,
            "component_count": model.n_components, "chain_steps": store.steps, "chain_blocks": 2,
            "exported_sample_rows": len(samples), "files_checked": files,
            "png_dimensions": image_sizes, "svg_xml_valid": True,
            "area_conservation_passed": True, "result_quantiles_match_fixture": True,
            "primary_report_generated": True, "organized_layout": True}


def run_export_self_check() -> dict:
    """Write and read every supported model's outputs in a temporary folder.

    Failure raises an exception so the GUI/frozen self-check reports a failure.
    Temporary fixture files are removed on success and on failure.
    """
    with TemporaryDirectory(prefix="ce3d_export_check_") as temporary:
        root = Path(temporary)
        cases = [_check_case(root / f"{shape}_sat_{int(satellites)}", shape, satellites)
                 for shape in ("LA", "Voigt") for satellites in (False, True)]
    if root.exists():
        raise AssertionError("Temporary output-check files were not removed")
    return {"status": "passed", "fixture": "Deterministic synthetic export fixture; not MCMC samples or measured data",
            "sampling_performed": False, "formats": ["HTML", "PNG", "SVG", "CSV", "JSON", "NPZ"],
            "temporary_files_removed": True, "cases": cases}
