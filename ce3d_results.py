"""A result-first, offline report; full scientific records remain available."""
from __future__ import annotations

import base64
import csv
import html
import json
import math
from pathlib import Path
import re
import shutil
from urllib.parse import quote


FIGURES = {"background.png", "fit_ce3d.png", "fit_ce3d.svg", "posterior_Ce3.png",
           "posterior_Ce3.svg", "posterior_parameters.png", "posterior_parameters.svg"}
TABLES = {"peak_table.csv", "fit_curves.csv", "posterior_parameters.csv"}
DETAILS = {"experiment_metadata.json", "run_configuration.json", "fit_result.json",
           "posterior_samples_ce3d.csv", "prepared_spectrum.csv", "initialization.json",
           "frozen_proposal.json", "progress.json", "error.log", "initial_positions.npy",
           "export_recovery.json", "batch_summary.json"}


def _number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def _e(value):
    return html.escape(str(value), quote=True)


def _format(value, places=3):
    n = _number(value)
    return "—" if n is None else f"{n:.{places}f}"


def _group(name):
    if name in FIGURES:
        return "figures"
    if name in TABLES:
        return "tables"
    return "details"


def _find(out, name):
    root = Path(out)
    for path in (root / name, root / _group(name) / name):
        if path.is_file():
            return path
    return None


def _json(out, name):
    path = _find(out, name)
    if path is None:
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return {}


def _rows(out, name):
    path = _find(out, name)
    if path is None:
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _link(out, name, label):
    if _find(out, name) is None:
        return ""
    target = quote(_group(name) + "/" + name, safe="/")
    return f'<a href="{target}">{_e(label)}</a>'


def _image(out, name, alt):
    path = _find(out, name)
    if path is None:
        return '<p class="muted">Figure unavailable.</p>'
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'<img alt="{_e(alt)}" src="data:image/png;base64,{data}">'


def result_headline(result):
    """Mean ± posterior SD, both on a percentage scale (SD in points)."""
    mean = _number(result.get("Ce3_fraction_mean"))
    sd = _number(result.get("Ce3_fraction_std"))
    if result.get("posterior_available") and mean is not None and sd is not None and sd >= 0:
        return f"Ce³⁺ = {100 * mean:.2f} ± {100 * sd:.2f}%"
    median = _number(result.get("Ce3_percent", {}).get("median"))
    if result.get("posterior_available") and median is not None:
        return f"Ce³⁺ = {median:.2f}% (posterior median; SD unavailable)"
    return "Ce³⁺ result unavailable"


def _status(result):
    return {"quick_complete": "Quick fit complete", "converged": "Engineering convergence criteria met",
            "cancelled": "Stopped; available samples only", "budget_exhausted": "Step budget reached; convergence criteria not met",
            "failed": "Run failed", "export_check_fixture": "Export-test fixture"}.get(
                result.get("status"), str(result.get("status", "Status unavailable")))


def _conditions(metadata):
    """Only the selected Ce 3d scan; never substitute survey settings."""
    values = metadata.get("selected_ce3d_conditions")
    if not isinstance(values, dict):
        values = metadata.get("experiment_conditions", {})
    if not isinstance(values, dict):
        values = {}
    labels = [("source_gun_type", "X-ray source", ""), ("photon_energy_eV", "Photon energy", "eV"),
              ("pass_energy_eV", "Pass energy", "eV"), ("energy_step_eV", "Energy step", "eV"),
              ("total_acquisition_time_s", "Acquisition time", "s"), ("scan_count", "Scans", "")]
    pairs = []
    for key, label, unit in labels:
        value = values.get(key)
        if value is None or str(value).strip() == "":
            continue
        number = _number(value)
        value = f"{number:g}" if number is not None else str(value)
        pairs.append((label, value + (" " + unit if unit else "")))
    return pairs


def _posterior_plots(out, result, flat, names):
    import numpy as np
    import matplotlib.pyplot as plt

    samples = np.asarray(flat, dtype=float)
    if samples.ndim != 2 or samples.shape[0] == 0 or samples.shape[1] != len(names):
        raise ValueError("Posterior plot samples must be a nonempty matrix matching parameter names.")
    if not np.isfinite(samples).all():
        raise ValueError("Posterior plot samples must be finite.")
    if "x" not in names:
        raise ValueError("Posterior samples have no Ce(III) fraction parameter x.")
    folder = Path(out) / "figures"
    folder.mkdir(exist_ok=True)
    values = 100 * samples[:, names.index("x")]
    q = result.get("Ce3_percent", {})
    fig, ax = plt.subplots(figsize=(7.5, 2.8), layout="constrained")
    try:
        ax.hist(values, bins=55, density=True, color="#356c99", alpha=.82)
        if _number(q.get("q025")) is not None and _number(q.get("q975")) is not None:
            ax.axvspan(float(q["q025"]), float(q["q975"]), color="#d9e6f0", alpha=.5, zorder=0,
                       label="95% credible interval")
        mean = _number(result.get("Ce3_fraction_mean"))
        if mean is not None:
            ax.axvline(100 * mean, color="#cc6539", lw=1.4, label="Posterior mean")
        ax.set_xlabel("Ce(III) fraction (%)"); ax.set_ylabel("Posterior density")
        ax.legend(frameon=False, fontsize=8)
        for ext in ("png", "svg"):
            fig.savefig(folder / ("posterior_Ce3." + ext), dpi=160)
    finally:
        plt.close(fig)
    cols = 3
    rows = math.ceil(len(names) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(10, 2.05 * rows), layout="constrained", squeeze=False)
    labels = {"dE": "Energy shift (eV)", "sigma_ce": "Gaussian sigma (eV)", "x": "Ce(III) fraction (%)",
              "A_ce_total": "Total observed-window area", "log_a_ce": "log(a): noise variance",
              "log_b_ce": "log(b): signal-dependent noise"}
    try:
        for index, ax in enumerate(axes.flat):
            if index >= len(names):
                ax.set_visible(False)
                continue
            name = names[index]; column = samples[:, index] * (100 if name == "x" else 1)
            ax.hist(column, bins=45, density=True, color="#356c99", alpha=.8)
            title = labels.get(name, name)
            if re.fullmatch(r"g[34]_[123]", name):
                is_la = result.get("model", {}).get("lineshape") == "LA" and name == "g4_1"
                title = name + (": LA kernel scale (eV)" if is_la else ": Lorentzian HWHM (eV)")
            ax.set_title(title, fontsize=8.5)
            ax.tick_params(labelsize=7); ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
            ax.yaxis.offsetText.set_fontsize(7)
        for ext in ("png", "svg"):
            fig.savefig(folder / ("posterior_parameters." + ext), dpi=145)
    finally:
        plt.close(fig)


STYLE = """
:root{font-family:Segoe UI,Arial,sans-serif;color:#20313f;background:#f2f5f7;line-height:1.5}
body{max-width:1080px;margin:30px auto;padding:0 22px}header,section{background:white;border:1px solid #dde4e9;border-radius:12px;padding:24px 28px;margin:0 0 18px}
h1{font-size:38px;line-height:1.2;color:#124969;margin:4px 0 8px;font-weight:700}h2{font-size:21px;margin:0 0 16px}h3{font-size:16px}p{margin:8px 0}.eyebrow{font-size:12px;text-transform:uppercase;letter-spacing:.08em;color:#65737d}.muted{color:#65737d;font-size:13px}.sample{font-size:17px;overflow-wrap:anywhere}.badge{display:inline-block;background:#edf3f7;border-radius:5px;padding:3px 8px;font-size:12px;margin:3px 0}img{width:100%;height:auto;display:block;margin:10px auto}a{color:#17638d;text-decoration:none}a:hover{text-decoration:underline}.links{display:flex;gap:18px;flex-wrap:wrap;margin-top:10px;font-size:13px}dl{display:grid;grid-template-columns:160px 1fr 160px 1fr;gap:8px 16px;margin:0}dt{color:#65737d}dd{margin:0;overflow-wrap:anywhere}.table-wrap{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:12px;font-variant-numeric:tabular-nums}th{background:#eef3f6;font-weight:600;text-align:left}th,td{padding:9px 10px;border-bottom:1px solid #e2e8ec;white-space:nowrap}tbody tr:nth-child(even){background:#fafcfd}summary{cursor:pointer;font-weight:600;padding:8px 0}details{margin-top:16px}footer{padding:5px 20px 28px}.status-note{font-size:13px;color:#73551a}.empty{padding:20px;background:#f5f7f9;border-radius:6px}
@media(max-width:700px){body{margin:12px auto;padding:0 10px}header,section{padding:18px}h1{font-size:29px}dl{grid-template-columns:140px 1fr}}
@media print{body{max-width:none;margin:0;background:white;padding:0}header,section{break-inside:avoid;border:0;padding:12px 0}details:not([open]){display:none}.links{display:none}h1{font-size:30px}}
"""


def _write_page(path, title, body):
    text = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{_e(title)}</title><style>{STYLE}</style></head><body>{body}</body></html>')
    temporary = path.with_suffix(".tmp.html")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)
    return path


def write_result_report(out, result, flat=None, param_names=None):
    """Create Results.html; all its images are embedded for offline sharing."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    config = result.get("configuration", {}) or _json(out, "run_configuration.json")
    metadata = _json(out, "experiment_metadata.json")
    source = str(config.get("source") or metadata.get("source", {}).get("path") or out.name)
    sample = re.split(r"[\\/]", source)[-1]
    q = result.get("Ce3_percent", {})
    have_mean = _number(result.get("Ce3_fraction_mean")) is not None
    have_sd = _number(result.get("Ce3_fraction_std")) is not None
    available = bool(result.get("posterior_available"))
    subtitle = "Posterior mean ± 1 posterior SD (percentage points)" if available and have_mean and have_sd else "No complete posterior summary available"
    median_line = (f'Posterior median: {_format(q.get("median"), 2)}%; '
                   f'95% credible interval: {_format(q.get("q025"), 2)}–{_format(q.get("q975"), 2)}%.') if available else ""
    model = config.get("model") or result.get("model", {}).get("lineshape", "Model unavailable")
    background = "Fixed Shirley background" if config.get("background") == "fixed_shirley" else "Provided background"
    satellites = "Satellites included" if config.get("satellites") else "Satellites excluded"
    state = _status(result)
    note = ""
    if result.get("status") == "quick_complete":
        note = '<p class="status-note">Requested quick-fit budget completed; this status does not establish convergence.</p>'
    elif result.get("status") in {"cancelled", "budget_exhausted", "failed"}:
        note = '<p class="status-note">Use the saved run status when interpreting these available results.</p>'
    body = (f'<header><div class="eyebrow">Ce 3d XPS result</div><h1>{_e(result_headline(result))}</h1>'
            f'<p class="muted">{_e(subtitle)}</p><p>{_e(median_line)}</p><p class="sample">{_e(sample)}</p>'
            f'<p class="muted">{_e(model)} · {_e(background)} · {_e(satellites)}</p>'
            f'<span class="badge">{_e(state)}</span>{note}</header>')
    if not available:
        error = result.get("error")
        body += '<section><h2>Available output</h2><p>Sampling did not provide retained posterior samples.</p>'
        if error:
            body += f'<p>{_e(error)}</p>'
        body += '<div class="links">' + _link(out, "experiment_metadata.json", "Experimental conditions JSON") + _link(out, "fit_result.json", "Saved run status") + '</div></section>'
        return _write_page(out / "Results.html", "Ce 3d XPS result — " + sample, body)
    body += '<section id="fit"><h2>Peak fit</h2>' + _image(out, "fit_ce3d.png", "Measured Ce 3d spectrum, fitted components and residual")
    body += '<div class="links">' + _link(out, "fit_ce3d.svg", "Vector figure (SVG)") + _link(out, "fit_curves.csv", "Fit curves (CSV)") + '</div></section>'
    conditions = _conditions(metadata)
    body += '<section id="conditions"><h2>Experimental conditions</h2>'
    if conditions:
        body += '<dl>' + ''.join(f'<dt>{_e(k)}</dt><dd>{_e(v)}</dd>' for k, v in conditions) + '</dl>'
    else:
        body += '<p class="muted">Acquisition conditions were not recorded in the selected input scan.</p>'
    body += '<div class="links">' + _link(out, "experiment_metadata.json", "Full experimental conditions (JSON)") + '</div></section>'
    peaks = _rows(out, "peak_table.csv")
    body += '<section id="peaks"><h2>Component fit results</h2>'
    if peaks:
        body += '<div class="table-wrap"><table><thead><tr>' + ''.join('<th>' + h + '</th>' for h in
                 ['Component', 'Model', 'Position / eV', 'Gaussian σ / eV', 'Lorentzian HWHM / eV', 'LA scale γ / eV', 'Area', 'Area / %']) + '</tr></thead><tbody>'
        for peak in peaks:
            fraction = _number(peak.get("fraction_of_total"))
            cells = [_e(peak.get("component", "")), _e(peak.get("lineshape", "")), _format(peak.get("center_eV")),
                     _format(peak.get("sigma_eV")), _format(peak.get("Lorentzian_HWHM_eV")), _format(peak.get("LA_kernel_gamma_eV")),
                     _format(peak.get("area"), 2), _format(None if fraction is None else 100 * fraction, 2)]
            body += '<tr>' + ''.join('<td>' + c + '</td>' for c in cells) + '</tr>'
        body += '</tbody></table></div><p class="muted">At the posterior median parameter vector; areas refer to the fitted energy window. Gaussian σ is shared. LA kernel scale is not a Lorentzian HWHM.</p>'
    else:
        body += '<p class="muted">Component table unavailable.</p>'
    body += '<div class="links">' + _link(out, "peak_table.csv", "Component table (CSV)") + '</div></section>'
    names = list(param_names) if param_names is not None else None
    if flat is None:
        path = _find(out, "posterior_samples_ce3d.csv")
        if path is not None:
            import pandas as pd
            frame = pd.read_csv(path)
            flat = frame.to_numpy(dtype=float)
            names = list(frame.columns)
    if flat is not None:
        if names is None:
            names = result.get("model", {}).get("parameter_names")
        if not names:
            raise ValueError("Parameter names are required to plot the posterior samples.")
        _posterior_plots(out, result, flat, list(names))
    body += '<section id="posterior"><h2>Posterior distributions</h2>'
    body += _image(out, "posterior_Ce3.png", "Ce(III) posterior distribution and credible interval")
    posterior_rows = _rows(out, "posterior_parameters.csv")
    body += '<details><summary>All fitted parameters: distributions and numerical summaries</summary>'
    if _find(out, "posterior_parameters.png") is not None:
        body += _image(out, "posterior_parameters.png", "Marginal posterior distributions of fitted parameters")
    if posterior_rows:
        body += '<div class="table-wrap"><table><thead><tr><th>Parameter</th><th>Median</th><th>Mean</th><th>SD</th><th>2.5%</th><th>97.5%</th></tr></thead><tbody>'
        for row in posterior_rows:
            body += '<tr><td>' + _e(row.get("parameter", "")) + '</td>' + ''.join('<td>' + _format(row.get(key), 5) + '</td>' for key in ['median', 'mean', 'std', 'q025', 'q975']) + '</tr>'
        body += '</tbody></table></div>'
    body += '<div class="links">' + _link(out, "posterior_parameters.csv", "Posterior summary (CSV)") + _link(out, "posterior_samples_ce3d.csv", "Full posterior samples (CSV)") + '</div></details></section>'
    body += '<section><details><summary>Full records and diagnostics</summary><div class="links">'
    for name, label in [("fit_result.json", "Result and diagnostics JSON"), ("run_configuration.json", "Run configuration JSON"),
                        ("background.png", "Background figure"), ("prepared_spectrum.csv", "Prepared spectrum CSV")]:
        body += _link(out, name, label)
    body += '</div><p class="muted">Raw chain blocks, checkpoints and other technical records are retained in the details folder.</p></details></section>'
    body += '<footer class="muted">Uncertainty is conditional on the chosen background, peak model, priors and noise likelihood. It does not include all model or experimental uncertainty.</footer>'
    return _write_page(out / "Results.html", "Ce 3d XPS result — " + sample, body)


def organize_outputs(out):
    """Move only known records under this run; never remove unknown user files."""
    out = Path(out).resolve()
    if not out.is_dir():
        raise NotADirectoryError(out)
    known = FIGURES | TABLES | DETAILS
    known |= {p.name for p in out.glob("checkpoint_*.json") if re.fullmatch(r"checkpoint_\d+\.json", p.name)}
    entries = [(name, _group(name)) for name in sorted(known)] + [("chain", "details")]
    for name, group in entries:
        source = out / name
        if not source.exists():
            continue
        if source.is_symlink() or not source.resolve().is_relative_to(out):
            raise ValueError(f"Refusing to move a linked output: {name}")
        if name == "chain":
            if not source.is_dir():
                raise ValueError("The saved chain must be a directory.")
        elif not source.is_file():
            raise ValueError(f"Expected an output file: {name}")
        folder = out / group
        if folder.is_symlink():
            raise ValueError(f"Refusing a linked result folder: {group}")
        folder.mkdir(exist_ok=True)
        destination = folder / name
        if destination.exists():
            raise FileExistsError(f"Output exists in both old and organized locations: {name}")
        shutil.move(str(source), str(destination))


def write_batch_report(root, results, status):
    """One batch landing page with per-sample summaries and report links."""
    root = Path(root).resolve(); root.mkdir(parents=True, exist_ok=True)
    body = '<header><h1>Ce³⁺ results</h1><p>Posterior mean ± 1 posterior SD (percentage points)</p>'
    body += f'<span class="badge">{_e(status)}</span></header><section><h2>Samples</h2><div class="table-wrap"><table><thead><tr><th>Sample</th><th>Ce³⁺ result</th><th>Run status</th><th>Report</th></tr></thead><tbody>'
    for entry in results:
        sample = re.split(r"[\\/]", str(entry.get("source", "Unknown sample")))[-1]
        directory = entry.get("output_dir")
        result = {}
        link = ""
        if directory:
            directory = Path(directory).resolve()
            result = _json(directory, "fit_result.json")
            if directory.is_relative_to(root) and (directory / "Results.html").is_file():
                href = quote((directory.relative_to(root) / "Results.html").as_posix(), safe="/")
                link = f'<a href="{href}">Open report</a>'
        body += '<tr><td>' + _e(sample) + '</td><td>' + _e(result_headline(result)) + '</td><td>' + _e(_status(result) if result else entry.get("status", "unknown")) + '</td><td>' + link + '</td></tr>'
    body += '</tbody></table></div></section>'
    return _write_page(root / "Results.html", "Ce 3d XPS batch results", body)
