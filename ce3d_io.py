"""Input discovery and measured-data preparation for the Ce 3d fitter.

Import metadata describes the selected scan. Missing acquisition facts stay null;
no survey setting, transmission correction or charge reference is inferred.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math
import re

import numpy as np

from ce3d_vamas import json_value, read_vamas, source_record
from shirley_v7 import shirley_v7


@dataclass
class SpectrumInput:
    energy: np.ndarray
    raw: np.ndarray
    background: np.ndarray
    metadata: dict
    notes: list
    source: Path


def _unknown_conditions():
    return dict(photon_energy_eV=None, analyzer_mode=None, pass_energy_eV=None,
                energy_step_eV=None, dwell_time_s=None, total_acquisition_time_s=None,
                scan_count=None, source_gun_type=None, spot_size_um=None, lens_mode=None,
                intensity_label=None, intensity_unit=None, energy_unit=None,
                transmission_correction_applied=None, charge_reference=None)


def _number(value):
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _column_letter(index):
    out = ''
    while index:
        index, rem = divmod(index - 1, 26)
        out = chr(65 + rem) + out
    return out


def _nonempty(value):
    return value is not None and not (isinstance(value, float) and math.isnan(value)) and str(value).strip('\ufeff \t\r\n') != ''


def _blocks(rows, columns):
    """Find contiguous measured rows; isolated numbers in headers are not data."""
    blocks, active = [], []
    for index, row in enumerate(rows):
        values = [_number(row[c]) if c < len(row) else None for c in columns]
        if all(v is not None for v in values):
            active.append((index, values))
        else:
            if active:
                blocks.append(active)
                active = []
    if active:
        blocks.append(active)
    blocks = [b for b in blocks if len(b) >= 20]
    if len(blocks) != 1:
        raise ValueError('Expected one contiguous spectrum with at least 20 numeric points; found ' + str(len(blocks)) + '. Select a single Ce 3d scan.')
    selected = blocks[0]
    # A partly missing/non-numeric row touching the spectrum is an error, not
    # a reason to silently shorten it. Unit/title rows contain no numeric fields.
    for adjacent in (selected[0][0] - 1, selected[-1][0] + 1):
        if 0 <= adjacent < len(rows):
            row = rows[adjacent]
            if any(c < len(row) and _number(row[c]) is not None for c in columns):
                raise ValueError('Incomplete numeric spectrum row: every energy/intensity/background value must be finite and numeric.')
    return selected


def _quantity(text, unit_pattern):
    match = re.search(r'([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*(' + unit_pattern + r')\b', str(text), re.I)
    return float(match.group(1)) if match else None


def _time_seconds(text):
    """Parse explicit time units only (e.g. the instrument's I4 export)."""
    found = re.findall(r'([-+]?\d+(?:\.\d*)?)\s*(hours?|hrs?|h|minutes?|mins?|min|m|seconds?|secs?|sec|s)\b', str(text), re.I)
    if not found:
        return None
    factors = {'h': 3600., 'm': 60., 's': 1.}
    return sum(float(value) * factors[unit.lower()[0]] for value, unit in found)


def _sheet_metadata(name, rows):
    # A/C are the instrument export's measured energy/intensity columns. The
    # whole metadata area is read, including H:I and J and experiment comments.
    try:
        first_data = _blocks(rows, (0, 2))[0][0]
    except ValueError:
        first_data = 20
    header_count = min(20, first_data)
    raw_cells, pairs = [], []
    for i, row in enumerate(rows[:header_count]):
        for j, value in enumerate(row):
            if not _nonempty(value):
                continue
            cell = f'{_column_letter(j + 1)}{i + 1}'
            raw_cells.append(dict(cell=cell, row=i + 1, column=j + 1, value=value))
            if isinstance(value, str) and _number(value) is None and j + 1 < len(row) and _nonempty(row[j + 1]):
                adjacent = row[j + 1]
                # Adjacency, values and explicit units are evidence, not a claim
                # that every adjacent text cell is an acquisition parameter.
                unit_match = re.search(r'(eV|µm|μm|um|Counts\s*/\s*s|cps|mins?\b|secs?\b)', str(adjacent), re.I)
                pairs.append(dict(label=value.strip(), value=adjacent,
                                  label_cell=cell, value_cell=f'{_column_letter(j + 2)}{i + 1}',
                                  unit=unit_match.group(1) if unit_match else None))
    conditions = _unknown_conditions()
    for p in pairs:
        label = p['label'].strip().lower()
        val = p['value']
        if 'total acquisition time' in label or 'total acquire time' in label:
            conditions['total_acquisition_time_s'] = _time_seconds(val)
        elif label in ('number of scans', 'scans', 'scan count'):
            value = _number(val)
            conditions['scan_count'] = int(value) if value is not None and math.isfinite(value) and value.is_integer() else None
        elif label == 'source gun type':
            conditions['source_gun_type'] = str(val)
        elif label == 'spot size':
            conditions['spot_size_um'] = _quantity(val, r'µm|μm|um')
        elif label == 'lens mode':
            conditions['lens_mode'] = str(val)
        elif label in ('analyser mode', 'analyzer mode'):
            conditions['analyzer_mode'] = str(val)
            match = re.search(r'pass\s+energy\s*[:=]?\s*([-+]?\d+(?:\.\d*)?)\s*eV\b', str(val), re.I)
            if match:
                conditions['pass_energy_eV'] = float(match.group(1))
        elif label in ('energy step size', 'energy step'):
            conditions['energy_step_eV'] = _quantity(val, r'eV')
        elif 'dwell time' in label:
            conditions['dwell_time_s'] = _time_seconds(val)
        elif 'photon energy' in label or 'source energy' in label:
            conditions['photon_energy_eV'] = _quantity(val, r'eV')
    # Units are only read from actual unit cells, never deduced from magnitudes.
    for row in rows[:header_count]:
        if len(row) > 0 and isinstance(row[0], str) and row[0].strip().lower() == 'ev':
            conditions['energy_unit'] = 'eV'
            if len(row) > 2 and _nonempty(row[2]):
                conditions['intensity_unit'] = str(row[2]).strip()
                conditions['intensity_label'] = str(row[2]).strip()
    return json_value(dict(sheet=name, raw_metadata_cells=raw_cells, paired_metadata=pairs,
                           experiment_conditions=conditions))


def _read_excel(path, sheet, provided):
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError('Excel input requires pandas, xlrd (.xls) and openpyxl (.xlsx).') from exc
    book = pd.ExcelFile(path)
    if sheet is None:
        candidates = [s for s in book.sheet_names if re.search(r'ce\s*3d', s, re.I)]
        if len(candidates) == 1:
            selected = candidates[0]
        elif len(book.sheet_names) == 1:
            selected = book.sheet_names[0]
        else:
            raise ValueError(f'Select an Excel sheet explicitly; Ce 3d sheet is ambiguous. Sheets: {book.sheet_names}')
    elif isinstance(sheet, int) and not isinstance(sheet, bool):
        if not 0 <= sheet < len(book.sheet_names):
            raise ValueError(f'Invalid Excel sheet index: {sheet}')
        selected = book.sheet_names[sheet]
    elif sheet in book.sheet_names:
        selected = sheet
    else:
        raise ValueError(f'Unknown Excel sheet: {sheet}; sheets: {book.sheet_names}')
    sheets, selected_rows = {}, None
    for name in book.sheet_names:
        frame = pd.read_excel(book, sheet_name=name, header=None, dtype=object)
        rows = frame.where(pd.notna(frame), None).values.tolist()
        sheets[name] = _sheet_metadata(name, rows)
        if name == selected:
            selected_rows = rows
    book.close()
    columns = (0, 2, 4) if provided else (0, 2)
    block = _blocks(selected_rows, columns)
    values = np.asarray([row for _, row in block], dtype=float)
    metadata = dict(input_format=path.suffix.lower().lstrip('.'),
                    parser=dict(package='pandas', version=pd.__version__),
                    sheet_inventory=list(sheets), sheets=sheets, selected_sheet=selected,
                    raw_metadata_cells=sheets[selected]['raw_metadata_cells'],
                    paired_metadata=sheets[selected]['paired_metadata'],
                    experiment_conditions=sheets[selected]['experiment_conditions'].copy(),
                    selected_ce3d_conditions=sheets[selected]['experiment_conditions'].copy(),
                    data_columns={'energy': 'A', 'raw': 'C', 'background': 'E' if provided else None},
                    data_rows_1based=[block[0][0] + 1, block[-1][0] + 1])
    return values[:, 0], values[:, 1], values[:, 2] if provided else None, metadata


def _read_text(path, provided):
    raw = path.read_bytes()
    try:
        text = raw.decode('utf-8-sig')
        encoding = 'utf-8-sig'
    except UnicodeDecodeError:
        text = raw.decode('gb18030')
        encoding = 'gb18030'
    lines = text.splitlines()
    needed = 3 if provided else 2
    possible = []
    for delimiter in ('\t', ',', ';', None):
        rows = []
        for line in lines:
            if delimiter is None:
                row = re.split(r'\s+', line.strip()) if line.strip() else []
            else:
                row = next(csv.reader([line], delimiter=delimiter))
            rows.append([v.strip() for v in row])
        try:
            block = _blocks(rows, tuple(range(needed)))
        except ValueError:
            continue
        # Only a 2/3-column table is accepted, avoiding instrument header numbers
        # or arbitrary many-column parameter/output tables being silently loaded.
        if not all(needed <= len(rows[i]) <= 3 for i, _ in block):
            continue
        possible.append((len(block), delimiter, rows, block))
    if not possible:
        raise ValueError('No contiguous 2-column raw or 3-column spectrum found (at least 20 numeric points required). Provided mode requires column 3.')
    _, delimiter, rows, block = max(possible, key=lambda item: item[0])
    values = np.asarray([row for _, row in block], dtype=float)
    header_lines = lines[:block[0][0]]
    conditions = _unknown_conditions()
    # A stated energy/intensity unit can be retained; absent units remain null.
    header = '\n'.join(header_lines)
    if re.search(r'\beV\b', header, re.I):
        conditions['energy_unit'] = 'eV'
    unit = re.search(r'counts\s*/\s*s|counts\s+per\s+second|\bcps\b', header, re.I)
    if unit:
        conditions['intensity_unit'] = unit.group(0)
        conditions['intensity_label'] = unit.group(0)
    metadata = dict(input_format=path.suffix.lower().lstrip('.'),
                    parser=dict(text_encoding=encoding, delimiter='whitespace' if delimiter is None else delimiter),
                    raw_header_lines=header_lines, experiment_conditions=conditions,
                    data_columns={'energy': 1, 'raw': 2, 'background': 3 if provided else None},
                    data_rows_1based=[block[0][0] + 1, block[-1][0] + 1])
    return values[:, 0], values[:, 1], values[:, 2] if provided else None, metadata


def load_spectrum(path, *, background_mode='fixed_shirley', vamas_block=None,
                  vamas_ordinate=None, kinetic_offset_ev=None, sheet=None):
    """Auto-detect a measured Ce 3d scan and prepare its fixed background.

    Energy is returned ascending and unique. The original lower bound is kept;
    only the common 925 eV upper limit is applied before the V7 background.
    Explicit ``provided`` mode uses text column 3 or Excel column E unchanged.
    """
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if background_mode not in ('fixed_shirley', 'provided'):
        raise ValueError("background_mode must be 'fixed_shirley' or 'provided'.")
    provided = background_mode == 'provided'
    suffix = path.suffix.lower()
    notes = []
    if suffix in ('.vms', '.vamas'):
        if provided:
            raise ValueError('VAMAS contains measured signal; no ordinate is automatically a background. Use fixed_shirley mode.')
        energy, signal, metadata = read_vamas(path, block_index=vamas_block,
                                             ordinate_index=vamas_ordinate,
                                             kinetic_offset_ev=kinetic_offset_ev)
        supplied = None
    elif suffix in ('.xls', '.xlsx'):
        if any(v is not None for v in (vamas_block, vamas_ordinate, kinetic_offset_ev)):
            raise ValueError('VAMAS selection and kinetic offset apply only to VAMAS input.')
        energy, signal, supplied, metadata = _read_excel(path, sheet, provided)
    elif suffix in ('.csv', '.txt', '.dat', '.xy'):
        if any(v is not None for v in (vamas_block, vamas_ordinate, kinetic_offset_ev, sheet)):
            raise ValueError('Sheet/VAMAS options do not apply to text input.')
        energy, signal, supplied, metadata = _read_text(path, provided)
    else:
        raise ValueError(f'Unsupported input extension: {suffix}; use .xls/.xlsx/.csv/.txt/.dat/.xy/.vms/.vamas.')
    energy, signal = np.asarray(energy, float), np.asarray(signal, float)
    if energy.ndim != 1 or signal.shape != energy.shape or energy.size < 20:
        raise ValueError('A Ce 3d spectrum requires at least 20 paired points.')
    if not np.isfinite(energy).all() or not np.isfinite(signal).all():
        raise ValueError('Energy and measured intensity must contain only finite numbers.')
    order = np.argsort(energy)
    was_sorted = bool(np.array_equal(order, np.arange(energy.size)))
    energy, signal = energy[order], signal[order]
    if np.any(np.diff(energy) <= 0):
        raise ValueError('Duplicate energy values are not allowed; choose a single scan without repeated energies.')
    if energy[-1] - energy[0] > 150:
        raise ValueError('Survey or overly wide input rejected: select a narrow Ce 3d scan (width at most 150 eV).')
    if energy[0] > 879 or energy[-1] < 918:
        raise ValueError('The selected Ce 3d scan must span at least 879–918 eV.')
    original_range = [float(energy[0]), float(energy[-1])]
    original_count = energy.size
    keep = energy <= 925
    if provided:
        supplied = np.asarray(supplied, float)
        if supplied.shape != order.shape or not np.isfinite(supplied).all():
            raise ValueError('Provided background requires one finite value for every measured point.')
        background = supplied[order][keep]
    energy, signal = energy[keep], signal[keep]
    if energy.size < 20:
        raise ValueError('Fewer than 20 points remain after the 925 eV crop.')
    if not was_sorted:
        notes.append('Energy and all signal columns sorted together into ascending binding energy.')
    cropped = int(original_count - energy.size)
    if cropped:
        notes.append(f'Removed {cropped} points above 925 eV; original lower energy retained.')
    if provided:
        bg_info = {'name': 'provided background', 'method': 'input values unchanged', 'jointly_fitted': False}
    else:
        background, info = shirley_v7(energy, signal, endpoint_mode='ends', smoothing_window=1)
        bg_info = {'name': 'fixed Shirley background', 'jointly_fitted': False, **info}
        notes.append('Computed fixed Shirley background (V7; endpoint_mode=ends; smoothing_window=1).')
    metadata.update(schema_version='2.0', source=source_record(path))
    metadata.setdefault('experiment_conditions', _unknown_conditions())
    metadata['preprocessing'] = dict(background_mode=background_mode, background=bg_info,
                                     input_range_eV=original_range,
                                     fit_range_eV=[float(energy[0]), float(energy[-1])],
                                     lower_energy_cut_eV=None, upper_energy_cut_eV=925.0,
                                     original_point_count=int(original_count), fitted_point_count=int(energy.size),
                                     points_removed_above_925_eV=cropped,
                                     energy_sorted_ascending=True, intensity_transform='none')
    metadata['notes'] = notes
    return SpectrumInput(energy=energy, raw=signal, background=background,
                         metadata=json_value(metadata), notes=notes, source=path)
