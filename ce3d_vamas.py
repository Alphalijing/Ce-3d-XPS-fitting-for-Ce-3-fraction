"""VAMAS adapter. Keep acquisition facts separate from fitting assumptions."""
from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import date, datetime
import hashlib
import math
from pathlib import Path
import re

import numpy as np


def json_value(value):
    """Make strict JSON values; absent/nonfinite and VAMAS 1e37 become null."""
    if is_dataclass(value):
        value = asdict(value)
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        if isinstance(value, np.ndarray) and value.ndim == 0:
            return json_value(value.item())
        return [json_value(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return None if not math.isfinite(value) or abs(value) >= 1e36 else float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if isinstance(value, str):
        return value if value.strip() else None
    return value


def source_record(path):
    path = Path(path)
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size}


def _parse(path):
    try:
        from vamas import Vamas
    except ImportError as exc:
        raise ImportError('VAMAS input requires vamas==0.2.0; install requirements.txt.') from exc
    raw = Path(path).read_bytes()
    # A BOM is not part of the identifier; undecodable identifiers are rejected.
    text = raw.decode('utf-8-sig')
    try:
        data = Vamas(text.encode('utf-8'))
    except Exception as exc:
        raise ValueError(f'Cannot parse VAMAS input {Path(path).name}: {exc}') from exc
    if data.header.scan_mode != 'REGULAR':
        raise ValueError('Only REGULAR VAMAS energy scans are supported.')
    return data


def _inventory(data):
    inventory, candidates = [], []
    for i, block in enumerate(data.blocks):
        label = ' '.join(str(getattr(block, key, '')) for key in
                         ('block_identifier', 'species_label', 'transition_or_charge_state_label')).strip()
        variables = block.corresponding_variables
        n = len(variables[0].y_values) if variables else 0
        axis_end = block.x_start + block.x_step * max(0, n - 1)
        inventory.append(dict(index=i, label=label, sample_identifier=block.sample_identifier,
                              technique=block.technique, axis=block.x_label, units=block.x_units,
                              start=block.x_start, step=block.x_step, end=axis_end, point_count=n,
                              ordinates=[dict(index=j, label=v.label, unit=v.unit)
                                         for j, v in enumerate(variables)]))
        if block.technique.strip().upper() == 'XPS' and re.search(r'ce\s*3d', label, re.I):
            candidates.append(i)
    return inventory, candidates


def inspect_vamas(path):
    """List scan/ordinate choices for a GUI without selecting or transforming data."""
    data = _parse(path)
    inventory, candidates = _inventory(data)
    return json_value(dict(source=source_record(path), header=asdict(data.header),
                           block_inventory=inventory, ce3d_candidate_indices=candidates,
                           parser=dict(package='vamas', version='0.2.0', text_encoding='utf-8-sig')))


def read_vamas(path, *, block_index=None, ordinate_index=None, kinetic_offset_ev=None):
    """Return ascending binding energy, unchanged signal, acquisition metadata.

    Indices are zero based. KE requires explicit BE = offset - KE. Automatic
    selection requires exactly one explicitly labelled Ce 3d XPS block; multiple
    signal channels always require explicit selection, including intensity+TXFN.
    """
    path = Path(path)
    data = _parse(path)
    inventory, candidates = _inventory(data)
    if block_index is None:
        if len(candidates) != 1:
            raise ValueError(f'Select a VAMAS block explicitly (zero based); Ce 3d candidate indices: {candidates}. Blocks: {inventory}')
        block_index = candidates[0]
    if isinstance(block_index, bool) or not isinstance(block_index, (int, np.integer)) or not 0 <= block_index < len(data.blocks):
        raise ValueError(f'Invalid VAMAS block index: {block_index}; blocks: {inventory}')
    block_index = int(block_index)
    block = data.blocks[block_index]
    if block.technique.strip().upper() != 'XPS':
        raise ValueError('Selected VAMAS block is not XPS.')
    variables = block.corresponding_variables
    if ordinate_index is None:
        if len(variables) != 1:
            raise ValueError('Multiple VAMAS ordinates: select an ordinate explicitly; no transmission correction is inferred.')
        ordinate_index = 0
    if isinstance(ordinate_index, bool) or not isinstance(ordinate_index, (int, np.integer)) or not 0 <= ordinate_index < len(variables):
        raise ValueError('Invalid VAMAS ordinate index.')
    ordinate_index = int(ordinate_index)
    variable = variables[ordinate_index]
    if 'transmission' in variable.label.lower() or re.search(r'\btxfn\b', variable.label, re.I):
        raise ValueError('The selected ordinate is a transmission function, not a measured intensity channel.')
    y = np.asarray(variable.y_values, dtype=float)
    x = block.x_start + block.x_step * np.arange(y.size, dtype=float)
    if y.size < 20 or not np.isfinite(y).all() or not np.isfinite(x).all() or block.x_step == 0:
        raise ValueError('VAMAS spectrum requires at least 20 finite points and a nonzero energy step.')
    if block.x_units.strip().lower() != 'ev':
        raise ValueError('Only VAMAS energy axes in eV are supported.')
    axis = block.x_label.lower().replace('_', ' ')
    if 'binding' in axis:
        if kinetic_offset_ev is not None:
            raise ValueError('Kinetic offset was supplied for a binding-energy block.')
        transform = dict(formula='BE = stored energy', kinetic_offset_eV=None)
    elif 'kinetic' in axis:
        if kinetic_offset_ev is None or not math.isfinite(kinetic_offset_ev):
            raise ValueError('Kinetic-energy VAMAS needs an explicit kinetic_offset_ev; BE = offset - KE. Confirm the instrument convention first.')
        x = float(kinetic_offset_ev) - x
        transform = dict(formula='BE = explicit offset - stored kinetic energy', kinetic_offset_eV=float(kinetic_offset_ev))
    else:
        raise ValueError(f'Unrecognized VAMAS energy axis: {block.x_label}')
    order = np.argsort(x)
    if np.any(np.diff(x[order]) <= 0):
        raise ValueError('VAMAS binding energies must be unique.')
    acquisition = asdict(block)
    acquisition.pop('corresponding_variables', None)
    metadata = {
        'schema_version': '2.0', 'input_format': 'VAMAS', 'source': source_record(path),
        'parser': dict(package='vamas', version='0.2.0', text_encoding='utf-8-sig'),
        'header': asdict(data.header), 'block_inventory': inventory,
        'ce3d_candidate_indices': candidates,
        'selected_block_index': block_index, 'selected_ordinate_index': ordinate_index,
        'acquisition': acquisition,
        'experiment_conditions': {
            'photon_energy_eV': block.analysis_source_characteristic_energy,
            'analyzer_mode': block.analyzer_mode,
            'pass_energy_or_retard_ratio_or_mass_resolution': block.analyzer_pass_energy_or_retard_ratio_or_mass_res,
            'pass_energy_eV': (block.analyzer_pass_energy_or_retard_ratio_or_mass_res
                               if block.analyzer_mode.strip().upper() == 'FAT' else None),
            'energy_step_eV': abs(block.x_step), 'dwell_time_s': block.signal_collection_time,
            'scan_count': block.num_scans_to_compile_block,
            'total_acquisition_time_s': None,
            'intensity_label': variable.label, 'intensity_unit': variable.unit, 'energy_unit': block.x_units,
            'transmission_correction_applied': None, 'charge_reference': None,
        },
        'energy_transform': transform, 'intensity_transform': 'none',
        'warnings': ['Transmission correction and charge referencing are not inferred from VAMAS presence.',
                     'Missing fields and VAMAS 1e37 sentinel values are represented as null.'],
    }
    return x[order], y[order], json_value(metadata)
