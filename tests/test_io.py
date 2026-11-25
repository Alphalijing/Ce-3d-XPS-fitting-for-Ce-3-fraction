"""Measured import regressions; these tests do not run MCMC."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ce3d_io import load_spectrum
from ce3d_vamas import inspect_vamas, json_value, read_vamas


@pytest.fixture
def text_scan(tmp_path):
    def make(*, extension='.csv', delimiter=',', background=False, descending=False, energy=None, header='Binding energy (eV),Counts / s'):
        x = np.arange(870., 931.) if energy is None else np.asarray(energy, float)
        y = 10. + 0.1 * (x - 870.) + 50. * np.exp(-((x - 892.) / 5.) ** 2)
        bg = np.full(x.size, 7.25)
        rows = np.column_stack((x, y, bg)) if background else np.column_stack((x, y))
        if descending:
            rows = rows[::-1]
        p = tmp_path / ('spectrum' + extension)
        with p.open('w', encoding='utf-8') as f:
            f.write(header + '\n')
            for row in rows:
                f.write(delimiter.join(str(v) for v in row) + '\n')
        return p, x, y
    return make


def test_real_xls_conditions_I4_and_no_survey_leak():
    p = ROOT / 'examples/shang_yong_CeO2_raw.xls'
    s = load_spectrum(p)
    assert s.metadata['selected_sheet'] == 'Ce3d Scan'
    assert len(s.energy) == 1092
    assert s.energy[0] == pytest.approx(870.45)
    assert s.energy[-1] == pytest.approx(925.)
    assert s.raw[0] == pytest.approx(53272.3)
    assert s.raw[-1] == pytest.approx(104307.)
    conditions = s.metadata['experiment_conditions']
    assert conditions['total_acquisition_time_s'] == pytest.approx(480.4)
    assert conditions['scan_count'] == 8
    assert conditions['pass_energy_eV'] == 50.
    assert conditions['energy_step_eV'] == 0.05
    assert conditions['intensity_unit'] == 'Counts / s'
    survey = s.metadata['sheets']['XPS Survey']['experiment_conditions']
    assert survey['pass_energy_eV'] == 150.
    assert survey['total_acquisition_time_s'] == pytest.approx(136.1)
    assert conditions['photon_energy_eV'] is None  # Al K Alpha alone is not an explicit numeric acquisition field.
    assert conditions['transmission_correction_applied'] is None
    assert conditions['charge_reference'] is None
    cells = {r['cell']: r['value'] for r in s.metadata['raw_metadata_cells']}
    assert cells['H4'] == 'Total acquisition time'
    assert cells['I4'] == '8 mins 0.4 secs'
    pair = next(r for r in s.metadata['paired_metadata'] if r['value_cell'] == 'I4')
    assert pair['label_cell'] == 'H4'
    assert pair['unit'] == 'mins'
    assert s.metadata['source']['sha256'] == hashlib.sha256(p.read_bytes()).hexdigest()
    json.dumps(s.metadata, allow_nan=False)


def test_real_xls_provided_E_is_used():
    s = load_spectrum(ROOT / 'examples/shang_yong_CeO2_raw.xls', background_mode='provided')
    assert s.metadata['data_columns']['background'] == 'E'
    assert s.metadata['preprocessing']['background']['name'] == 'provided background'
    assert np.isfinite(s.background).all()
    assert s.background.shape == s.raw.shape


def test_excel_Ce3d_without_survey_is_supported(tmp_path):
    import openpyxl
    book = openpyxl.Workbook()
    ws = book.active
    ws.title = 'Ce3d Scan'
    ws['H4'], ws['I4'] = 'Total acquisition time', '2 mins 3.5 secs'
    ws['H9'], ws['I9'] = 'Analyser Mode', 'CAE : Pass Energy 25.0 eV'
    ws['A15'], ws['C15'] = 'eV', 'Counts / s'
    for i, x in enumerate(np.arange(930., 869., -1.), 17):
        ws.cell(i, 1, x)
        ws.cell(i, 3, x * 2)
    p = tmp_path / 'single_sheet.xlsx'
    book.save(p)
    s = load_spectrum(p)
    assert s.metadata['sheet_inventory'] == ['Ce3d Scan']
    assert s.metadata['experiment_conditions']['total_acquisition_time_s'] == 123.5
    assert s.metadata['experiment_conditions']['pass_energy_eV'] == 25.
    assert np.allclose(s.raw, 2 * s.energy)


@pytest.mark.parametrize('extension,delimiter', [('.csv', ','), ('.txt', '\t'), ('.dat', ' '), ('.xy', '  ')])
def test_auto_text_formats_and_numeric_header(text_scan, extension, delimiter):
    p, x, y = text_scan(extension=extension, delimiter=delimiter,
                         header='Acquired in 2026\n1\nEnergy (eV) Intensity (cps)')
    s = load_spectrum(p)
    assert np.array_equal(s.energy, x[x <= 925])
    assert np.allclose(s.raw, y[x <= 925])
    assert s.metadata['experiment_conditions']['energy_unit'] == 'eV'
    bg = s.metadata['preprocessing']['background']
    assert bg['name'] == 'fixed Shirley background'
    assert bg['endpoint_mode'] == 'ends'
    assert bg['smoothing_window_points'] == 1
    assert s.background[0] == s.raw[0]
    assert s.background[-1] == s.raw[-1]
    assert np.all(s.background <= s.raw + 1e-8)
    assert s.metadata['preprocessing']['lower_energy_cut_eV'] is None


def test_descending_three_columns_sort_together_crop_and_explicit_mode(text_scan):
    p, x, y = text_scan(background=True, descending=True)
    s = load_spectrum(p, background_mode='provided')
    keep = x <= 925
    assert np.array_equal(s.energy, x[keep])
    assert np.allclose(s.raw, y[keep])
    assert np.array_equal(s.background, np.full(sum(keep), 7.25))
    default = load_spectrum(p)
    assert default.background[0] == default.raw[0]
    assert default.background[-1] == default.raw[-1]
    assert default.metadata['preprocessing']['background']['name'] == 'fixed Shirley background'


def test_unknown_units_stay_null(text_scan):
    p, _, _ = text_scan(header='Energy Intensity')
    conditions = load_spectrum(p).metadata['experiment_conditions']
    assert conditions['energy_unit'] is None
    assert conditions['intensity_unit'] is None
    assert conditions['pass_energy_eV'] is None


@pytest.mark.parametrize('case,match', [('duplicates', 'Duplicate'), ('survey', 'Survey'), ('incomplete', 'span'), ('nonfinite', 'finite')])
def test_rejects_invalid_or_wrong_scan(text_scan, case, match):
    x = np.arange(870., 931.)
    if case == 'duplicates':
        x[10] = x[11]
    elif case == 'survey':
        x = np.arange(0., 1200., 2.)
    elif case == 'incomplete':
        x = np.arange(885., 925.)
    elif case == 'nonfinite':
        x[30] = np.nan
    p, _, _ = text_scan(energy=x)
    with pytest.raises(ValueError, match=match):
        load_spectrum(p)


def test_provided_mode_requires_background_column(text_scan):
    p, _, _ = text_scan()
    with pytest.raises(ValueError, match='column 3'):
        load_spectrum(p, background_mode='provided')


def test_xls_survey_selection_is_rejected():
    with pytest.raises(ValueError, match='Survey'):
        load_spectrum(ROOT / 'examples/shang_yong_CeO2_raw.xls', sheet='XPS Survey')


VMS = ROOT / 'examples/GDC10_published.vms'


def test_real_vamas_inventory_and_ambiguities():
    inventory = inspect_vamas(VMS)
    assert len(inventory['block_inventory']) == 61
    assert inventory['ce3d_candidate_indices'] == [2, 12, 17, 22, 27, 32, 37, 41, 46, 51, 60]
    block = inventory['block_inventory'][12]
    assert block['axis'] == 'Kinetic Energy'
    assert block['point_count'] == 1001
    assert [o['label'] for o in block['ordinates']] == ['Intensity', 'Transmission']
    with pytest.raises(ValueError, match='Select a VAMAS block explicitly'):
        load_spectrum(VMS)
    with pytest.raises(ValueError, match='Multiple VAMAS ordinates'):
        load_spectrum(VMS, vamas_block=12)
    with pytest.raises(ValueError, match='explicit kinetic_offset_ev'):
        load_spectrum(VMS, vamas_block=12, vamas_ordinate=0)
    with pytest.raises(ValueError, match='transmission function'):
        load_spectrum(VMS, vamas_block=12, vamas_ordinate=1, kinetic_offset_ev=1486.6)


def test_real_vamas_kinetic_conversion_and_metadata():
    s = load_spectrum(VMS, vamas_block=12, vamas_ordinate=0, kinetic_offset_ev=1486.6)
    assert s.energy.size == 1001
    assert np.all(np.diff(s.energy) > 0)
    assert s.energy[0] == pytest.approx(874.)
    assert s.energy[-1] == pytest.approx(924.)
    assert s.metadata['selected_block_index'] == 12
    assert s.metadata['experiment_conditions']['pass_energy_eV'] == 27.
    assert s.metadata['experiment_conditions']['dwell_time_s'] == 0.3
    assert s.metadata['experiment_conditions']['transmission_correction_applied'] is None
    assert s.metadata['energy_transform']['kinetic_offset_eV'] == 1486.6
    assert s.metadata['preprocessing']['background']['name'] == 'fixed Shirley background'
    assert s.metadata['acquisition']['sample_identifier'] == '+800 mV'
    json.dumps(s.metadata, allow_nan=False)


def test_binding_axis_does_not_accept_kinetic_offset(tmp_path):
    # Explicit synthetic edit of the real file's label, for axis-policy testing.
    p = tmp_path / 'binding_label_only_synthetic.vamas'
    p.write_bytes(VMS.read_bytes().replace(b'Kinetic Energy', b'Binding Energy'))
    with pytest.raises(ValueError, match='binding-energy block'):
        read_vamas(p, block_index=12, ordinate_index=0, kinetic_offset_ev=1486.6)
    x, _, metadata = read_vamas(p, block_index=12, ordinate_index=0)
    assert x[0] == pytest.approx(562.6)
    assert x[-1] == pytest.approx(612.6)
    assert metadata['energy_transform']['formula'] == 'BE = stored energy'


def test_json_value_is_strict_json_for_numpy_and_missing():
    result = json_value(dict(path=Path('sample'), flag=np.bool_(True), value=np.array(3.5),
                             missing=[np.nan, np.inf, -1e37, 1e37], count=np.int64(2)))
    assert result == dict(path='sample', flag=True, value=3.5, missing=[None] * 4, count=2)
    json.dumps(result, allow_nan=False)


def test_excel_missing_intensity_does_not_silently_truncate(tmp_path):
    import openpyxl
    book = openpyxl.Workbook()
    ws = book.active
    ws.title = 'Ce3d Scan'
    for i, x in enumerate(np.arange(870., 931.), 2):
        ws.cell(i, 1, x)
        if x != 925.:
            ws.cell(i, 3, 100. + x)
    p = tmp_path / 'missing_point.xlsx'
    book.save(p)
    with pytest.raises(ValueError, match='finite and numeric'):
        load_spectrum(p)
