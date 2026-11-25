"""Meaningful result/report regressions using a small deterministic fixture; no MCMC."""
import csv
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import ce3d_results as reports


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


@pytest.fixture
def result():
    return dict(posterior_available=True, status='quick_complete',
                configuration=dict(source=r'C:\spectra\sample & <script>.xls', model='LA', background='fixed_shirley', satellites=False),
                model=dict(lineshape='LA', parameter_names=['dE', 'sigma_ce', 'x', 'A_ce_total']),
                Ce3_fraction_mean=.2, Ce3_fraction_std=.0123,
                Ce3_percent=dict(median=19.9, q025=17.6, q975=22.6))


def test_headline_uses_mean_and_sd_in_percentage_points_without_changing_median(result):
    assert reports.result_headline(result) == 'Ce³⁺ = 20.00 ± 1.23%'
    assert result['Ce3_percent']['median'] == 19.9
    assert reports.result_headline({'posterior_available': False}) == 'Ce³⁺ result unavailable'
    assert 'SD unavailable' in reports.result_headline(dict(posterior_available=True, Ce3_percent={'median': 19.9}))


def test_full_report_priority_acquisition_scope_and_real_posterior_figures(tmp_path, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(); ax.plot([870, 925], [0, 1]); fig.savefig(tmp_path/'fit_ce3d.png'); plt.close(fig)
    (tmp_path/'fit_ce3d.svg').write_text('<svg/>', encoding='utf-8')
    (tmp_path/'fit_result.json').write_text(json.dumps(result), encoding='utf-8')
    metadata = dict(selected_ce3d_conditions=dict(pass_energy_eV=50, energy_step_eV=.1, scan_count=5),
                    experiment_conditions=dict(pass_energy_eV=150),
                    survey_conditions=dict(pass_energy_eV=200))
    (tmp_path/'experiment_metadata.json').write_text(json.dumps(metadata), encoding='utf-8')
    peaks = [dict(component='Ce4+ v', lineshape='LA', center_eV=882.6, sigma_eV=.5,
                  Lorentzian_HWHM_eV='', LA_kernel_gamma_eV=.3, area=100, fraction_of_total=.2),
             dict(component='Ce3+ v0', lineshape='Voigt', center_eV=880.6, sigma_eV=.5,
                  Lorentzian_HWHM_eV=.4, LA_kernel_gamma_eV='', area=50, fraction_of_total=.1)]
    write_csv(tmp_path/'peak_table.csv', peaks)
    write_csv(tmp_path/'posterior_parameters.csv', [dict(parameter='x', median=.199, mean=.2, std=.0123, q025=.176, q975=.226)])
    rng=np.random.default_rng(17)
    samples=rng.normal(size=(120,4))*.001+np.array([0,.5,.2,100])
    path=reports.write_result_report(tmp_path,result,flat=samples,param_names=['dE','sigma_ce','x','A_ce_total'])
    text=path.read_text(encoding='utf-8')
    for left,right in [('Ce³⁺ = 20.00 ± 1.23%','id="fit"'), ('id="fit"','id="conditions"'),
                       ('id="conditions"','id="peaks"'), ('id="peaks"','id="posterior"')]:
        assert text.index(left)<text.index(right)
    assert 'Posterior mean ± 1 posterior SD (percentage points)' in text
    assert 'Posterior median: 19.90%' in text and '17.60–22.60%' in text
    assert '<dd>50 eV</dd>' in text and '<dd>150 eV</dd>' not in text and '<dd>200 eV</dd>' not in text
    assert 'sample &amp; &lt;script&gt;.xls' in text and '<script>' not in text
    assert 'data:image/png;base64,' in text
    assert 'http://' not in text and 'https://' not in text
    assert '<details><summary>All fitted parameters' in text
    assert 'this status does not establish convergence' in text
    assert 'LA kernel scale is not a Lorentzian HWHM' in text
    assert 'figures/fit_ce3d.svg' in text and 'tables/peak_table.csv' in text
    for name in ['posterior_Ce3.png','posterior_Ce3.svg','posterior_parameters.png','posterior_parameters.svg']:
        p=tmp_path/'figures'/name
        assert p.stat().st_size>1000
    reports.organize_outputs(tmp_path)
    assert (tmp_path/'figures'/'fit_ce3d.png').exists()
    assert (tmp_path/'details'/'fit_result.json').exists()
    assert (tmp_path/'tables'/'peak_table.csv').exists()
    # A second report refresh reads the organized locations and remains usable.
    reports.write_result_report(tmp_path,result,flat=samples,param_names=['dE','sigma_ce','x','A_ce_total'])
    assert 'data:image/png;base64,' in path.read_text(encoding='utf-8')


def test_no_posterior_report_never_invents_values_or_conditions(tmp_path):
    result=dict(posterior_available=False,status='cancelled')
    text=reports.write_result_report(tmp_path,result).read_text(encoding='utf-8')
    assert 'Ce³⁺ result unavailable' in text
    assert 'Stopped; available samples only' in text
    assert '0.00 ±' not in text and '1 posterior SD' not in text
    assert not (tmp_path/'figures').exists()


def test_organization_is_repeatable_and_preserves_unknown_files_and_chain_bytes(tmp_path):
    (tmp_path/'Results.html').write_text('report',encoding='utf-8')
    (tmp_path/'fit_result.json').write_text('{}',encoding='utf-8')
    (tmp_path/'checkpoint_0008000.json').write_text('{}',encoding='utf-8')
    (tmp_path/'batch_summary.json').write_text('{}',encoding='utf-8')
    (tmp_path/'my_notes.txt').write_text('keep',encoding='utf-8')
    (tmp_path/'chain').mkdir(); (tmp_path/'chain'/'sample.npz').write_bytes(b'unchanged chain')
    reports.organize_outputs(tmp_path); reports.organize_outputs(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir())==['Results.html','details','my_notes.txt']
    assert (tmp_path/'details'/'chain'/'sample.npz').read_bytes()==b'unchanged chain'
    assert (tmp_path/'details'/'checkpoint_0008000.json').exists()
    assert (tmp_path/'details'/'batch_summary.json').exists()


def test_organization_refuses_duplicate_destinations_without_overwriting(tmp_path):
    (tmp_path/'details').mkdir()
    (tmp_path/'fit_result.json').write_text('new',encoding='utf-8')
    (tmp_path/'details'/'fit_result.json').write_text('old',encoding='utf-8')
    with pytest.raises(FileExistsError): reports.organize_outputs(tmp_path)
    assert (tmp_path/'details'/'fit_result.json').read_text(encoding='utf-8')=='old'
    assert (tmp_path/'fit_result.json').read_text(encoding='utf-8')=='new'


def test_batch_report_reads_organized_results_and_escapes_names(tmp_path,result):
    sample=tmp_path/'sample one'; (sample/'details').mkdir(parents=True)
    (sample/'details'/'fit_result.json').write_text(json.dumps(result),encoding='utf-8')
    (sample/'Results.html').write_text('sample report',encoding='utf-8')
    text=reports.write_batch_report(tmp_path,[dict(source='a<sample>.xls',output_dir=str(sample),status='quick_complete'),
        dict(source='failed.xls',output_dir=None,status='failed')],'completed_with_errors').read_text(encoding='utf-8')
    assert '20.00 ± 1.23%' in text and 'a&lt;sample&gt;.xls' in text
    assert 'sample%20one/Results.html' in text and 'Ce³⁺ result unavailable' in text
