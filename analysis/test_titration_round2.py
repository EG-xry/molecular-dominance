"""Tests for documented input fractions, ratio sensitivity and decision ranges."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import run_titration as T

HERE = Path(__file__).resolve().parent


def test_holm_step_down_handles_order_ties_and_cumulative_maximum():
    # Sorted .001,.01,.03,.04 gives .004,.03,.06,.06.
    assert T.holm_adjust([.04, .001, .03, .01]) == pytest.approx([.06, .004, .06, .03])
    assert T.holm_adjust([.01,.01,.01]) == pytest.approx([.03,.03,.03])


def test_all_candidates_kept_in_both_families_even_if_profile_fails():
    df = pd.read_csv(HERE / 'titration_candidate_tests.csv')
    assert len(df) == 52
    assert df.groupby('coding').size().tolist() == [26,26]
    assert (df.Holm_P <= df.Bonferroni_P + 1e-15).all()
    for coding, group in df.groupby('coding'):
        assert set(group.series) == set(df.series)
        assert group.Holm_P.values == pytest.approx(T.holm_adjust(group.nominal_P.values))


def test_grid_range_rejects_favorable_island(monkeypatch):
    monkeypatch.setattr(T,'FRACTION_SHIFTS',(-.05,-.02,0,.02,.05))
    scenarios={}
    for shift in T.FRACTION_SHIFTS:
        for endpoint in (0,.05,.10):
            # At -.02 the decision flips, although it returns at -.05.
            p=.2 if shift==-.02 or endpoint==.05 else .001
            scenarios[T.scenario_key('c',shift,endpoint)]={'one_profile_fit':{'p_bonferroni_7':1},'candidates':{'g':{'p_bonferroni_26':p}}}
    r=T.decision_ranges(scenarios,'c','g','bonferroni',(0,.05,.10))
    assert r['fraction_shift_low_at_zero_endpoint']==0
    assert r['fraction_shift_high_at_zero_endpoint']==.05
    assert r['largest_endpoint_increment_at_zero_fraction_shift']==0
    assert r['joint_rectangle_frontier'][1]['max_endpoint_increment_preserving_every_included_grid_point'] is None


def test_each_retained_record_has_source_n_unit_and_error_bar_status():
    df=pd.read_csv(HERE/'titration_replication.csv')
    assert df.n_recordings.gt(0).all()
    assert df.unit.notna().all()
    assert len(df[df.series.str.startswith('Hosy')])==8
    r=df[(df.series=='Clatot et al.: WT + R104W')]
    assert r.n_recordings.tolist()==[55,19,39,13,13,15]
    assert .1 not in r.plotted_fraction.tolist()
    assert not r.iloc[-1].reported_sem_available
    baseline=pd.read_csv(HERE/'data/titration/katp_baseline_currents.csv')
    assert baseline.nucleotide_free_current_nA.tolist()==[4.4,1,1.3,5.3,1.4,1.6,3.5,1.5]


def test_source_design_supports_all_fits_and_literal_rejections_are_endpoint_sensitive():
    out=json.loads((HERE/'titration.json').read_text())
    assert out['primary_coding']=='plotted_positions'
    for r in out['series'].values():
        fit=r['codings']['plotted_positions']['one_profile_fit']
        assert fit['p_bonferroni_7'] >= .0796
        assert fit['p_holm_7'] >= .0796
    l325r=out['series']['Clatot et al.: WT + L325R']['sensitivity']
    difo=out['series']['Clatot et al.: WT + L325R + difopein']['sensitivity']
    for correction in ('bonferroni','holm'):
        key='p_'+correction+'_7'
        assert l325r[T.scenario_key('stated_cDNA_ratios',0,0)]['one_profile_fit'][key]<.05
        assert l325r[T.scenario_key('stated_cDNA_ratios',0,.01)]['one_profile_fit'][key]>.05
        assert difo[T.scenario_key('stated_cDNA_ratios',0,.02)]['one_profile_fit'][key]<.05
        assert difo[T.scenario_key('stated_cDNA_ratios',0,.05)]['one_profile_fit'][key]>.05


def test_r104w_partner_effect_persists_across_ratio_sensitivity_and_corrections():
    out=json.loads((HERE/'titration.json').read_text())['series']
    for coding in T.CODINGS:
        wt=out['Clatot et al.: WT + R104W']['codings'][coding]
        s460a=out['Clatot et al.: S460A + R104W']['codings'][coding]
        assert .83<wt['profile_estimate'][0]<.88
        assert .40<s460a['profile_estimate'][0]<.48
        for correction in ('bonferroni','holm'):
            key='p_'+correction+'_26'
            assert wt['candidates']['independent subunits (g1 = 0.5)'][key]<.0001
            assert s460a['candidates']['heterodimer inactive (g1 = 1)'][key]<.00042
