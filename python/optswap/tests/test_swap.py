from optswap import *
import cobra
from numpy.testing import assert_approx_equal

def test_setup_model():
    model, biomass_reaction = setup_model('iJO1366', aerobic=False)
    assert isinstance(model, cobra.Model)
    assert str(model) == 'iJO1366'
    assert model.reactions.get_by_id('CAT').upper_bound == 0

    model, biomass_reaction = setup_model('iMM904', aerobic=False)
    assert isinstance(model, cobra.Model)
    assert str(model) == 'Saccharomyces cerevisiae model iMM904'

def test_swap_reaction():
    model, biomass_reaction = setup_model('iJO1366', aerobic=False)
    model_swap = swap_reaction(model.copy(), 'GAPD', turn_off_old_reaction=True)
    assert model is not model_swap
    assert model_swap.reactions.GAPD.upper_bound == 0
    assert model_swap.reactions.GAPD.lower_bound == 0
    assert model_swap.reactions.GAPD_swap.lower_bound == \
        model.reactions.GAPD.lower_bound
    assert model_swap.reactions.GAPD_swap.upper_bound == \
        model.reactions.GAPD.upper_bound
    nadp = model_swap.metabolites.get_by_id('nadp_c')
    nadph = model_swap.metabolites.get_by_id('nadph_c')
    assert model_swap.reactions.GAPD_swap._metabolites[nadp] == -1
    assert model_swap.reactions.GAPD_swap._metabolites[nadph] == 1
    
    model, biomass_reaction = setup_model('iMM904', aerobic=False)
    model_swap = swap_reaction(model.copy(), 'ALCD2x', turn_off_old_reaction=False)
    assert model_swap.reactions.ALCD2x.upper_bound == \
        model.reactions.ALCD2x.upper_bound
    assert model_swap.reactions.ALCD2x.lower_bound == \
        model.reactions.ALCD2x.lower_bound
    nadp = model_swap.metabolites.get_by_id('nadp[c]')
    nadph = model_swap.metabolites.get_by_id('nadph[c]')
    assert model_swap.reactions.ALCD2x_swap._metabolites[nadp] == -1
    assert model_swap.reactions.ALCD2x_swap._metabolites[nadph] == 1

def test_max_growth_rate():
    model, biomass_reaction = setup_model('iJO1366', aerobic=True)
    model.optimize(new_objective=model.reactions.GAPD)
    gr = max_growth_rate(model, biomass_reaction)
    assert 0.7 < gr and gr < 0.8
    
def test_yield_for_product():
    model, biomass_reaction = setup_model('iJO1366', aerobic=True)
    y = yield_for_product(model, 'EX_etoh_e', 10, 2, 6)
    assert_approx_equal(y, 0.66, significant=2)

def test_swap_yield():
    model, biomass_reaction = setup_model('iJO1366', aerobic=False)
    yields = swap_yield(model, 'GAPD', 'EX_cys__L_e', biomass_reaction,
                        10, 3, 6, print_results=True)
    for y, t in zip(yields, [0.10, 0.18, 0.76]):
        assert_approx_equal(y, t, significant=2)
