from optswap import *
import cobra
from numpy.testing import assert_approx_equal

def test_setup_model():
    model, biomass_reaction = setup_model('iJO1366', aerobic=False)
    assert isinstance(model, cobra.Model)
    assert str(model) == 'iJO1366'
    assert model.reactions.get_by_id('CAT').upper_bound == 0
    
    model, biomass_reaction = setup_model('iJO1366-heterologous')
    assert isinstance(model, cobra.Model)
    assert str(model) == 'MODELID_4188667'
    assert isinstance(model.metabolites.lyco_c, cobra.Metabolite)
    assert isinstance(model.reactions.FPS, cobra.Reaction)
    assert model.reactions.get_by_id('FPS').lower_bound == 0
    assert model.reactions.get_by_id('FPS').upper_bound == 0
    model = turn_on_subsystem(model, 'Lycopene production')
    assert yield_for_product(model, 'EX_lyco_LPAREN_e_RPAREN_', 'EX_glc_LPAREN_e_RPAREN_') > 0
    model = turn_on_subsystem(model, 'Caprolactone production')
    model.optimize(new_objective=model.reactions.get_by_id('CMHO'))
    assert model.solution.f > 0

    model, biomass_reaction = setup_model('iMM904', aerobic=False)
    assert isinstance(model, cobra.Model)
    assert str(model) == 'Saccharomyces cerevisiae model iMM904'

def test_turn_on_subsystem():
    model, biomass_reaction = setup_model('iJO1366-heterologous')
    model = turn_on_subsystem(model, '1,3-Propanediol production')
    assert model.reactions.get_by_id('13PPDH2').upper_bound == 1000
    assert model.reactions.get_by_id('13PPDH2').lower_bound == -1000
    model = turn_on_subsystem(model, '1,4-Butanediol production')
    assert model.reactions.get_by_id('BTDP2').upper_bound == 1000
    assert model.reactions.get_by_id('BTDP2').lower_bound == 0
    
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
    y = yield_for_product(model, 'EX_etoh_e', 'EX_glc_e')
    assert_approx_equal(y, 0.66, significant=2)

    # TODO fix
    # model, biomass_reaction = setup_model('iJO1366', aerobic=False, substrate='EX_xyl__D_e')
    # y = yield_for_product(model, 'EX_cys__L_e', 'EX_xyl__D_e')
    # assert_approx_equal(y, 0.11, significant=2)
    
    model, biomass_reaction = setup_model('iJO1366', aerobic=False, substrate='EX_xyl__D_e')
    model.reactions.get_by_id(biomass_reaction).lower_bound = 10
    y = yield_for_product(model, 'EX_cys__L_e', 'EX_xyl__D_e')
    assert y == 0

def test_swap_yield():
    model, biomass_reaction = setup_model('iJO1366', aerobic=False)
    yields = swap_yield(model, 'GAPD', 'EX_cys__L_e', 'EX_glc_e', biomass_reaction,
                        print_results=True)
    for y, t in zip(yields, [0.10, 0.18, 0.76]):
        assert_approx_equal(y, t, significant=2)

    yields = swap_yield(model, 'GAPD', 'EX_cys__L_e', 'EX_glc_e', biomass_reaction,
                        print_results=True, min_biomass = lambda x: 10)
    assert yields == (0, 0, 0)
            
def test_carbons_for_exchange_reaction():
    model, biomass_reaction = setup_model('iJO1366')
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_glc_e')) == 6
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_cys__L_e')) == 3
    assert carbons_for_exchange_reaction(model.reactions.get_by_id('EX_nh4_e')) == 0
    
