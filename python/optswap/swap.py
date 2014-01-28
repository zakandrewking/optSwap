# -*- coding: utf-8 -*-
import cobra
import cobra.io
from os.path import join
import re

model_directory = '/Users/zaking/models/'

def setup_model(model_name, aerobic=True, min_biomass=0.1,
                sur=10, our=10, substrate=None):
    if model_name=='iJO1366':
        path = join(model_directory, 'iJO1366_cobrapy.mat')
        o2 = 'EX_o2_e'
        def_substrate = 'EX_glc_e'
        biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
    elif model_name=='iMM904':
        path = join(model_directory, 'iMM904_cobrapy.mat')
        o2 = 'EX_o2(e)'
        def_substrate = 'EX_glc(e)'
        biomass_reaction = 'biomass_SC5_notrace'
    else:
        raise Exception('Unrecognized model name %s' % model_name)
    
    model = cobra.io.load_matlab_model(path)
            
    substrate = def_substrate if substrate is None else substrate
    model.reactions.get_by_id(def_substrate).lower_bound = 0
    model.reactions.get_by_id(substrate).lower_bound = -sur
    if aerobic:
        model.reactions.get_by_id(o2).lower_bound = -our
    else:
        model.reactions.get_by_id(o2).lower_bound = 0
        
    # model specific setup
    if model_name=='iJO1366' and aerobic==False:
        for r in ['CAT', 'SPODM', 'SPODMpp']:
            model.reactions.get_by_id(r).lower_bound = 0
            model.reactions.get_by_id(r).upper_bound = 0

    if model_name=='iMM904' and aerobic==False:
        necessary_ex = ['EX_ergst(e)', 'EX_zymst(e)', 'EX_hdcea(e)',
                        'EX_ocdca(e)', 'EX_ocdcea(e)', 'EX_ocdcya(e)']
        for r in necessary_ex:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

    return model, biomass_reaction

def swap_reaction(model, reaction, turn_off_old_reaction=True):
    """Swap a reaction.

    Only works for cytosolic NAD(P)H right now.

    TODO replace the replace dictionary (lol) with a regex.
    
    """
    replace = { 'nad[c]': 'nadp[c]',
                'nadh[c]': 'nadph[c]',
                'nadp[c]': 'nad[c]',
                'nadph[c]': 'nadh[c]',
                'nad_c': 'nadp_c',
                'nadh_c': 'nadph_c',
                'nadp_c': 'nad_c',
                'nadph_c': 'nadh_c' }
    old_rxn = model.reactions.get_by_id(str(reaction))
    swap = cobra.Reaction(name=str(old_rxn)+'_swap')
    mets = old_rxn._metabolites
    mets = dict([(replace[str(k)], v) if str(k) in replace else (str(k), v)
                 for k, v in mets.iteritems()])
    mets = dict([(model.metabolites.get_by_id(k), v)
                 for k, v in mets.iteritems()])
    swap.add_metabolites(mets)
    swap.reversibility = old_rxn.reversibility
    swap.lower_bound = old_rxn.lower_bound
    swap.upper_bound = old_rxn.upper_bound
    model.add_reaction(swap)
    if turn_off_old_reaction:
        old_rxn.lower_bound = 0
        old_rxn.upper_bound = 0
    return model

def max_growth_rate(model, biomass_reaction):
    model.optimize(new_objective=model.reactions.get_by_id(str(biomass_reaction)))
    return model.solution.f

def yield_for_product(model, target, substrate):
    target_rxn = model.reactions.get_by_id(str(target))
    substrate_rxn = model.reactions.get_by_id(str(substrate))
    model.optimize(new_objective=target_rxn)
    sur = -substrate_rxn.lower_bound
    target_c = carbons_for_exchange_reaction(target_rxn)
    substrate_c = carbons_for_exchange_reaction(substrate_rxn)
    return model.solution.f / sur * target_c / substrate_c

def swap_yield(model, reaction_to_swap, target, substrate, biomass,
               min_biomass=lambda mu_max: 0.1, print_results=False):
    
    gr = max_growth_rate(model, biomass)
    biomass = model.reactions.get_by_id(str(biomass))
    biomass.lower_bound = min_biomass(gr)
    model.optimize(new_objective=model.reactions.get_by_id(str(target)))
    yield_wt = yield_for_product(model, target, substrate)

    # swap model
    model_swap = swap_reaction(model.copy(), reaction_to_swap)
    yield_swap = yield_for_product(model_swap, target, substrate)
    yield_change = (yield_swap - yield_wt) / yield_wt
    
    # print
    if print_results:
        print 'µ  = %.5g' % model.solution.f
        print 'wt yield\t%.5g' % yield_wt
        print 'swap yield\t%.5g' % yield_swap
        print 'change\t%.5g' % yield_change
    return yield_wt, yield_swap, yield_change

def carbons_for_exchange_reaction(reaction):
    if len(reaction._metabolites) > 1:
        raise Exception('%s not an exchange reaction' % str(reaction))

    metabolite = reaction._metabolites.iterkeys().next()
    match = re.match(r'C([0-9]+)', str(metabolite.formula))
    try:
        return int(match.group(1))
    except AttributeError:
        return 0
