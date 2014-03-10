# -*- coding: utf-8 -*-
import cobra
import cobra.io
from os.path import join
import re
import cPickle as pickle

model_directory = '/Users/zaking/models/'

def add_martinez_pathways(model):
    """Add pathways for lycopene and e-caprolactone

    From:
    1. Martínez I, Zhu J, Lin H, et al. Replacing Escherichia coli
    NAD-dependent glyceraldehyde 3-phosphate dehydrogenase (GAPDH) with a
    NADP-dependent enzyme from Clostridium acetobutylicum facilitates NADPH
    dependent pathways. Metab. Eng. 2008;10(6):352–9.

    """

    new_metabolites = { 'ggpp_c': 'C20H33O7P2',
                        'phyto_c': 'C40H64',
                        'lyco_c': 'C40H56',
                        'lyco_e': 'C40H56' }
    for k, v in new_metabolites.iteritems():
        m = cobra.Metabolite(id=k, formula=v)
        model.add_metabolites([m])
                        
    new_reactions = { #'FPS': { 'ipdp_c': -2,
                      #         'ppi_c': 1,
                      #         'grdp_c': 1 },
                      'CRTE': { 'ipdp_c': -1,
                                'frdp_c': -1,
                                'ggpp_c': 1,
                                'ppi_c': 1 },
                      'CRTB': { 'ggpp_c': -2,
                                'phyto_c': 1,
                                'ppi_c': 2 },
                      'CRTI': { 'phyto_c': -1,
                                'nadp_c': -8,
                                'lyco_c': 1,
                                'nadph_c': 8 },
                      'lycotex': { 'lyco_c': -1,
                                   'lyco_e': 1},
                      'EX_lyco_LPAREN_e_RPAREN_': { 'lyco_e': -1 },
                      'CMHO': { #'o2_c': -1,
                                'nadph_c': -1,
                                'h_c': 1,
                                #'h2o_c': 1,
                                'nadp_c': 1 }
                    }
    # subsytems
    subsystems = { #'FPS': 'Lycopene production',
                   'CRTE': 'Lycopene production', 
                   'CRTB': 'Lycopene production',
                   'CRTI': 'Lycopene production',
                   'lycotex': 'Lycopene production',
                   'EX_lyco_LPAREN_e_RPAREN_': 'Lycopene production',
                   'CMHO': 'Caprolactone production' }    
    # all irreversible
    reversibility = { #'FPS': 0,
                      'CRTE': 0, 
                      'CRTB': 0,
                      'CRTI': 0,
                      'lycotex': 0,
                      'EX_lyco_LPAREN_e_RPAREN_': 0,
                      'CMHO': 0 }
    
    for name, mets in new_reactions.iteritems():
        r = cobra.Reaction(name=name)
        m_obj = {}
        for k, v in mets.iteritems():
            m_obj[model.metabolites.get_by_id(k)] = v
        r.add_metabolites(m_obj)
        r.reversibility = reversibility[name]
        r.subsystem = subsystems[name]
        r.lower_bound = 0
        r.upper_bound = 0
        model.add_reaction(r)
    return model

def add_cap_yeast(model):
    """Add pathway for e-caprolactone to yeast model

    From:
    1. Martínez I, Zhu J, Lin H, et al. Replacing Escherichia coli
    NAD-dependent glyceraldehyde 3-phosphate dehydrogenase (GAPDH) with a
    NADP-dependent enzyme from Clostridium acetobutylicum facilitates NADPH
    dependent pathways. Metab. Eng. 2008;10(6):352–9.

    """
                        
    new_reactions = { 'CMHO': { 'nadph[c]': -1,
                                'h[c]': 1,
                                'nadp[c]': 1 }
                    }
    # subsytems
    subsystems = { 'CMHO': 'Caprolactone production' }    
    # all irreversible
    reversibility = { 'CMHO': 0 }
    
    for name, mets in new_reactions.iteritems():
        r = cobra.Reaction(name=name)
        m_obj = {}
        for k, v in mets.iteritems():
            m_obj[model.metabolites.get_by_id(k)] = v
        r.add_metabolites(m_obj)
        r.reversibility = reversibility[name]
        r.subsystem = subsystems[name]
        r.lower_bound = 0
        r.upper_bound = 0
        model.add_reaction(r)
    return model

def setup_model(model_name, aerobic=True, sur=10, our=10, substrate=None):
    if model_name=='iJO1366':
        path = join(model_directory, 'iJO1366_cobrapy.mat')
        model = cobra.io.load_matlab_model(path)
        o2 = 'EX_o2_e'
        def_substrate = 'EX_glc_e'
        biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
    elif model_name=='iJO1366-heterologous':
        # sbml import for subsytems
        path = join(model_directory, 'iJO1366-heterologous-pathways_cobrapy.pickle')
        with open(path, 'r') as f:
            model = pickle.load(f)
        model = add_martinez_pathways(model)
        o2 = 'EX_o2_LPAREN_e_RPAREN_'
        def_substrate = 'EX_glc_LPAREN_e_RPAREN_'
        biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
    elif model_name=='iAF1260':
        path = join(model_directory, 'Ec_iAF1260_flux1.xml')
        model = cobra.io.read_sbml_model(path)
        o2 = 'EX_o2_e_'
        def_substrate = 'EX_glc_e_'
        biomass_reaction = 'Ec_biomass_iAF1260_core_59p81M'
    elif model_name=='iMM904':
        path = join(model_directory, 'iMM904_cobrapy.mat')
        model = cobra.io.load_matlab_model(path)
        o2 = 'EX_o2(e)'
        def_substrate = 'EX_glc(e)'
        biomass_reaction = 'biomass_SC5_notrace'
    elif model_name=='iMM904-h':
        path = join(model_directory, 'iMM904_cobrapy.mat')
        model = cobra.io.load_matlab_model(path)
        model = add_cap_yeast(model)
        o2 = 'EX_o2(e)'
        def_substrate = 'EX_glc(e)'
        biomass_reaction = 'biomass_SC5_notrace'
    elif model_name=='e_coli_core':
        path = join(model_directory, 'ecoli_core_model_cobrapy.mat')
        model = cobra.io.load_matlab_model(path)
        o2 = 'EX_o2(e)'
        def_substrate = 'EX_glc(e)'
        biomass_reaction = 'Biomass_Ecoli_core_N(w/GAM)-Nmet2'
    else:
        raise Exception('Unrecognized model name %s' % model_name)
    
            
    substrate = def_substrate if substrate is None else substrate
    model.reactions.get_by_id(def_substrate).lower_bound = 0
    model.reactions.get_by_id(substrate).lower_bound = -sur
    if aerobic:
        model.reactions.get_by_id(o2).lower_bound = -our
    else:
        model.reactions.get_by_id(o2).lower_bound = 0
        
    # model specific setup
    if (model_name=='iJO1366' or model_name=='iJO1366-heterologous') and aerobic==False:
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

def turn_on_subsystem(model, subsytem):
    for reaction in model.reactions:
        if reaction.subsystem.strip('_') == subsytem.strip('_'):
            reaction.lower_bound = -1000 if reaction.reversibility else 0
            reaction.upper_bound = 1000
    return model

def swap_reaction(model, reaction, turn_off_old_reaction=True):
    """Swap a reaction.

    Only works for cytosolic NAD(P)H right now.

    TODO accept list of reactions
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
    if sur==0:
        raise Exception('Substrate uptake rate is zero')
    if model.solution.f is None:
        return 0
    return model.solution.f / sur * target_c / substrate_c

def swap_yield(model, reaction_to_swap, target, substrate, biomass,
               min_biomass=lambda mu_max: 0.1, print_results=False, value_for_inf=1000):
    
    gr = max_growth_rate(model, biomass)
    biomass = model.reactions.get_by_id(str(biomass))
    biomass.lower_bound = min_biomass(gr)
    model.optimize(new_objective=model.reactions.get_by_id(str(target)))
    yield_wt = yield_for_product(model, target, substrate)

    # swap model
    model_swap = swap_reaction(model.copy(), reaction_to_swap)
    yield_swap = yield_for_product(model_swap, target, substrate)
    if yield_wt == 0:
        if yield_swap == 0:
            yield_change = 0
        else:
            yield_change = value_for_inf
    else:
        yield_change = (yield_swap - yield_wt) / yield_wt
    
    # print
    if print_results:
        if gr is not None:
            print 'µ  = %.5g' % gr
            print 'min_biomass = %.5g' % biomass.lower_bound
        else:
            print 'INFEASIBLE'
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
