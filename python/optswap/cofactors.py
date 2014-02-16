from cobra.oven.aliebrahim.fluxAnalysis import optimize_minimal_flux

class NoFluxError(Exception):
    pass

def track_cofactor(model, target, cofactor):
    # get the metabolite
    cofactor_met = model.metabolites.get_by_id(cofactor)
    # get the reaction
    target_reaction = model.reactions.get_by_id(target)

    #change objective
    for r in model.reactions:
        r.objective_coefficient = 0
    target_reaction.objective_coefficient = 1
    pfba_model = model.copy()
    optimize_minimal_flux(pfba_model)
    flux = pfba_model.solution.x_dict
    if flux is None:
        raise NoFluxError

    flux = fix_flux(flux)
    
    # find the reactions using that cofactor
    reactions_using_cofactor = [(str(r), r._metabolites[cofactor_met], flux[str(r)]) for r in model.reactions
                                if cofactor_met in r._metabolites
                                and str(r) in flux
                                and abs(flux[str(r)]) > 0]
    cofactor_production = sum([f * st for (r, st, f) in reactions_using_cofactor if f * st > 0])
    cofactor_useage = sum([f * st for (r, st, f) in reactions_using_cofactor if f * st < 0])
    return cofactor_production, cofactor_useage, reactions_using_cofactor, flux


def fix_flux(flux):
    return dict([(k.replace('_reverse', ''), -v) if '_reverse' in k else (k,v)
                 for k, v in flux.iteritems() if abs(v) > 1e-5])
