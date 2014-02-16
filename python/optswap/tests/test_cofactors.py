import cobra
import cobra.io

from optswap import track_cofactor, fix_flux

def test_track_cofactor():
    model = cobra.io.load_matlab_model('/Users/zaking/models/iJO1366_cobrapy.mat')

    # anaerobic
    model.reactions.get_by_id('EX_o2_e').lower_bound = 0

    # run
    production, consumption, reactions, flux = track_cofactor(model, 'EX_cys__L_e', 'nadph_c')
    assert type(float(production)) is float
    assert production > 0
    assert type(float(consumption)) is float
    assert consumption < 0

def test_fix_flux():
    flux = fix_flux({'a':1, 'b':0, 'b_reverse':2, 'a_reverse':0, 'c': 1e-20})
    assert flux['a'] == 1
    assert flux['b'] == -2
    assert 'a_reverse' not in flux
    assert 'b_reverse' not in flux
    assert 'c_reverse' not in flux
    
